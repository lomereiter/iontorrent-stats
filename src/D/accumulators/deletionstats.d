module accumulators.deletionstats;

import constants;

import std.stdio;
import std.exception;
import std.conv : to;
import std.algorithm : max;

class DeletionStatsAccumulator
{
    private
    {
        size_t _n_deletions;
        size_t _n_resolved; // likely to be caused solely by undercalls

        size_t _n_undercalls;

        uint[] _intensities; // undercall intensities
    }

    this(size_t maximum_intensity_value=MAX_INTENSITY_VALUE)
    {
        _intensities = new uint[maximum_intensity_value];
    }

    static auto merge(DeletionStatsAccumulator acc1,
                      DeletionStatsAccumulator acc2)
    {
        auto max_intensity_value = max(acc1._intensities.length,
                                       acc2._intensities.length);

        auto acc = new DeletionStatsAccumulator(max_intensity_value);
        acc._n_deletions = acc1._n_deletions + acc2._n_deletions;
        acc._n_resolved = acc1._n_resolved + acc2._n_resolved;
        acc._n_undercalls = acc1._n_undercalls + acc2._n_undercalls;
        acc._intensities[0 .. acc1._intensities.length] = acc1._intensities[];
        acc._intensities[0 .. acc2._intensities.length] += acc2._intensities[];

        return acc;
    }

    void updateStatistics(Deletion)(Deletion deletion)
    {
        ++_n_deletions;

        with (deletion)
        {
            if (number_of_homopolymers_in_deleted_sequence != deleted_base_intensities.length)
                return;

            _n_undercalls += number_of_homopolymers_in_deleted_sequence;
            ++_n_resolved;

            foreach (intensity; deleted_base_intensities)
            {
                if (intensity < 0) intensity = 0; // FIXME

                enforce(intensity < _intensities.length,
                        "Unexpectedly large intensity value (" ~ to!string(intensity) ~ ")");

                _intensities[intensity] += 1;
            }
        }
    }

    void printSummary(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("######## DELETIONS #######");
        _file.writeln("# Number of deletions: ", _n_deletions);
        _file.writeln("# Likely to be caused solely by undercalls: ", _n_resolved,
                      " (", cast(double)_n_resolved * 100.0 / _n_deletions, "%)");

        _file.writeln("# Number of undercalls: ", _n_undercalls);
    }

    void printUndercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# intensity: signal intensity value for undercalled homopolymer");
        _file.writeln("# count: number of undercalls with this intensity value");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }
}
