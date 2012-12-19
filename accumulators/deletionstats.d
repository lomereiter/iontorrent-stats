module accumulators.deletionstats;

import std.stdio;

class DeletionStatsAccumulator
{
    private
    {
        size_t _n_deletions;
        size_t _n_undercalls;

        uint[] _intensities; // undercall intensities
    }

    this(size_t maximum_intensity_value=1536)
    {
        _intensities = new uint[maximum_intensity_value];
    }

    void updateStatistics(Deletion)(Deletion deletion)
    {
        ++_n_deletions;

        with (deletion)
        {
            if (number_of_homopolymers_in_deleted_sequence != deleted_base_intensities.length)
                return;

            _n_undercalls += number_of_homopolymers_in_deleted_sequence;

            foreach (intensity; deleted_base_intensities)
                _intensities[intensity] += 1;
        }
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
