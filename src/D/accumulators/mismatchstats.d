module accumulators.mismatchstats;

import constants;

import std.stdio;
import std.algorithm : max;

// TODO: classify mismatches into types and collect per-type statistics

class MismatchStatsAccumulator
{
    private 
    {
        size_t _n_mismatches;
        uint[] _intensities;
    }

    this(uint max_intensity_value=MAX_INTENSITY_VALUE)
    {
        _intensities = new uint[max_intensity_value];
    }

    static auto merge(MismatchStatsAccumulator acc1, MismatchStatsAccumulator acc2)
    {
        auto acc = new MismatchStatsAccumulator();

        acc._n_mismatches = acc1._n_mismatches + acc2._n_mismatches;
        acc._intensities.length = max(acc1._intensities.length, acc2._intensities.length);
        acc._intensities[0 .. acc1._intensities.length] = acc1._intensities[];
        acc._intensities[0 .. acc2._intensities.length] += acc2._intensities[];

        return acc;
    }

    void updateStatistics(Mismatch)(Mismatch mismatch)
    {
        ++_n_mismatches;

        with (mismatch)
        {
            auto int_v = base.flow_call.intensity_value;

            // FIXME
            if (int_v > 0) 
                _intensities[int_v] += 1; 
        }
    }

    void printMismatchesReport(string filename)
    {
        auto _file = File(filename, "w+");

        // TODO: description

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }
}
