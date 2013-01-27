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
        uint[] _swap_start_intensities;
        uint[] _swap_end_intensities;

        size_t _n_swaps;
    }

    this(uint max_intensity_value=MAX_INTENSITY_VALUE)
    {
        _intensities = new uint[max_intensity_value];
        _swap_start_intensities = new uint[max_intensity_value];
        _swap_end_intensities = new uint[max_intensity_value];
    }

    static auto merge(MismatchStatsAccumulator acc1, MismatchStatsAccumulator acc2)
    {
        auto acc = new MismatchStatsAccumulator();

        acc._n_mismatches = acc1._n_mismatches + acc2._n_mismatches;

        static auto mergeArrays(in uint[] a, in uint[] b)
        {
            auto c = new uint[max(a.length, b.length)];
            c[0 .. a.length] = a[];
            c[0 .. b.length] += b[];
            return c;
        }

        acc._intensities = mergeArrays(acc1._intensities, acc2._intensities);
        acc._swap_start_intensities = mergeArrays(acc1._swap_start_intensities, 
                                                  acc2._swap_start_intensities);
        acc._swap_end_intensities = mergeArrays(acc1._swap_end_intensities, 
                                                acc2._swap_end_intensities);

        acc._n_swaps = acc1._n_swaps + acc2._n_swaps;

        return acc;
    }

    void updateStatistics(Mismatch)(Mismatch mismatch)
    {
        ++_n_mismatches;

        with (mismatch)
        {
            auto int_v = base_info.flow_call.intensity_value;

            // FIXME
            if (int_v > 0) 
                _intensities[int_v] += 1; 

            if (starts_swap)
            {
                ++_n_swaps;
                _swap_start_intensities[int_v] += 1;
            }

            if (ends_swap)
            {
                _swap_end_intensities[int_v] += 1;
            }
        }
    }

    void printMismatchesReport(string filename)
    {
        auto _file = File(filename, "w+");

        // TODO: description

        _file.writeln("intensity\tcount\tswap.start\tswap.end");

        foreach (size_t intensity, count; _intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count,
                                         '\t', _swap_start_intensities[intensity],
                                         '\t', _swap_end_intensities[intensity]);
        }
    }
}
