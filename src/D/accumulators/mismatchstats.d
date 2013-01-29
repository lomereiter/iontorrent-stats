module accumulators.mismatchstats;

import constants;

import std.stdio;
import std.algorithm : max;

class MismatchStatsAccumulator
{
    private 
    {
        size_t _n_mismatches;
        uint[] _intensities;
        uint[] _swap_start_intensities;
        uint[] _swap_end_intensities;
        uint[] _uo_under_intensities;
        uint[] _uo_over_intensities;
        uint[] _ou_over_intensities;
        uint[] _ou_under_intensities;
        uint[] _sandwich_intensities;

        size_t _n_swaps;
        size_t _n_under_over;
        size_t _n_over_under;
        size_t _n_sandwich;
    }

    this(uint max_intensity_value=MAX_INTENSITY_VALUE)
    {
        _intensities = new uint[max_intensity_value];
        _swap_start_intensities = new uint[max_intensity_value];
        _swap_end_intensities = new uint[max_intensity_value];
        _uo_under_intensities = new uint[max_intensity_value];
        _uo_over_intensities = new uint[max_intensity_value];
        _ou_over_intensities = new uint[max_intensity_value];
        _ou_under_intensities = new uint[max_intensity_value];
        _sandwich_intensities = new uint[max_intensity_value];
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

        acc._uo_under_intensities = mergeArrays(acc1._uo_under_intensities, 
                                                acc2._uo_under_intensities);
        acc._uo_over_intensities = mergeArrays(acc1._uo_over_intensities, 
                                               acc2._uo_over_intensities);
        acc._ou_over_intensities = mergeArrays(acc1._ou_over_intensities, 
                                               acc2._ou_over_intensities);
        acc._ou_under_intensities = mergeArrays(acc1._ou_under_intensities, 
                                                acc2._ou_under_intensities);

        acc._sandwich_intensities = mergeArrays(acc1._sandwich_intensities,
                                                acc2._sandwich_intensities);

        acc._n_swaps = acc1._n_swaps + acc2._n_swaps;
        acc._n_under_over = acc1._n_under_over + acc2._n_under_over;
        acc._n_over_under = acc1._n_over_under + acc2._n_over_under;
        acc._n_sandwich = acc1._n_sandwich + acc2._n_sandwich;

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
            {
                _intensities[int_v] += 1; 

                if (starts_swap)
                {
                    auto iv2 = next_base_info.flow_call.intensity_value;

                    if (iv2 > 0)
                    {
                        ++_n_swaps;
                        _swap_start_intensities[int_v] += 1;
                        _swap_end_intensities[iv2] += 1;
                    }
                }
                else if (is_under_over)
                {
                    auto iv1 = previous_base_info.flow_call.intensity_value;
                    auto iv2 = base_info.flow_call.intensity_value;
                    if (iv1 > 0 && iv2 > 0)
                    {
                        ++_n_under_over;
                        _uo_under_intensities[iv1] += 1;
                        _uo_over_intensities[iv2] += 1;
                    }
                }
                else if (is_over_under)
                {
                    auto iv1 = base_info.flow_call.intensity_value;
                    auto iv2 = next_base_info.flow_call.intensity_value;
                    if (iv1 > 0 && iv2 > 0)
                    {
                        ++_n_over_under;
                        _ou_over_intensities[iv1] += 1;
                        _ou_under_intensities[iv2] += 1;
                    }
                } else if (is_between_matches)
                {
                    auto iv = base_info.flow_call.intensity_value;
                    if (iv > 0)
                    {
                        ++_n_sandwich;
                        _sandwich_intensities[iv] += 1;
                    }
                }
            }
        }
    }

    void printSummary(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("######## MISMATCHES #######");
        _file.writeln("# Number of mismatches: ", _n_mismatches);
        _file.writeln("# Swaps (XY -> YX): ", _n_swaps, 
                      " (", cast(double)_n_swaps * 200.0 / _n_mismatches, "%)");
        _file.writeln("# Undercall + overcall (XYY -> XXY): ", _n_under_over,
                      " (", cast(double)_n_under_over * 100.0 / _n_mismatches, "%)");
        _file.writeln("# Overcall + undercall (XXY -> XYY): ", _n_over_under,
                      " (", cast(double)_n_over_under * 100.0 / _n_mismatches, "%)");
        _file.writeln("# Between two matched bases (XYZ -> XWZ), excluding above two cases: ", _n_sandwich,
                      " (", cast(double)_n_sandwich * 100.0 / _n_mismatches, "%)");
    }

    void printMismatchesReport(string filename)
    {
        auto _file = File(filename, "w+");

        // TODO: description

        _file.writeln("intensity\tcount\tswap.start\tswap.end\tuo.under\tuo.over\tou.over\tou.under\tsandwich");

        foreach (size_t intensity, count; _intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count,
                                         '\t', _swap_start_intensities[intensity],
                                         '\t', _swap_end_intensities[intensity],
                                         '\t', _uo_under_intensities[intensity],
                                         '\t', _uo_over_intensities[intensity],
                                         '\t', _ou_over_intensities[intensity],
                                         '\t', _ou_under_intensities[intensity],
                                         '\t', _sandwich_intensities[intensity]);
        }
    }
}
