module accumulators.flowstats;

import constants;
import bio.core.base;

import std.stdio;
import std.algorithm;
import std.range;
import std.exception;
import std.conv;

class FlowStatsAccumulator
{
    private 
    {
        // nucleotide -> called length -> intensity value -> counts
        // counts:
        //          total
        //          overcalls (approx.)
        //          undercalls (approx.)
        //          mismatches (approx.)
        uint[4][][][4] _distributions;

        uint _max_len;
        uint _max_int;
    }

    this(uint max_length=48, uint max_intensity_value=MAX_INTENSITY_VALUE)
    {
        _max_len = max_length;
        _max_int = max_intensity_value;

        foreach (nuc; 0 .. 4)
        {
            _distributions[nuc] = new uint[4][][](max_length, max_intensity_value);
        }
    }

    static auto merge(FlowStatsAccumulator acc1, FlowStatsAccumulator acc2)
    {
        auto max_length = max(acc1._max_len, acc2._max_len);
        auto max_intensity_value = max(acc1._max_int, acc2._max_int);
        auto acc = new FlowStatsAccumulator(max_length, max_intensity_value);

        foreach (i; 0 .. 4)
            foreach (len; 0 .. max_length)
                foreach (intensity; 0 .. max_intensity_value)
                {
                    uint[4] result;

                    if (len < acc1._max_len && intensity < acc1._max_int)
                       result[] += acc1._distributions[i][len][intensity][];

                    if (len < acc2._max_len && intensity < acc2._max_int)
                       result[] += acc2._distributions[i][len][intensity][];

                    acc._distributions[i][len][intensity][] = result[];
                }

        return acc;
    }

    void updateStatistics(BaseInfo)(BaseInfo[] bases)
    {
        foreach (size_t offset, baseinfo; bases)
        {
            if (offset == 0)
                goto process;

            auto prev_int = bases[offset - 1].flow_call.intensity_value;

            if (prev_int != baseinfo.flow_call.intensity_value)
                goto process;

            auto prev_base = bases[offset - 1].flow_call.base;

            if (prev_base != baseinfo.flow_call.base)
                goto process;

            continue;
process:
            auto flow_call = baseinfo.flow_call;
            auto flow_call_base = cast(Base5)flow_call.base;
            auto len = flow_call.length;
            auto intensity = flow_call.intensity_value;
            increment(flow_call_base.internal_code, len, intensity, 0);
        }
    }

    void printReport(string filename)
    {
        auto _out = File(filename, "w+");

        _out.writeln("# base: called base");
        _out.writeln("# length: called length");
        _out.writeln("# intensity: intensity value");
        _out.writeln("# count: number of flow calls with these base, length, and intensity");
        _out.writeln("# ins: approximate number of overcalls");
        _out.writeln("# del: approximate number of undercalls");
        _out.writeln("# mismatch: approximate number of mismatches");

        _out.writeln("base\tlength\tintensity\tcount\tins\tdel\tmismatch");

        foreach (size_t base, intensity_distributions; _distributions)
            foreach (size_t length, distribution; intensity_distributions)
                foreach (size_t intensity_value, counts; distribution)
                    if (counts[0] != 0 || counts[1] != 0 || counts[2] != 0 || counts[3] != 0) 
                    {
                        _out.writeln(Base5.fromInternalCode(cast(ubyte)base), '\t', 
                                     length, '\t', 
                                     intensity_value, '\t', 
                                     counts[0], '\t',
                                     counts[1], '\t',
                                     counts[2], '\t',
                                     counts[3]);
                    }
    }

    void updateInsertionStats(I)(I insertion)
    {
        foreach (flow_call; insertion.flowcalls)
        {
            auto flow_call_base = cast(Base5)flow_call.base;
            auto len = flow_call.length;
            auto intensity = flow_call.intensity_value;
            increment(flow_call_base.internal_code, len, intensity, 1);
        }
    }

    void updateDeletionStats(D)(D deletion)
    {
        auto intensities = deletion.deleted_base_intensities.data;
        foreach (base_and_intensity; zip(uniq(deletion.bases), uniq(intensities)))
        {
            auto base5 = Base5(base_and_intensity[0]);
            auto intensity = base_and_intensity[1];
            auto len = (intensity + 50) / 100; // FIXME
            increment(base5.internal_code, len, intensity, 2);
        }
    }

    void updateMismatchStats(M)(M mismatch)
    {
        with (mismatch)
        {
            auto intensity = base.flow_call.intensity_value;
            auto len = base.flow_call.length;
            auto base5 = cast(Base5)base.base;
            auto code = base5.internal_code;
            increment(code, len, intensity, 3);
        }
    }

    private void increment(size_t internal_code, size_t len, int intensity, int index)
    {
       // FIXME
       if (intensity < 0) intensity = 0;

       enforce(len < _distributions[0].length, 
               "Unexpectedly long homopolymer (length " ~ to!string(len) ~ ")");

       enforce(len < _distributions[0][0].length, 
               "Unexpectedly large intensity value (" ~ to!string(intensity) ~ ")");

       _distributions[internal_code][len][intensity][index] += 1;
    }
}
