module accumulators.flowstats;

import bio.core.base;

import std.stdio;
import std.algorithm;
import std.range;

class FlowStatsAccumulator
{
    private 
    {
        // nucleotide -> called length -> intensity value -> counts
        // counts:
        //          total
        //          overcalls (approx.)
        //          undercalls (approx.)
        uint[3][][][4] _distributions;

        uint _max_len;
        uint _max_int;
    }

    this(uint max_length=16, uint max_intensity_value=1536)
    {
        _max_len = max_length;
        _max_int = max_intensity_value;

        foreach (nuc; 0 .. 4)
        {
            _distributions[nuc] = new uint[3][][](max_length, max_intensity_value);
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
                    uint[3] result;

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
            _distributions[flow_call_base.internal_code][len][intensity][0] += 1;
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

        _out.writeln("base\tlength\tintensity\tcount\tins\tdel");

        foreach (size_t base, intensity_distributions; _distributions)
            foreach (size_t length, distribution; intensity_distributions)
                foreach (size_t intensity_value, counts; distribution)
                    if (counts[0] != 0 || counts[1] != 0 || counts[2] != 0) 
                    {
                        _out.writeln(Base5.fromInternalCode(cast(ubyte)base), '\t', 
                                     length, '\t', 
                                     intensity_value, '\t', 
                                     counts[0], '\t',
                                     counts[1], '\t',
                                     counts[2]);
                    }
    }

    void updateInsertionStats(I)(I insertion)
    {
        foreach (flow_call; insertion.flowcalls)
        {
            auto flow_call_base = cast(Base5)flow_call.base;
            auto len = flow_call.length;
            auto intensity = flow_call.intensity_value;
            _distributions[flow_call_base.internal_code][len][intensity][1] += 1;
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
            _distributions[base5.internal_code][len][intensity][2] += 1;
        }
    }
}
