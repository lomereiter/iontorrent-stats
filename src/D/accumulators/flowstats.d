module accumulators.flowstats;

import bio.core.base;

import std.stdio;
import std.algorithm : max;

class FlowStatsAccumulator
{
    private 
    {
        // nucleotide -> called length -> intensity value -> count
        uint[][][4] _distributions;

        uint _max_len;
        uint _max_int;
    }

    this(uint max_length=16, uint max_intensity_value=1536)
    {
        _max_len = max_length;
        _max_int = max_intensity_value;

        foreach (nuc; 0 .. 4)
        {
            _distributions[nuc] = new uint[][](max_length, max_intensity_value);
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
                    uint result = 0;

                    if (len < acc1._max_len && intensity < acc1._max_int)
                       result += acc1._distributions[i][len][intensity];

                    if (len < acc2._max_len && intensity < acc2._max_int)
                       result += acc2._distributions[i][len][intensity];

                    acc._distributions[i][len][intensity] = result;
                }

        return acc;
    }

    void updateStatistics(BaseInfo)(BaseInfo[] bases)
    {
        foreach (size_t offset, baseinfo; bases)
        {
            if (offset == 0 || bases[offset - 1].flow_call != baseinfo.flow_call)
            {
                auto flow_call = baseinfo.flow_call;
                auto flow_call_base = cast(Base5)flow_call.base;
                auto len = flow_call.length;
                auto intensity = flow_call.intensity_value;
                _distributions[flow_call_base.internal_code][len][intensity] += 1;
            }
        }
    }

    void printReport(string filename)
    {
        auto _out = File(filename, "w+");

        _out.writeln("# base: called base");
        _out.writeln("# length: called length");
        _out.writeln("# intensity: intensity value");
        _out.writeln("# count: number of flow calls with these base, length, and intensity");

        _out.writeln("base\tlength\tintensity\tcount");

        foreach (size_t base, intensity_distributions; _distributions)
            foreach (size_t length, distribution; intensity_distributions)
                foreach (size_t intensity_value, count; distribution)
                    if (count != 0) 
                    {
                        _out.writeln(Base5.fromInternalCode(cast(ubyte)base), '\t', 
                                     length, '\t', 
                                     intensity_value, '\t', 
                                     count);
                    }
    }
}
