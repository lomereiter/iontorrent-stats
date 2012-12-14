module printers.insertioninfo;

import bio.bam.fz.flowcall;

import std.stdio;
import std.range;
import std.algorithm;
import std.conv;

class InsertionInfoPrinter
{
    private
    {
        File _out;
    }

    this(string filename)
    {
        _out = File(filename, "w+");

        writeComments();
        writeHeader();
    }

    private void writeComments()
    {
        _out.writeln("# ref.pos: 1-based position on the reference");
        _out.writeln("# offset: 0-based offset from the read beginning");

        _out.writeln("# prev: base called before the insertion");
        _out.writeln("# prev.int: flow signal intensity of the previous base");
        _out.writeln("# prev.len: how many bases were called from the flowcall");
        _out.writeln("#           the previous base was called from");

        _out.writeln("# len: insertion length");
        _out.writeln("# sequence: inserted sequence");

        _out.writeln("# next: base called after the insertion");
        _out.writeln("# next.int: flow signal intensity of the next base");
        _out.writeln("# next.len: how many bases were called from the flowcall");
        _out.writeln("#           the next base was called from");

        _out.writeln("# fc.num: number of flow calls in the inserted sequence");
        _out.writeln("# intensities: comma-delimited list of flow call intensities");
    }

    private void writeHeader()
    {
        _out.writeln("ref.pos\toffset\tprev\tprev.int\tprev.len\tlen\tsequence\t" ~
                     "next\tnext.int\tnext.len\tfc.num\tintensities");
    }

    void printInsertions(BaseInfo)(BaseInfo[] bases)
    {
        static auto getFlowCall(const ref BaseInfo info)
        {
            return info.flow_call;
        }

        static auto getIntensityValueStr(ReadFlowCall fc)
        {
            return to!string(fc.intensity_value);
        }

        foreach (size_t offset, baseinfo; bases)
        {
            if (baseinfo.cigar_operation.type == 'I' &&
                baseinfo.cigar_operation_offset == 0)
            {
                auto start_offset = offset;
                auto end_offset = offset + baseinfo.cigar_operation.length;

                assert(end_offset <= bases.length);

                auto insertion_length = end_offset - start_offset;
                auto inserted_sequence = bases[start_offset .. end_offset];

                auto flowcalls = inserted_sequence.map!getFlowCall().uniq();

                auto number_of_flowcalls = flowcalls.walkLength();

                auto reference_position = baseinfo.position + 1;

                auto previous_base = "NA";
                auto previous_call_intensity = "NA";
                auto previous_call_length = "NA";
                auto next_base = "NA";
                auto next_call_intensity = "NA";
                auto next_call_length = "NA";

                if (start_offset > 0)
                {
                    auto previous_baseinfo = bases[start_offset - 1];
                    auto previous_flow_call = previous_baseinfo.flow_call;
                    previous_base = to!string(previous_baseinfo);
                    previous_call_intensity = to!string(previous_flow_call.intensity_value);
                    previous_call_length = to!string(previous_flow_call.length);
                }

                if (end_offset < bases.length)
                {
                    auto next_baseinfo = bases[end_offset];
                    auto next_flow_call = next_baseinfo.flow_call;
                    next_base = to!string(next_baseinfo);
                    next_call_intensity = to!string(next_flow_call.intensity_value);
                    next_call_length = to!string(next_flow_call.length);
                }

                _out.writeln(reference_position, '\t',
                             offset, '\t',
                             previous_base, '\t',
                             previous_call_intensity, '\t',
                             previous_call_length, '\t',
                             insertion_length, '\t',
                             inserted_sequence, '\t',
                             next_base, '\t',
                             next_call_intensity, '\t',
                             next_call_length, '\t',
                             number_of_flowcalls, '\t',
                             flowcalls.map!getIntensityValueStr().joiner(","));
            }
        }
    }
}
