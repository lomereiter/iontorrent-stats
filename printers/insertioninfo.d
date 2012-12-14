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

        uint[] _neighbour_intensities; // inserted bases are from a single flow call
        ushort[][] _intensity_pairs;       // inserted bases are from two flow calls
        uint _total_insertions;
        uint _insertions_of_type_1; // homopolymer-related of length 1
        uint _insertions_of_type_2; // homopolymer-related of length 2
        uint _insertions_of_type_3; // homopolymer-related of length 3
        uint[] _orphan_intensities;     // calls of unexistent bases

        uint[] _intensities; // unifies all three types of errors
    }

    this(string filename, uint max_intensity_value=1536)
    {
        _out = File(filename, "w+");

        _neighbour_intensities = new uint[max_intensity_value];
        _orphan_intensities = new uint[max_intensity_value];
        _intensities = new uint[max_intensity_value];
        _intensity_pairs = new ushort[][](max_intensity_value, max_intensity_value);

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
                ++_total_insertions;

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

                bool type3 = true;

                if (start_offset > 0)
                {
                    auto previous_baseinfo = bases[start_offset - 1];
                    auto previous_flow_call = previous_baseinfo.flow_call;
                    previous_base = to!string(previous_baseinfo);
                    previous_call_intensity = to!string(previous_flow_call.intensity_value);
                    previous_call_length = to!string(previous_flow_call.length);

                    // check if insertion is related to a single homopolymer
                    if (number_of_flowcalls == 1 && previous_flow_call == flowcalls.front)
                    {
                        _neighbour_intensities[previous_flow_call.intensity_value] += 1;
                        _intensities[previous_flow_call.intensity_value] += 1;
                        _insertions_of_type_1 += 1;
                        type3 = false;
                    }
                }

                if (end_offset < bases.length)
                {
                    auto next_baseinfo = bases[end_offset];
                    auto next_flow_call = next_baseinfo.flow_call;
                    next_base = to!string(next_baseinfo);
                    next_call_intensity = to!string(next_flow_call.intensity_value);
                    next_call_length = to!string(next_flow_call.length);

                    // ditto
                    if (number_of_flowcalls == 1 && next_flow_call == flowcalls.front)
                    {
                        _neighbour_intensities[next_flow_call.intensity_value] += 1;
                        _intensities[next_flow_call.intensity_value] += 1;
                        _insertions_of_type_1 += 1;
                        type3 = false;
                    }
                }

                // check for type 2 (call of unexistent base)
                if (type3 && number_of_flowcalls == 1)
                {
                    _insertions_of_type_2 += 1;
                    _orphan_intensities[flowcalls.front.intensity_value] += 1;
                    _intensities[flowcalls.front.intensity_value] += 1;
                }

                // finally, check the case of two homopolymers glued together (type 3)
                if (number_of_flowcalls == 2 &&
                    start_offset > 0 && end_offset < bases.length)
                {
                    auto left = bases[start_offset].flow_call;
                    auto right = bases[end_offset - 1].flow_call;

                    if (left == bases[start_offset - 1].flow_call &&
                        right == bases[end_offset].flow_call)
                    {
                        _intensity_pairs[left.intensity_value][right.intensity_value] += 1;
                        _insertions_of_type_3 += 1;

                        _intensities[left.intensity_value] += 1;
                        _intensities[right.intensity_value] += 1;
                    }
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

    void printNeighbourSummary(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# total number of insertions: ", _total_insertions);

        _file.writeln("# number of one-sided overcalls: ",
                      _insertions_of_type_1, 
                      " (",
                      cast(double)_insertions_of_type_1 * 100.0 / _total_insertions,
                      "%)");

        _file.writeln("# number of 'orphans' (single flowcalls different from both neighbors): ",
                      _insertions_of_type_2, 
                      " (",
                      cast(double)_insertions_of_type_2 * 100.0 / _total_insertions,
                      "%)");

        _file.writeln("# number of two-sided overcalls: ",
                      _insertions_of_type_3,
                      " (",
                      cast(double)_insertions_of_type_3 * 100.0 / _total_insertions,
                      "%)");
    }

    void printOrphansReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# There are a lot of insertions that consist of a single flow call");
        _file.writeln("# which is different from both neighbour flow calls. As it turns out");
        _file.writeln("# most of these insertions have signal intensity close to 50.");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _orphan_intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }

    void printOneSidedOvercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# Typically, most insertion errors are due to homopolymer overcall");
        _file.writeln("# Below is the signal intensity distribution of these flow calls");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _neighbour_intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }

    void printTwoSidedOvercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# Typically, aligner merges two neighboring overcalls into one insertion");
        _file.writeln("# because one insertion is 'cheaper' than two, in its view.");
        _file.writeln("# ");
        _file.writeln("# first: intensity of the first flow call");
        _file.writeln("# second: intensity of the second flow call");

        _file.writeln("first\tsecond");

        foreach (size_t first, dist; _intensity_pairs)
        {
            foreach (size_t second, count; dist)
                for (ushort i = 0; i < count; ++i)
                    _file.writeln(first, '\t', second);
        }
    }

    void printOvercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# intensity: round(normalized signal intensity * 100.0)");
        _file.writeln("# count: number of overcalls with this intensity");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }
}
