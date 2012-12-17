module printers.deletioninfo;

import std.stdio;
import std.algorithm;
import std.conv;

class DeletionInfoPrinter
{
    private 
    {
        File _out;
        string _flow_order;
    }

    this(string filename, string flow_order)
    {
        _out = File(filename, "w+");
        _flow_order = flow_order;

        writeComments();
        writeHeader();
    }

    private void writeComments()
    {
        // FIXME: clarify meaning of these numbers
        _out.writeln("# ref.pos: 1-based position on the reference");
        _out.writeln("# offset: 0-based offset from the read beginning");
        
        _out.writeln("# prev: base called before the deletion");
        _out.writeln("# prev.int: flow signal intensity of the previous base");
        _out.writeln("# prev.len: how many bases were called from the flowcall");
        _out.writeln("#           the previous base was called from");

        _out.writeln("# len: deletion length");
        _out.writeln("# bases: deleted sequence");

        _out.writeln("# next: base called after the deletion");
        _out.writeln("# next.int: flow signal intensity of the next base");
        _out.writeln("# next.len: how many bases were called from the flowcall");
        _out.writeln("#           the next base was called from");

        _out.writeln("# hom: number of homopolymers in the deleted sequence");
        _out.writeln("# ints: comma-delimited list of intensities for these homopolymers");
    }

    private void writeHeader()
    {
        _out.writeln("ref.pos\toffset\tprev\tprev.int\tprev.len\tlen\tbases\t" ~
                     "next\tnext.int\tnext.len\thom\tints");
    }

    void printDeletion(Deletion)(/*const*/ ref Deletion deletion, ushort[] intensities)
    {
        with (deletion)
        {
            _out.write(start_position, '\t', read_offset, '\t',
                       previous_flow_call.base, '\t', 
                       previous_flow_call.intensity_value, '\t',
                       previous_flow_call.length, '\t',
                       length, '\t',
                       bases, '\t',
                       next_flow_call.base, '\t', 
                       next_flow_call.intensity_value, '\t', 
                       next_flow_call.length, '\t',
                       number_of_homopolymers_in_deleted_sequence, '\t');

            auto unique_bases = bases.uniq();
            assert(!unique_bases.empty);
            auto current_base = unique_bases.front;
            size_t j = 0;

            ushort[1024] buf;
            ushort[] deleted_base_intensities = void;
            auto n = number_of_homopolymers_in_deleted_sequence;
            if (n > 1024) // unlikely to happen
                deleted_base_intensities = new ushort[n];
            else
                deleted_base_intensities = buf[0 .. n];

            for (size_t i = previous_flow_call.flow_index;
                        i <= next_flow_call.flow_index;
                        i++)
            {
                if (_flow_order[i] == current_base)
                {
                    deleted_base_intensities[j++] = intensities[i];
                    if (j == n)
                    {
                        break;
                    }
                    unique_bases.popFront();
                    if (!unique_bases.empty)
                    {
                        current_base = unique_bases.front;
                    }
                    else
                    {
                        break;
                    }
                }
            }

            if (j == n)
            {
                _out.write(deleted_base_intensities.map!(to!string)().joiner(","));
            }
            else
            {
                _out.write("NA");
            }

            _out.writeln();
        }
    }
}
