module printers.deletioninfo;

import std.stdio;
import std.algorithm;
import std.conv;

class DeletionInfoPrinter
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
        _out.writeln("# ref.pos: 1-based position of the first deleted base on the reference");
        _out.writeln("# offset: 0-based offset of the read base called before the deletion");
        
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

    void printDeletion(Deletion)(/*const*/ ref Deletion deletion)
    {
        with (deletion)
        {
            _out.write(start_position + 1, '\t', read_offset, '\t',
                       previous_flow_call.base, '\t', 
                       previous_flow_call.intensity_value, '\t',
                       previous_flow_call.length, '\t',
                       length, '\t',
                       bases, '\t',
                       next_flow_call.base, '\t', 
                       next_flow_call.intensity_value, '\t', 
                       next_flow_call.length, '\t',
                       number_of_homopolymers_in_deleted_sequence, '\t');

            if (deleted_base_intensities.length == number_of_homopolymers_in_deleted_sequence)
            {
                _out.write(deleted_base_intensities.data.map!(to!string)().joiner(","));
            }
            else
            {
                _out.write("NA");
            }

            _out.writeln();
        }
    }
}
