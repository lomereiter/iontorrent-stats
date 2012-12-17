module printers.insertioninfo;

import bio.bam.fz.flowcall;
import std.typecons;

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
        // FIXME: clarify meaning of these numbers
        _out.writeln("# ref.pos: 1-based position on the reference");
        _out.writeln("# offset: 0-based offset from the read beginning");

        _out.writeln("# prev: base called before the insertion");
        _out.writeln("# prev.int: flow signal intensity of the previous base");
        _out.writeln("# prev.len: how many bases were called from the flowcall");
        _out.writeln("#           the previous base was called from");

        _out.writeln("# len: insertion length");
        _out.writeln("# bases: inserted sequence");

        _out.writeln("# next: base called after the insertion");
        _out.writeln("# next.int: flow signal intensity of the next base");
        _out.writeln("# next.len: how many bases were called from the flowcall");
        _out.writeln("#           the next base was called from");

        _out.writeln("# hom: number of homopolymers in the inserted sequence");
        _out.writeln("# ints: comma-delimited list of flow call intensities for these homopolymers");
    }

    private void writeHeader()
    {
        _out.writeln("ref.pos\toffset\tprev\tprev.int\tprev.len\tlen\tbases\t" ~
                     "next\tnext.int\tnext.len\thom\tints");
    }
   
    private static void printNullableFlowcall(Nullable!ReadFlowCall flowcall, ref File file)
    {
        if (flowcall.isNull)
            file.write("NA\tNA\tNA\t");
        else
            file.write(flowcall.base, '\t', 
                       flowcall.intensity_value, '\t', 
                       flowcall.length, '\t');
    }

    void printInsertion(Insertion)(const ref Insertion insertion)
    {
        with (insertion)
        {
            _out.write(reference_position, '\t', start_offset, '\t');
            printNullableFlowcall(previous_flow_call, _out);
            _out.write(length, '\t', bases, '\t');
            printNullableFlowcall(next_flow_call, _out);
            _out.write(number_of_flowcalls, '\t');
            _out.writeln(flowcalls.map!"a.intensity_value"().map!(to!string)().joiner(","));
        }
    }
}
