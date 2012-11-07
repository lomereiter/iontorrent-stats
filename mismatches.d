/**
  rdmd -Ipath/to/sambamba -O -release -inline mismatches.d B7-295.bam > mismatches.dat
*/
import bamfile;

import baseinfo;

import BioD.Base;
import BioD.Fasta;

import std.stdio;
import std.ascii;
import std.parallelism;
import std.algorithm;
import std.conv;

class MismatchFinder {
    this(string filename, string fastafilename, string read_group, TaskPool taskpool) {

        bam = BamFile(filename, taskpool);
        reference = to!string(map!toUpper(fastaRecords(fastafilename).front.sequence));
        flow_order = bam.header.read_groups[read_group].flow_order;
    }

    void process() {

        foreach (read; bam.alignments) {
            if (read.is_unmapped)
                continue;

            size_t i = 0;
            foreach (base; basesWith!"FZ"(read, arg!"FZ"(flow_order))) {
                if (base.cigar_operation == 'M' && reference[base.position] != base)
                    writeln(base.position, '\t', read.strand, '\t', 
                            reference[base.position], '\t', to!string(base), '\t', 
                            base.flow_call.intensity, '\t', i);

                ++i;
            }
        }
    }

    private {
        BamFile bam;
        string reference;
        string flow_order;
    }
}

void main(string[] args) {
    auto taskpool = new TaskPool(totalCPUs);
    scope(exit) taskpool.finish();
    auto mismatchFinder = new MismatchFinder(args[1], "U00096.2.fas", "9IKNG", taskpool);
    writeln("pos\tstrand\tref.nuc\tnuc\tintensity\toffset");
    mismatchFinder.process();
}
