/**
  rdmd -Ipath/to/sambamba -O -release -inline flowindices.d B7-295.bam > flowindices.dat
*/
import bamfile;

import fz.flowcall;

import std.stdio;

void main(string[] args) {
    auto bam = BamFile(args[1]);

    uint[1024] counts;

    auto flow_order = bam.header.read_groups["9IKNG"].flow_order;

    foreach (read; bam.alignments) {
        if (read.is_unmapped)
            continue;

        foreach (call; readFlowCalls(read, flow_order)) {
            counts[call.flow_index] += 1; 
        }
    }

    writeln("flow\tcount");

    foreach (size_t i, count; counts) {
        if (count != 0) 
            writeln(i, '\t', count);
    }
}
