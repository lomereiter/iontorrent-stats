import bamfile;
import fz.flowcall;

import std.stdio;

size_t asInt(char c) {
    switch (c) { 
        case 'A': return 0; 
        case 'C': return 1; 
        case 'G': return 2; 
        case 'T': return 3; 
        default: assert(0); 
    }
    assert(0);
}

char asChar(size_t i) {
    return "ACGT"[i];
}


void main(string[] args) {

    // nucleotide -> called length -> intensity value -> count
    ulong[2048][16][4] distributions;

    auto bam = BamFile(args[1]);
    auto flow_order = bam.header.read_groups["9IKNG"].flow_order;

    foreach (read; bam.alignments) {
        if (read.is_unmapped)
            continue;

        foreach (call; readFlowCalls(read, flow_order)) {
            distributions[asInt(call.base)][call.length][call.intensity_value] += 1;
        }
    }

    writeln("base\tlength\tintensity\tcount");

    foreach (size_t base, intensity_distributions; distributions) {
        foreach (size_t length, distribution; intensity_distributions) {
            foreach (size_t intensity_value, count; distribution) {
                if (count != 0) {
                    writeln(asChar(base), '\t', length, '\t', intensity_value, '\t', count, '\t');
                }
            }
        }
    }
}
