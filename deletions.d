/**
    rdmd -Ipath/to/sambamba -O -release -inline deletions.d B7-295.bam > deletions.dat
*/
import bamfile;

import fz.flowcall;

import BioD.Fasta;
import BioD.Base;

import std.algorithm;
import std.range;
import std.stdio;

__gshared string reference;

void searchDeletions(Cigar, F, R)(Cigar cigar, F flow_calls, R read) {
    auto offset = 0;
    auto ref_offset = 0;
    ElementType!F prev_call;

    // loop invariant:
    // flow_calls.front.offset >= offset

    char[] deleted_sequence;

    foreach (op; cigar) {
        if (op == 'D') {
            if (offset == 0) continue;

            auto len = op.length;
            auto curr_call = flow_calls.empty ? prev_call : flow_calls.front;
           
            auto strand = read.is_reverse_strand ? '-' : '+';
            auto pos = read.is_reverse_strand ? (read.position + read.basesCovered() - ref_offset - 1) 
                                              : (read.position + ref_offset);

            deleted_sequence.length = len;

            Base prev_base = prev_call.base;
            Base curr_base = curr_call.base;

            if (read.is_reverse_strand) {
                foreach_reverse(size_t i, c; reference[(pos - len + 1) .. pos + 1])
                    deleted_sequence[i] = Base(c).complement;
            } else {
                foreach(size_t i, c; reference[pos .. pos + len])
                    deleted_sequence[i] = Base(c);
            }

            writeln(pos, '\t', strand, '\t', offset, '\t', len, '\t',
                    prev_call.offset, '\t', prev_base, '\t', 
                    prev_call.length, '\t', prev_call.intensity, '\t',
                    curr_call.offset, '\t', curr_base, '\t',
                    curr_call.length, '\t', curr_call.intensity, '\t',
                    deleted_sequence, '\t', prev_call.flow_index, '\t',
                    curr_call.flow_index);
        }

        if (op.is_query_consuming) {
            offset += op.length;
            while (!flow_calls.empty && flow_calls.front.offset < offset) {
                prev_call = flow_calls.front;
                flow_calls.popFront();
            }
        }

        if (op.is_reference_consuming) {
            ref_offset += op.length;
        }
    }
}

void process(R)(R read, string flow_order) {
    auto flow_calls = readFlowCalls(read, flow_order);
    if (read.is_reverse_strand) {
        searchDeletions(retro(read.cigar), flow_calls, read);
    } else {
        searchDeletions(read.cigar, flow_calls, read);
    }
}

void main(string[] args) {
    auto bam = BamFile(args[1]);

    auto flow_order = bam.header.read_groups["9IKNG"].flow_order;
    auto reads = filter!"!a.is_unmapped"(bam.alignments);

    reference = fastaRecords("U00096.2.fas").front.sequence;


    write("pos\tstrand\toffset\tlen\tprev.offset\tprev.nuc\tprev.len\tprev.int");
    write("\tnext.offset\tnext.nuc\tnext.len\tnext.int\tseq");
    write("\tprev.flow\tnext.flow");
    writeln();

    foreach (read; reads) {
        process(read, flow_order);
    }
}
