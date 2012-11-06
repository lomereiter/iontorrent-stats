/**
  rdmd -Ipath/to/sambamba -O -release -inline counts.d B7-295.bam > counts.dat
*/
import bamfile;

import BioD.Fasta;

import std.stdio;
import std.ascii;
import std.parallelism;
import std.algorithm;
import std.conv;

ulong[1024] deletions;
ulong[1024] insertions;
ulong[1024] mismatches;
ulong[1024] reads;

void update(alias a)(size_t offset, bool reverse, size_t len) {
    if (reverse)
        a[len - 1 - offset] += 1;
    else
        a[offset] += 1;
}

void main(string[] args) {
    auto reference = to!string(map!toUpper(fastaRecords("U00096.2.fas").front.sequence));

    auto taskpool = new TaskPool(2);
    scope(exit) taskpool.finish();
    auto reads = BamFile(args[1], taskpool).alignments;

    foreach (read; reads) {
        if (read.is_unmapped) continue;
       
        auto pos = read.position;

        auto seq = read.sequence;
        auto offset = 0;
        auto rev = read.is_reverse_strand;
        auto len = read.sequence_length;

        foreach (j; 0 .. len)
            .reads[j] += 1;

        foreach (op; read.cigar) {
            if (op == 'X') {
                foreach (j; pos .. pos + op.length)
                    update!mismatches(offset + j - pos, rev, len);
            } else if (op == 'M') {
                foreach (j; pos .. pos + op.length)
                    if (reference[j] != seq[offset + j - pos])
                        update!mismatches(offset + j - pos, rev, len);
            } else if (op == 'D') {
                update!deletions(offset, rev, len);
            } else if (op == 'I') {
                update!insertions(offset, rev, len);
            }

            if (op.is_reference_consuming)
                pos += op.length;
            if (op.is_query_consuming)
                offset += op.length;
        }
    }
          
    writeln("offset\tmismatches\tinsertions\tdeletions\treads");
    foreach (i; 0 .. 1024) {
        if (.reads[i] > 0)
            writeln(i, '\t', mismatches[i], '\t', insertions[i], '\t', deletions[i], '\t', .reads[i]);
    }
}
