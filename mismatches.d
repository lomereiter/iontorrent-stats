/**
  rdmd -Ipath/to/sambamba -O -release -inline mismatches.d B7-295.bam > mismatches.dat
*/
import bamfile;

import baseinfo;
import fz.flowcall;

import BioD.Base;
import BioD.Fasta;

import std.stdio;
import std.ascii;
import std.parallelism;
import std.typecons;
import std.algorithm;
import std.conv;

immutable MAX_LEN = 1024;

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

            // copy everything into a static array, because having random access
            // makes read information much easier to work with
            auto base_range = basesWith!"FZ"(read, arg!"FZ"(flow_order));
            typeof(base_range.front)[MAX_LEN * 2] bases_buf = void;
            size_t n_bases;
            while (!base_range.empty) {
                bases_buf[n_bases++] = base_range.front;
                base_range.popFront();
            }
            auto bases = bases_buf[0 .. n_bases];

            // buffer to store mismatching base indices
            size_t[MAX_LEN] indices = void;
            size_t n_indices;

            foreach (size_t i, base; bases) {
                if (base.cigar_operation == 'M' && reference[base.position] != base)
                    indices[n_indices++] = i;
            }

            for (size_t i = 0; i < n_indices; ++i) {
                auto index = indices[i];
                auto base = bases[index];
                auto pos = base.position;
                auto strand = read.strand;

                auto ref_nuc = Base(reference[base.position]);
                auto curr_flow_call = base.flow_call;

                Nullable!ReadFlowCall prev_flow_call;
                Nullable!ReadFlowCall next_flow_call;

                for (long j = index - 1; j >= 0; --j) {
                    if (bases[j].flow_call != base.flow_call) {
                        prev_flow_call = bases[j].flow_call;
                        break;
                    }
                }

                for (long j = index + 1; j < bases.length; ++j) {
                    if (bases[j].flow_call != base.flow_call) {
                        next_flow_call = bases[j].flow_call;
                        break;
                    }
                }

                if (strand == '-') {
                    ref_nuc = ref_nuc.complement;
                }

                writeln(pos, '\t', strand, '\t', index, '\t',
                        ref_nuc, '\t', 
                        prev_flow_call.isNull ? "NA" : to!string(prev_flow_call.base), '\t',
                        prev_flow_call.isNull ? "NA" : to!string(prev_flow_call.intensity), '\t',
                        prev_flow_call.isNull ? "NA" : to!string(prev_flow_call.length), '\t',
                        curr_flow_call.base, '\t', 
                        curr_flow_call.intensity, '\t',
                        curr_flow_call.length, '\t',
                        next_flow_call.isNull ? "NA" : to!string(next_flow_call.base), '\t',
                        next_flow_call.isNull ? "NA" : to!string(next_flow_call.intensity), '\t',
                        next_flow_call.isNull ? "NA" : to!string(next_flow_call.length));
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
    writeln("pos\tstrand\toffset\tref.nuc\tprev.nuc\tprev.intensity\tprev.len\tnuc\tintensity\tlen\tnext.nuc\tnext.intensity\tnext.len");
    mismatchFinder.process();
}
