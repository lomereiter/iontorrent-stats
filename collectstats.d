import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.baseinfo;

import bio.core.base;
import bio.core.tinymap;

import std.stdio;
import std.typetuple;

import events.insertion;
import events.deletion;
import printers.columnstats;
import printers.insertioninfo;
import printers.deletioninfo;
import accumulators.offsetstats;
import accumulators.flowstats;
import accumulators.insertionstats;
import accumulators.deletionstats;

immutable collect_insertion_stats = false;
immutable collect_deletion_stats = true;

void printUsage() {
    stderr.writeln("usage: ./collectstats <input.bam>");
    stderr.writeln();
    stderr.writeln("       BAM file must provide FZ, ZF, and MD tags");
}

void main(string[] args) {
    if (args.length < 2) {
        printUsage();
    }

    auto filename = args[1];

    auto bam = new BamReader(filename);
    auto rg = bam.header.read_groups.values.front;
    auto flow_order = rg.flow_order;
    auto key_sequence = rg.key_sequence;

    auto column_stats_printer = new ColumnStatsPrinter("columns.dat");
    auto insertion_info_printer = new InsertionInfoPrinter("insertions.dat");
    auto deletion_info_printer = new DeletionInfoPrinter("deletions.dat", flow_order);

    auto offset_stats_accumulator = new OffsetStatsAccumulator();
    auto flow_stats_accumulator = new FlowStatsAccumulator();
    auto insertion_stats_accumulator = new InsertionStatsAccumulator();
    auto deletion_stats_accumulator = new DeletionStatsAccumulator();

    foreach (column; makePileup(bam.reads, true))
    {
        column_stats_printer.printColumn(column);

        foreach (read; column.reads_starting_here)
        {
            alias TypeTuple!("FZ", "MD", Option.mdNextOp) Options;
            auto bases = basesWith!Options(read, arg!"flowOrder"(flow_order), 
                                                 arg!"keySequence"(key_sequence));

            typeof(bases.front)[1024] baseinfo_buf = void;
            size_t i;
            foreach (info; bases)
                baseinfo_buf[i++] = info;

            auto baseinfo = baseinfo_buf[0 .. i];

            if (collect_insertion_stats)
            {
                foreach (insertion; insertionEvents(baseinfo))
                {
                    insertion_info_printer.printInsertion(insertion);
                    insertion_stats_accumulator.updateStatistics(insertion);
                }
            }

            if (collect_deletion_stats)
            {
                foreach (deletion; deletionEvents(baseinfo))
                {
                    auto tmp = read["FZ"];
                    auto intensities = *(cast(ushort[]*)(&tmp));
                    deletion_info_printer.printDeletion(deletion, intensities);
                    deletion_stats_accumulator.updateStatistics(deletion);
                }
            }

            offset_stats_accumulator.updateStatistics(baseinfo);
            flow_stats_accumulator.updateStatistics(baseinfo);
        }
    }

    offset_stats_accumulator.printReport("offsets.dat");
    flow_stats_accumulator.printReport("flows.dat");

    insertion_stats_accumulator.printNeighbourSummary("/dev/stdout");
    insertion_stats_accumulator.printOvercallsReport("overcall.intensities.dat");

    deletion_stats_accumulator.printOneSidedUndercallsReport("undercall.intensities.dat");
}
