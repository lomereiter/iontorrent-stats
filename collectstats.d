import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.baseinfo;

import bio.core.base;
import bio.core.tinymap;

import std.stdio;

import columnstats;
import offsetstats;

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

    auto column_stats_printer = new ColumnStatsPrinter("columns.dat");
    auto offset_stats_accumulator = new OffsetStatsAccumulator();

    foreach (column; makePileup(bam.reads, true))
    {
        column_stats_printer.printColumn(column);
        offset_stats_accumulator.updateStatistics(column);
    }

    offset_stats_accumulator.printReport("offsets.dat");
}
