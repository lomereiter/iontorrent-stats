import bio.bam.reader;
import bio.bam.pileup;

import std.stdio;
import std.getopt;
import std.file;

import processor;

__gshared bool no_insertion_stats = false;
__gshared bool no_deletion_stats = false;
__gshared bool no_flow_stats = false;
__gshared bool no_offset_stats = false;
__gshared bool no_column_stats = false;

__gshared string flow_order;
__gshared string key_sequence;

auto createPileupProcessor(Pileup)(Pileup pileup)
{
    auto processor = new PileupProcessor!Pileup(pileup);

    processor.flow_order = flow_order;
    processor.key_sequence = key_sequence;

    processor.settings.collect_column_stats    = !no_column_stats;
    processor.settings.collect_deletion_stats  = !no_deletion_stats;
    processor.settings.collect_insertion_stats = !no_insertion_stats;
    processor.settings.collect_flow_stats      = !no_flow_stats;
    processor.settings.collect_offset_stats    = !no_offset_stats;

    return processor;
}

void printUsage(string prg_name) {
    stderr.writeln("usage: ", prg_name, " [OPTIONS] <input.bam>");
    stderr.writeln();
    stderr.writeln("       BAM file must provide FZ, ZF, and MD tags");
    stderr.writeln();
    
    stderr.writeln("OPTIONS:    -d, --output-dir=DIR");
    stderr.writeln("                 Directory to which results will be output (created if doesn't exist).");
    stderr.writeln("            -I, --no-insertion-stats");
    stderr.writeln("                 Don't collect insertion statistics");
    stderr.writeln("            -D, --no-deletion-stats");
    stderr.writeln("                 Don't collect deletion statistics");
    stderr.writeln("            -F, --no-flow-stats");
    stderr.writeln("                 Don't collect flow signal intensity statistics");
    stderr.writeln("            -O, --no-offset-stats");
    stderr.writeln("                 Don't collect read offset statistics");
    stderr.writeln("            -C, --no-column-stats");
    stderr.writeln("                 Don't print column statistics");
}

int main(string[] args) {

    try
    {
        string dir;

        getopt(args,
               std.getopt.config.caseSensitive,
               "output-dir|d", &dir,
               "no-insertion-stats|I", &no_insertion_stats,
               "no-deletion-stats|D",  &no_deletion_stats,
               "no-flow-stats|F",      &no_flow_stats,
               "no-offset-stats|O",    &no_offset_stats,
               "no-column-stats|C",    &no_column_stats);

        if (args.length < 2) {
            printUsage(args[0]);
            return 0;
        }

        auto filename = args[1];

        auto bam = new BamReader(filename);
 
        ////////////////////////////////////////////////////////////////////////////////////////////
        ///                         get flow order and key sequence                                 
        ////////////////////////////////////////////////////////////////////////////////////////////
        auto rg = bam.header.read_groups.values.front;
        flow_order = rg.flow_order;
        key_sequence = rg.key_sequence;

        ////////////////////////////////////////////////////////////////////////////////////////////
        ///                             change working directory                                    
        ////////////////////////////////////////////////////////////////////////////////////////////
        if (dir !is null)
        {
            std.file.mkdirRecurse(dir);
            std.file.chdir(dir);
        }

        auto pileup = makePileup(bam.reads, true);

        ////////////////////////////////////////////////////////////////////////////////////////////
        ///                                 process pileup                                          
        ////////////////////////////////////////////////////////////////////////////////////////////
        auto processor = createPileupProcessor(pileup);
        processor.run();

        ////////////////////////////////////////////////////////////////////////////////////////////
        ///                                 output results                                          
        ////////////////////////////////////////////////////////////////////////////////////////////
        if (!no_offset_stats)
        {
            processor.offset_stats_accumulator.printReport("offsets.dat");
        }

        if (!no_flow_stats)
        {
            processor.flow_stats_accumulator.printReport("flows.dat");
        }

        if (!no_insertion_stats)
        {
            processor.insertion_stats_accumulator.printSummary("/dev/stdout");
            processor.insertion_stats_accumulator.printOvercallsReport("overcall.intensities.dat");
        }

        if (!no_deletion_stats)
        {
            processor.deletion_stats_accumulator.printSummary("/dev/stdout");
            processor.deletion_stats_accumulator.printUndercallsReport("undercall.intensities.dat");
        }
    }
    catch (Exception e)
    {
        stderr.writeln(e);
        return 1;
    }

    return 0;
}
