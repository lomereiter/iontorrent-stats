import bio.bam.reader;
import bio.bam.pileup;

import std.stdio;
import std.getopt;
import std.range;
import std.file;
import std.path;

version(parallel)
{
    import std.parallelism;
}

import processor;

__gshared bool no_insertion_stats = false;
__gshared bool no_deletion_stats = false;
__gshared bool no_flow_stats = false;
__gshared bool no_offset_stats = false;
__gshared bool no_column_stats = false;

__gshared string flow_order;
__gshared string key_sequence;

__gshared string dir; // where to output reports

auto createPileupProcessor(Pileup)(Pileup pileup, ulong id, string column_stats_fn)
{
    auto processor = new PileupProcessor!Pileup(pileup);

    processor.id = id;
    processor.flow_order = flow_order;
    processor.key_sequence = key_sequence;

    processor.settings.column_stats_filename   = column_stats_fn;
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
        version (parallel)
        {
            int chunk_size = 32_000_000; 

            getopt(args,
                   std.getopt.config.caseSensitive,
                   "output-dir|d", &dir,
                   "no-insertion-stats|I", &no_insertion_stats,
                   "no-deletion-stats|D",  &no_deletion_stats,
                   "no-flow-stats|F",      &no_flow_stats,
                   "no-offset-stats|O",    &no_offset_stats,
                   "no-column-stats|C",    &no_column_stats,
                   "chunk-size|s",         &chunk_size); // for testing only
        }
        else
        {
            getopt(args,
                   std.getopt.config.caseSensitive,
                   "output-dir|d", &dir,
                   "no-insertion-stats|I", &no_insertion_stats,
                   "no-deletion-stats|D",  &no_deletion_stats,
                   "no-flow-stats|F",      &no_flow_stats,
                   "no-offset-stats|O",    &no_offset_stats,
                   "no-column-stats|C",    &no_column_stats);
        }

        if (args.length < 2) {
            printUsage(args[0]);
            return 0;
        }

        auto filename = args[1];

        version (parallel)
        {
            auto task_pool = new TaskPool(totalCPUs);
            scope(exit) task_pool.finish();

            auto bam = new BamReader(filename, task_pool);
        }
        else // default settings
        {
            auto bam = new BamReader(filename);
        }
 
        ////////////////////////////////////////////////////////////////////////////////////////////
        ///                         get flow order and key sequence                                 
        ////////////////////////////////////////////////////////////////////////////////////////////
        auto rg = bam.header.read_groups.values.front;
        flow_order = rg.flow_order;
        key_sequence = rg.key_sequence;

        ////////////////////////////////////////////////////////////////////////////////////////////
        ///                             change working directory                                    
        ////////////////////////////////////////////////////////////////////////////////////////////
        if (dir !is null && !std.file.exists(dir))
        {
            std.file.mkdirRecurse(dir);
        }

        version (parallel) // parallel version
        {
            auto pileups = pileupChunks(bam.reads, true, chunk_size);

            static auto createAndRunProcessor(T)(T pileup_and_id)
            {
                auto id = pileup_and_id[1];
                auto column_stats_fn = buildPath(dir, "columns"~to!string(id)~".dat");

                auto processor = createPileupProcessor(pileup_and_id[0], id, column_stats_fn);

                ////////////////////////////////////////////////////////////////////////////////////
                ///                                 process pileup                                  
                ////////////////////////////////////////////////////////////////////////////////////
                processor.run();

                return processor;
            }

            auto labeled_pileups = zip(pileups, sequence!"n"());
            auto processors = task_pool.map!createAndRunProcessor(labeled_pileups, totalCPUs, 1);

            ////////////////////////////////////////////////////////////////////////////////////////
            ///                                 merge processor results                             
            ////////////////////////////////////////////////////////////////////////////////////////
            auto processor = processors.front;
           
            if (!no_column_stats)
            {
                processor.column_stats_printer.closeFileHandle();
                std.file.rename(processor.column_stats_filename, buildPath(dir, "columns.dat"));
            }

            processors.popFront();

            foreach (p; processors)
            {
                processor.mergeResultsWith(p);

                if (!no_column_stats)
                {
                    p.column_stats_printer.closeFileHandle();
                    auto contents = std.file.read(p.column_stats_filename);
                    std.file.append(buildPath(dir, "columns.dat"), contents);
                    std.file.remove(p.column_stats_filename);
                }
            }
        }
        else // serial version
        {
            auto pileup = makePileup(bam.reads, true);
            auto processor = createPileupProcessor(pileup, 0, buildPath(dir, "columns.dat"));
            processor.run();
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        ///                                 output results                                          
        ////////////////////////////////////////////////////////////////////////////////////////////
        if (!no_offset_stats)
        {
            processor.offset_stats_accumulator.printReport(buildPath(dir, "offsets.dat"));
        }

        if (!no_flow_stats)
        {
            processor.flow_stats_accumulator.printReport(buildPath(dir, "flows.dat"));
        }

        if (!no_insertion_stats)
        {
            processor.insertion_stats_accumulator.printSummary("/dev/stdout");
            auto fn = buildPath(dir, "overcall.intensities.dat");
            processor.insertion_stats_accumulator.printOvercallsReport(fn);
        }

        if (!no_deletion_stats)
        {
            processor.deletion_stats_accumulator.printSummary("/dev/stdout");
            auto fn = buildPath(dir, "undercall.intensities.dat");
            processor.deletion_stats_accumulator.printUndercallsReport(fn);
        }
    }
    catch (Exception e)
    {
        stderr.writeln(e);
        throw e;
        return 1;
    }

    return 0;
}
