module processor;

import logger;

import bio.bam.baseinfo;

import events.insertion;
import events.deletion;
import events.mismatch;

import printers.columnstats;

import accumulators.offsetstats;
import accumulators.flowstats;
import accumulators.insertionstats;
import accumulators.deletionstats;
import accumulators.mismatchstats;

import std.typetuple;
import std.algorithm;
import std.exception;

class PileupProcessor(Pileup)
{
    private 
    {
        Pileup _pileup;
    }

    public 
    {
        bool collect_insertion_stats;
        bool collect_deletion_stats;
        bool collect_mismatch_stats;
        bool collect_flow_stats;
        bool collect_offset_stats;
        bool collect_column_stats;
            
        ColumnStatsPrinter column_stats_printer;
        OffsetStatsAccumulator offset_stats_accumulator;
        FlowStatsAccumulator flow_stats_accumulator;
        InsertionStatsAccumulator insertion_stats_accumulator;
        DeletionStatsAccumulator deletion_stats_accumulator;
        MismatchStatsAccumulator mismatch_stats_accumulator;

        string flow_order;
        string key_sequence;
        string column_stats_filename = "columns.dat";

        ulong id;
    }

    this(Pileup pileup)
    {
        _pileup = pileup;
    }

    void mergeResultsWith(PileupProcessor other)
    {
        offset_stats_accumulator = OffsetStatsAccumulator.merge(
                                       offset_stats_accumulator,
                                       other.offset_stats_accumulator);

        flow_stats_accumulator = FlowStatsAccumulator.merge(
                                     flow_stats_accumulator,
                                     other.flow_stats_accumulator);

        insertion_stats_accumulator = InsertionStatsAccumulator.merge(
                                          insertion_stats_accumulator,
                                          other.insertion_stats_accumulator);

        deletion_stats_accumulator = DeletionStatsAccumulator.merge(
                                         deletion_stats_accumulator,
                                         other.deletion_stats_accumulator);

        mismatch_stats_accumulator = MismatchStatsAccumulator.merge(
                                         mismatch_stats_accumulator,
                                         other.mismatch_stats_accumulator);
    }

    PileupProcessor settings() @property 
    { 
        return this; 
    }

    private bool collect_some_read_stats() @property const 
    {
        return collect_flow_stats      || collect_offset_stats ||
               collect_insertion_stats || collect_deletion_stats ||
               collect_mismatch_stats;
    }

    private bool collect_nothing() @property const 
    {
        return !collect_column_stats && !collect_some_read_stats;
    }

    void run()
    {
        if (collect_nothing)
            return;

        if (collect_column_stats)
        {
            column_stats_printer = new ColumnStatsPrinter(column_stats_filename, id == 0);
        }

        if (collect_flow_stats)
        {
            flow_stats_accumulator = new FlowStatsAccumulator();
        }

        if (collect_offset_stats)
        {
            offset_stats_accumulator = new OffsetStatsAccumulator();
        }

        if (collect_insertion_stats)
        {
            insertion_stats_accumulator = new InsertionStatsAccumulator();
        }

        if (collect_deletion_stats)
        {
            deletion_stats_accumulator = new DeletionStatsAccumulator();
        }

        if (collect_mismatch_stats)
        {
            mismatch_stats_accumulator = new MismatchStatsAccumulator();
        }

        foreach (column; _pileup)
        {
            if (collect_column_stats)
            {
                column_stats_printer.printColumn(column);
            }

            if (!collect_some_read_stats)
            {
                continue;
            }

            foreach (read; column.reads_starting_here)
            {
                try 
                {
                    alias TypeTuple!("FZ", "MD", Option.mdNextOp, Option.cigarExtra) Options;
                    auto bases = basesWith!Options(read, arg!"flowOrder"(flow_order), 
                                                         arg!"keySequence"(key_sequence));

                    typeof(bases.front)[1024] baseinfo_buf = void;
                    size_t i;
                    while (!bases.empty)
                    {
                        bases.constructFront(&baseinfo_buf[i++]);
                        bases.popFront();
                    }

                    auto baseinfo = baseinfo_buf[0 .. i];

                    if (collect_insertion_stats)
                    {
                        foreach (insertion; insertionEvents(baseinfo))
                        {
                            insertion_stats_accumulator.updateStatistics(insertion);
                            if (collect_flow_stats)
                            {
                                flow_stats_accumulator.updateInsertionStats(insertion);
                            }
                        }
                    }

                    if (collect_deletion_stats)
                    {
                        auto tmp = read["FZ"];
                        if (tmp.is_nothing)
                            tmp = read["ZM"];

                        enforce(!(tmp.is_nothing), "FZ or ZM tag must be presented in a mapped read");
                        auto intensities = *(cast(short[]*)(&tmp));

                        foreach (deletion; deletionEvents(baseinfo, flow_order, intensities))
                        {
                            deletion_stats_accumulator.updateStatistics(deletion);
                            if (collect_flow_stats)
                            {
                                flow_stats_accumulator.updateDeletionStats(deletion);
                            }
                        }
                    }

                    if (collect_mismatch_stats)
                    {
                        foreach (mismatch; mismatchEvents(baseinfo))
                        {
                            mismatch_stats_accumulator.updateStatistics(mismatch);
                            if (collect_flow_stats)
                            {
                                flow_stats_accumulator.updateMismatchStats(mismatch);
                            }
                        }
                    }

                    if (collect_offset_stats)
                    {
                        offset_stats_accumulator.updateStatistics(baseinfo);
                    }

                    if (collect_flow_stats)
                    {
                        flow_stats_accumulator.updateStatistics(baseinfo);
                    }
                }
                catch (Exception e) {
                    Logger.warn("couldn't process read with name '" ~
                                read.name ~ "' : " ~ e.msg);                       
                }
            }
        }
    }
}
