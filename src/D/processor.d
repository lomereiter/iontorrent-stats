module processor;

import bio.bam.baseinfo;

import events.insertion;
import events.deletion;

import printers.columnstats;

import accumulators.offsetstats;
import accumulators.flowstats;
import accumulators.insertionstats;
import accumulators.deletionstats;

import std.typetuple;
import std.algorithm;

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
        bool collect_flow_stats;
        bool collect_offset_stats;
        bool collect_column_stats;
            
        ColumnStatsPrinter column_stats_printer;
        OffsetStatsAccumulator offset_stats_accumulator;
        FlowStatsAccumulator flow_stats_accumulator;
        InsertionStatsAccumulator insertion_stats_accumulator;
        DeletionStatsAccumulator deletion_stats_accumulator;

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
    }

    PileupProcessor settings() @property 
    { 
        return this; 
    }

    private bool collect_some_read_stats() @property const 
    {
        return collect_flow_stats      || collect_offset_stats ||
               collect_insertion_stats || collect_deletion_stats;
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
                alias TypeTuple!("FZ", "MD", Option.mdNextOp, Option.cigarExtra) Options;
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
                        insertion_stats_accumulator.updateStatistics(insertion);
                    }
                }

                if (collect_deletion_stats)
                {
                    auto tmp = read["FZ"];
                    auto intensities = *(cast(ushort[]*)(&tmp));

                    foreach (deletion; deletionEvents(baseinfo, flow_order, intensities))
                    {
                        deletion_stats_accumulator.updateStatistics(deletion);
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
        }
    }
}
