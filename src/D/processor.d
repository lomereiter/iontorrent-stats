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
        ColumnStatsPrinter _column_stats_printer;

        bool _collect_column_stats;
        bool _collect_deletion_stats;
        bool _collect_insertion_stats;
        bool _collect_flow_stats;
        bool _collect_offset_stats;
    }

    public 
    {
        OffsetStatsAccumulator offset_stats_accumulator;
        FlowStatsAccumulator flow_stats_accumulator;
        InsertionStatsAccumulator insertion_stats_accumulator;
        DeletionStatsAccumulator deletion_stats_accumulator;

        string flow_order;
        string key_sequence;
    }

    this(Pileup pileup)
    {
        _pileup = pileup;
    }

    PileupProcessor settings() @property 
    { 
        return this; 
    }

    bool collect_column_stats() @property const
    {
        return _collect_column_stats;
    }

    void collect_column_stats(bool collect) @property 
    {
        _collect_column_stats = collect;

        if (collect && _column_stats_printer is null)
            _column_stats_printer = new ColumnStatsPrinter("columns.dat");
    }

    bool collect_offset_stats() @property const
    {
        return _collect_offset_stats;
    }

    void collect_offset_stats(bool collect) @property 
    {
        _collect_offset_stats = collect;

        if (collect && offset_stats_accumulator is null)
            offset_stats_accumulator = new OffsetStatsAccumulator();
    }

    bool collect_flow_stats() @property const
    {
        return _collect_flow_stats;
    }

    void collect_flow_stats(bool collect) @property 
    {
        _collect_flow_stats = collect;

        if (collect && flow_stats_accumulator is null)
            flow_stats_accumulator = new FlowStatsAccumulator();
    }

    bool collect_insertion_stats() @property const
    {
        return _collect_insertion_stats;
    }

    void collect_insertion_stats(bool collect) @property 
    {
        _collect_insertion_stats = collect;

        if (collect && insertion_stats_accumulator is null)
            insertion_stats_accumulator = new InsertionStatsAccumulator();
    }

    bool collect_deletion_stats() @property const
    {
        return _collect_deletion_stats;
    }

    void collect_deletion_stats(bool collect) @property 
    {
        _collect_deletion_stats = collect;

        if (collect && deletion_stats_accumulator is null)
            deletion_stats_accumulator = new DeletionStatsAccumulator();
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

        foreach (column; _pileup)
        {
            if (collect_column_stats)
            {
                _column_stats_printer.printColumn(column);
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
