module accumulators.offsetstats;

import std.stdio;
import std.array;
import std.math : abs;

class OffsetStatsAccumulator {

    private 
    {
        int[1024] _total_reads;
        int[1024] _mismatches_at;
        int[1024] _deletions_before;
        int[1024] _insertions_starting_at;
    }

    static auto merge(OffsetStatsAccumulator acc1, OffsetStatsAccumulator acc2)
    {
        auto acc = new OffsetStatsAccumulator();

        acc._total_reads[] = acc1._total_reads[] + acc2._total_reads[];
        acc._mismatches_at[] = acc1._mismatches_at[] + acc2._mismatches_at[];
        acc._deletions_before[] = acc1._deletions_before[] + acc2._deletions_before[];
        acc._insertions_starting_at[] = acc1._insertions_starting_at[] + 
                                        acc2._insertions_starting_at[];

        return acc;
    }

    void printReport(string filename)
    {
        auto _out = File(filename, "w+");

        _out.writeln("# offset: read sequence offset, taking strand into account");
        _out.writeln("# total: number of reads with length at least (offset + 1)");
        _out.writeln("# mismatches: number of mismatches at the offset");
        _out.writeln("# deletions: number of deletions ending right before the offset");
        _out.writeln("# insertions: number of insertions starting at the offset");

        _out.writeln("offset\ttotal\tmismatches\tdeletions\tinsertions");

        foreach (size_t i, count; _total_reads)
        {
            if (count > 0)
            {
                _out.writeln(i, '\t', 
                             count, '\t',
                             _mismatches_at[i], '\t',
                             _deletions_before[i], '\t',
                             _insertions_starting_at[i]);
            }
        }
    }

    void updateStatistics(BaseInfo)(BaseInfo[] bases)
    {
        for (size_t i = 0; i < bases.length - 1; i++)
            _total_reads[i] += 1;

        foreach (size_t offset, baseinfo; bases)
        {
            if (baseinfo.cigar_operation.is_match_or_mismatch &&
                baseinfo.reference_base != baseinfo.base)
            {
                _mismatches_at[offset] += 1; 
            }
            
            if (baseinfo.cigar_operation.type == 'I' && 
                baseinfo.cigar_operation_offset == 0)
            {
                _insertions_starting_at[offset] += 1;
            }

            if (offset > 0 && 
                (cast(long)baseinfo.position - cast(long)bases[offset - 1].position).abs > 1)
            {
                _deletions_before[offset] += 1;
            }
        }
    }
}
