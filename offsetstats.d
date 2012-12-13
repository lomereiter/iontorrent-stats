module offsetstats;

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
        _total_reads[0 .. bases.length - 1] += 1;

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
