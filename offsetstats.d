module offsetstats;

import std.stdio;
import std.array;

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
    
    void updateStatistics(Column)(Column column)
    {
        foreach (read; column.reads)
        {
            if (read.query_offset == 0 && read.cigar_operation.is_reference_consuming)
            {
                _total_reads[0 .. read.sequence_length - 1] += 1;
            }

            switch (read.cigar_operation.type)
            {
                case 'M', 'X':
                    if (read.current_base != column.reference_base)
                    {
                        auto index = read.query_offset;
                        if (read.is_reverse_strand)
                            index = read.sequence_length - index - 1;
                        _mismatches_at[index] += 1;
                    }
                    break;
                case 'D':
                    if (read.cigar_operation_offset == 0)
                    {
                        auto index = read.query_offset;
                        if (read.is_reverse_strand)
                        {
                            index = read.sequence_length - index;
                        }
                        _deletions_before[index] += 1;
                    }
                    break;
                default:
                    break;
            }

            if (!read.cigar_before.empty && read.cigar_before.back.type == 'I')
            {
                if (read.cigar_operation_offset == 0)
                {
                    auto index = read.query_offset - read.cigar_before.back.length;
                    if (read.is_reverse_strand)
                    {
                        index = read.sequence_length - read.query_offset;
                    }
                    _insertions_starting_at[index] += 1;
                }
            }
        }
    }
}

