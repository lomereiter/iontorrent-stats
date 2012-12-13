module columnstats;

import std.array;
import std.stdio;

class ColumnStatsPrinter 
{
    private 
    {
        File _out;
    }

    this(string filename) 
    {
        _out = File(filename, "w+");
        writeComments();
        writeHeader();
    }

    private void writeComments()
    {
        _out.writeln("# pos: 1-based position on the reference");
        _out.writeln("# ref: reference base at the position");
        _out.writeln("# coverage: number of mapped reads overlapping the site");
        _out.writeln("# mismatches: number of mismatches at the position");
        _out.writeln("# deletions: number of deletions starting at the position");
        _out.writeln("# insertions: number of insertions ending right before the position");
    }

    private void writeHeader()
    {
        _out.writeln("pos\tref\tcoverage\tmismatches\tdeletions\tinsertions");
    }

    void printColumn(Column)(ref Column column)
    {
        int n_mismatches_at;
        int n_deletions_starting_at;
        int n_insertions_before;

        foreach (read; column.reads)
        {
            switch (read.cigar_operation.type)
            {
                case 'M', 'X':
                    if (read.current_base != column.reference_base)
                        n_mismatches_at += 1;
                    break;
                case 'D':
                    // count each deletion only once
                    if (read.cigar_operation_offset == 0) 
                        n_deletions_starting_at += 1;
                    break;
                default:
                    break;
            }

            // if we're at the beginning of the current CIGAR operation
            // check for insertions before it
            if (read.cigar_operation_offset == 0)
            {
                if (!read.cigar_before.empty && read.cigar_before.back.type == 'I')
                {
                    n_insertions_before += 1;
                }
            }
        }

        _out.writeln(column.position + 1, '\t', // 0-based -> 1-based
                     column.reference_base, '\t',
                     column.coverage, '\t',
                     n_mismatches_at, '\t',
                     n_deletions_starting_at, '\t',
                     n_insertions_before);
    }
}
