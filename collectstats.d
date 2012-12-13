import bio.bam.reader;
import bio.bam.pileup;
import bio.bam.baseinfo;

import bio.core.base;
import bio.core.tinymap;

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
        _out.writeln("# ref: reference base at this position");
        _out.writeln("# coverage: number of mapped reads overlapping the site");
        _out.writeln("# mismatches: number of mismatches at this position");
        _out.writeln("# deletions: number of deletions starting at this position");
        _out.writeln("# insertions: number of insertions starting after this position");
    }

    private void writeHeader()
    {
        _out.writeln("pos\tref\tcoverage\tmismatches\tdeletions\tinsertions");
    }

    void printColumn(Column)(ref Column column)
    {
        int n_mismatches_at;
        int n_deletions_starting_at;
        int n_insertions_after;

        foreach (read; column.reads)
        {
            switch (read.current_cigar_operation.type)
            {
                case 'M', 'X':
                    if (read.current_base != column.reference_base)
                        n_matches_at += 1;
                    break;
                case 'D':
                    // count each deletion only once
                    if (read.cigar_operation_offset == 0) 
                        n_deletions_starting_at += 1;
                    break;
                default:
                    break;
            }

            if (!read.cigar_after.empty && read.cigar_after.front.type == 'I')
                n_insertions_after += 1;
        }

        _out.writeln(column.position + 1, // 0-based -> 1-based
                     '\t',
                     column.reference_base,
                     '\t',
                     column.coverage,
                     '\t',
                     n_mismatches_at,
                     '\t',
                     n_deletions_starting_at,
                     '\t',
                     n_insertions_after);
    }
}

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

    auto column_stats_printer = new ColumnStatsPrinter("columnstats.dat");

    foreach (column; makePileup(bam, true))
    {
        column_stats_printer.printColumn(column);
    }
}
