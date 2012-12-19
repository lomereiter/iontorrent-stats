module accumulators.insertionstats;

import std.stdio;

// TODO: rewrite using approach similar to that of DeletionStatsAccumulator
class InsertionStatsAccumulator
{
    private
    {
        uint[] _neighbour_intensities; // inserted bases are from a single flow call
        ushort[][] _intensity_pairs;       // inserted bases are from two flow calls
        uint _total_insertions;
        uint _insertions_of_type_1; // single overcall
        uint _insertions_of_type_2; // 'orphan' call
        uint _insertions_of_type_3; // two neighbor overcalls
        uint[] _orphan_intensities;     // calls of unexistent bases

        uint[] _intensities; // unifies all three types of errors
    }

    this(size_t max_intensity_value=1536)
    {
        _neighbour_intensities = new uint[max_intensity_value];
        _orphan_intensities = new uint[max_intensity_value];
        _intensities = new uint[max_intensity_value];
        _intensity_pairs = new ushort[][](max_intensity_value, max_intensity_value);
    }

    void updateStatistics(Insertion)(Insertion insertion)
    {
        ++_total_insertions;

        if (insertion.number_of_flowcalls == 1)
        {
            if (!insertion.previous_flow_call.isNull &&
                insertion.first_flow_call == insertion.previous_flow_call)
            {
                with (insertion)
                {
                    _neighbour_intensities[previous_flow_call.intensity_value] += 1;
                    _intensities[previous_flow_call.intensity_value] += 1;
                    _insertions_of_type_1 += 1;
                }
            }
            else if (!insertion.next_flow_call.isNull &&
                     insertion.last_flow_call == insertion.next_flow_call)
            {
                with (insertion)
                {
                    _neighbour_intensities[next_flow_call.intensity_value] += 1;
                    _intensities[next_flow_call.intensity_value] += 1;
                    _insertions_of_type_1 += 1;
                }
            }
            else
            {
                with (insertion)
                {
                    _insertions_of_type_2 += 1;
                    _orphan_intensities[flowcalls.front.intensity_value] += 1;
                    _intensities[flowcalls.front.intensity_value] += 1;
                }
            }
        }

        if (insertion.number_of_flowcalls == 2 &&
            !insertion.previous_flow_call.isNull &&
            !insertion.next_flow_call.isNull &&
            insertion.first_flow_call == insertion.previous_flow_call &&
            insertion.last_flow_call == insertion.next_flow_call)
        {
            with (insertion)
            {
                auto left = previous_flow_call.intensity_value;
                auto right = next_flow_call.intensity_value;

                _intensity_pairs[left][right] += 1;
                _insertions_of_type_3 += 1;

                _intensities[left] += 1;
                _intensities[right] += 1;
            }
        }
    }

    void printNeighbourSummary(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# total number of insertions: ", _total_insertions);

        _file.writeln("# number of one-sided overcalls: ",
                      _insertions_of_type_1, 
                      " (",
                      cast(double)_insertions_of_type_1 * 100.0 / _total_insertions,
                      "%)");

        _file.writeln("# number of 'orphans' (single flowcalls different from both neighbors): ",
                      _insertions_of_type_2, 
                      " (",
                      cast(double)_insertions_of_type_2 * 100.0 / _total_insertions,
                      "%)");

        _file.writeln("# number of two-sided overcalls: ",
                      _insertions_of_type_3,
                      " (",
                      cast(double)_insertions_of_type_3 * 100.0 / _total_insertions,
                      "%)");
    }

    void printOrphansReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# There are a lot of insertions that consist of a single flow call");
        _file.writeln("# which is different from both neighbour flow calls. As it turns out");
        _file.writeln("# most of these insertions have signal intensity close to 50.");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _orphan_intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }

    void printOneSidedOvercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# intensity: signal intensity value for overcalled homopolymer");
        _file.writeln("# count: number of overcalls with this intensity value");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _neighbour_intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }

    void printTwoSidedOvercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# Usually aligner merges two neighboring overcalls into one insertion");
        _file.writeln("# because one insertion is 'cheaper' than two, in its view.");
        _file.writeln("# ");
        _file.writeln("# first: intensity of the first flow call");
        _file.writeln("# second: intensity of the second flow call");

        _file.writeln("first\tsecond");

        foreach (size_t first, dist; _intensity_pairs)
        {
            foreach (size_t second, count; dist)
                for (ushort i = 0; i < count; ++i)
                    _file.writeln(first, '\t', second);
        }
    }

    void printOvercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# intensity: round(normalized signal intensity * 100.0)");
        _file.writeln("# count: number of overcalls with this intensity");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }
}
