module accumulators.insertionstats;

import std.stdio;
import std.algorithm;
import std.range;

class InsertionStatsAccumulator
{
    private
    {
        uint _n_insertions;
        uint _n_resolved; // likely to be caused solely by overcalls

        uint _n_overcalls;

        uint[] _intensities; // overcall intensities
    }

    this(size_t max_intensity_value=1536)
    {
        _intensities = new uint[max_intensity_value];
    }

    static auto merge(InsertionStatsAccumulator acc1,
                      InsertionStatsAccumulator acc2)
    {
        auto max_intensity_value = max(acc1._intensities.length,
                                       acc2._intensities.length);

        auto acc = new InsertionStatsAccumulator(max_intensity_value);
        acc._n_insertions = acc1._n_insertions + acc2._n_insertions;
        acc._n_resolved = acc1._n_resolved + acc2._n_resolved;
        acc._n_overcalls = acc1._n_overcalls + acc2._n_overcalls;
        acc._intensities[0 .. acc1._intensities.length] = acc1._intensities[];
        acc._intensities[0 .. acc2._intensities.length] += acc2._intensities[];

        return acc;
    }

    void updateStatistics(Insertion)(Insertion insertion)
    {
        ++_n_insertions;
        _n_overcalls += insertion.number_of_flowcalls;

        with (insertion)
        {
            auto has_left = !previous_flow_call.isNull;
            auto has_right = !next_flow_call.isNull;

            if (number_of_flowcalls == 1)
            {
                _n_resolved += 1;
                _intensities[first_flow_call.intensity_value] += 1;
            }
            else if (number_of_flowcalls == 2)
            {
                if ((has_left && previous_flow_call == first_flow_call) ||
                    (has_right && next_flow_call == last_flow_call))
                {
                    _n_resolved += 1;

                    _intensities[first_flow_call.intensity_value] += 1;
                    _intensities[last_flow_call.intensity_value] += 1;
                }
            }
            else if (number_of_flowcalls == 3)
            {
                if ((has_left && previous_flow_call == first_flow_call) &&
                    (has_right && next_flow_call == last_flow_call))
                {
                    _n_resolved += 1;

                    _intensities[first_flow_call.intensity_value] += 1;
                    _intensities[last_flow_call.intensity_value] += 1;
                    _intensities[flowcalls.drop(1).front.intensity_value] += 1;
                }
            }
            // There are also cases like 2.61, .56, .52, which are not counted here.
            // However, they are quite rare and don't influence the results
            // significantly. Overwhelming majority of insertions has length 1 or 2.
        }
    }

    void printSummary(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# Number of insertions: ", _n_insertions);
        _file.writeln("# Likely to be caused solely by overcalls: ", _n_resolved,
                      " (", cast(double)_n_resolved * 100.0 / _n_insertions, "%)");

        _file.writeln("# Number of overcalls: ", _n_overcalls);
    }


    void printOvercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# intensity: signal intensity value for overcalled homopolymer");
        _file.writeln("# count: number of overcalls with this intensity value");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }
}
