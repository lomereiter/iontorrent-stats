module accumulators.insertionstats;

import std.stdio;

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

    void updateStatistics(Insertion)(Insertion insertion)
    {
        ++_n_insertions;
        _n_overcalls += insertion.number_of_flowcalls;

        with (insertion)
        {
            if (number_of_flowcalls == 1)
            {
                _n_resolved += 1;
                _intensities[first_flow_call.intensity_value] += 1;
            }
            else if (number_of_flowcalls == 2)
            {
                auto has_left = !previous_flow_call.isNull;
                auto has_right = !next_flow_call.isNull;

                if ((has_left && previous_flow_call == first_flow_call) ||
                    (has_right && next_flow_call == last_flow_call))
                {
                    _n_resolved += 1;

                    _intensities[first_flow_call.intensity_value] += 1;
                    _intensities[last_flow_call.intensity_value] += 1;
                }
            }
            // There are also cases like 2.61, .56, 1.53, which are not counted here.
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
