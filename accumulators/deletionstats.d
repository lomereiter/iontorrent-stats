module accumulators.deletionstats;

import std.stdio;

class DeletionStatsAccumulator
{
    private
    {
        size_t _n_deletions;
        size_t _n_type1; // undercall of one homopolymer
        size_t _n_type2; // undercall of two neighbor homopolymers

        uint[] _neighbor_intensities;
        ushort[][] _intensity_pairs; // neighboring undercalls
    }

    this(size_t maximum_intensity_value=1536)
    {
        _neighbor_intensities = new uint[maximum_intensity_value];
        _intensity_pairs = new ushort[][](maximum_intensity_value, maximum_intensity_value);
    }

    void updateStatistics(Deletion)(Deletion deletion)
    {
        ++_n_deletions;

        switch (deletion.number_of_homopolymers_in_deleted_sequence)
        {
            case 1:
                handleType1(deletion);
                break;
            case 2:
                handleType2(deletion);
                break;
            default:
                break;
        }
    }

    private void handleType1(Deletion)(Deletion deletion)
    {
        with (deletion)
        {
            if (previous_flow_call.base == bases[0])
            {
                _neighbor_intensities[previous_flow_call.intensity_value] += 1;
            }
            else if (!next_flow_call.isNull && next_flow_call.base == bases[0])
            {
                _neighbor_intensities[next_flow_call.intensity_value] += 1;
            }
        }
    }

    private void handleType2(Deletion)(Deletion deletion)
    {
        with (deletion)
        {
            if (next_flow_call.isNull) 
                return;

            if (previous_flow_call.base == bases.front && 
                next_flow_call.base == bases.back)
            {
                auto left = previous_flow_call.intensity_value;
                auto right = next_flow_call.intensity_value;

                _intensity_pairs[left][right] += 1;
            }
        }
    }

    void printOneSidedUndercallsReport(string filename)
    {
        auto _file = File(filename, "w+");

        _file.writeln("# intensity: signal intensity value for undercalled homopolymer");
        _file.writeln("# count: number of undercalls with this intensity value");

        _file.writeln("intensity\tcount");

        foreach (size_t intensity, count; _neighbor_intensities)
        {
            if (count > 0)
                _file.writeln(intensity, '\t', count);
        }
    }
}
