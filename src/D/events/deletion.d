module events.deletion;

import bio.bam.iontorrent.flowcall;
import bio.core.sequence;

import std.typecons;
import std.algorithm;
import std.range;

struct SmallArray(T, ubyte N) {
    private 
    {
        ubyte _len;
        union {
            T[N] _buf;
            T[] _array;
        }
    }

    this(size_t n)
    {
        if (n > ubyte.max)
            n = ubyte.max;

        _len = cast(ubyte)n;
        if (n > N)
            _array = new T[n];
    }

    T[] data() @property {
        return _len > N ? _array[] : _buf[0 .. _len];
    }

    size_t length() @property const {
        return _len > N ? _array.length : _len;
    }

    void length(size_t size) @property {
        assert(size < N); // FIXME

        if (_len > N)
            _array.length = size;
        else
            _len = cast(ubyte)size;
    }

    alias data this;
}

struct Deletion(BaseInfo)
{
    // It is assumed that deletion is surrounded by two successful base calls.
    ReadFlowCall previous_flow_call;
    ReadFlowCall next_flow_call;

    size_t start_position; // 0-based position of first deleted base on the reference
    size_t end_position; // 1-based position of last deleted base on the reference
    size_t read_offset; // after which read base the deletion occurred

    size_t length() @property const {
        return end_position - start_position;
    }

    size_t number_of_homopolymers_in_deleted_sequence;

    NucleotideSequence bases; // deleted bases


    private SmallArray!(short, 8) _deleted_base_intensities;

    SmallArray!(short, 8) deleted_base_intensities()
    {
        return _deleted_base_intensities;
    }
}

auto deletionEvents(BaseInfo)(BaseInfo[] bases, string flow_order, in short[] intensities)
{
    static struct Result
    {
        private BaseInfo[] _bases;
        private bool _reversed;
        private string _fo;
        private const(short[]) _ints;

        this(BaseInfo[] bases, string flow_order, in short[] intensities)
        {
            _reversed = false;
            _bases = bases;

            if (_bases.length > 1)
                _reversed = bases[$ - 1].position < bases[0].position;

            _fo = flow_order;
            _ints = intensities;
        }

        int opApply(scope int delegate(/*const*/ ref Deletion!BaseInfo deletion) dg)
        {
            Deletion!BaseInfo result = void;

            foreach (size_t offset, base; _bases)
            {
                if (base.cigar_operation_offset < base.cigar_operation.length - 1)
                    continue;

                if (base.next_md_operation.isNull)
                    continue;

                if (base.cigar_after.empty)
                    continue;

                if (base.cigar_after.front.type != 'D')
                    continue;

                if (base.next_md_operation.is_deletion)
                {
                    result.bases = base.next_md_operation.deletion;

                    auto unique_bases = result.bases.uniq();
                    auto n = unique_bases.walkLength();
                    assert(n > 0);
                    result.number_of_homopolymers_in_deleted_sequence = n;

                    result.read_offset = offset;
                    result.start_position = base.position;
                    auto length = base.next_md_operation.deletion.length;

                    if (_reversed)
                    {
                        result.start_position -= (length - 1);
                    }

                    result.end_position = base.position + length;

                    result.previous_flow_call = _bases[offset].flow_call;

                    if (offset < _bases.length)
                    {
                        result.next_flow_call = _bases[offset + 1].flow_call;
                    }
                    else
                    {
                        continue;
                    }

                    result._deleted_base_intensities = SmallArray!(short, 8)(n);
                    short[] deleted_base_intensities = result._deleted_base_intensities.data;

                    auto current_base = unique_bases.front;
                    size_t j = 0; // number of consumed homopolymers

                    bool has_zero_intensity_calls = false;
                    bool has_nonzero_intensity_calls = false;

                    for (size_t i = result.previous_flow_call.flow_index;
                                i <= result.next_flow_call.flow_index;
                                i++)
                    {
                        if (_fo[i] == current_base)
                        {
                            if (_ints[i] == 0)
                            {
                                has_zero_intensity_calls = true;
                            }
                            else if (!has_nonzero_intensity_calls)
                            {
                                has_nonzero_intensity_calls = true;
                            }

                            deleted_base_intensities[j++] = _ints[i];
                            if (j == n)
                            {
                                break;
                            }
                            unique_bases.popFront();
                            assert(!unique_bases.empty);
                            current_base = unique_bases.front;
                        }
                    }

                    // If all intensity calls are 0, don't take this deletion
                    // into account.
                    // This heuristic works well on experimental data, i.e.
                    // 0.00 intensity has almost the same frequency as 0.01.
                    if (has_zero_intensity_calls && !has_nonzero_intensity_calls)
                    {
                        j = 0;
                    }

                    result._deleted_base_intensities.length = j;

                    auto res = dg(result);
                    if (res != 0)
                        return res;
                }
            }

            return 0;
        }
    }

    return Result(bases, flow_order, intensities);
}
