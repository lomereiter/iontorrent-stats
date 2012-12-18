module events.deletion;

import bio.bam.fz.flowcall;
import bio.core.sequence;

import std.typecons;
import std.algorithm;
import std.range;

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
}

auto deletionEvents(BaseInfo)(BaseInfo[] bases)
{
    static struct Result
    {
        private BaseInfo[] _bases;
        private bool _reversed;

        this(BaseInfo[] bases)
        {
            _reversed = false;
            _bases = bases;

            if (_bases.length > 1)
                _reversed = bases[$ - 1].position < bases[0].position;
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
                    auto n = result.bases.uniq().walkLength();
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

                    auto res = dg(result);
                    if (res != 0)
                        return res;
                }
            }

            return 0;
        }
    }

    return Result(bases);
}
