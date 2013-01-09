module events.mismatch;

import bio.bam.iontorrent.flowcall;
import bio.core.base;

import std.typecons;

struct Mismatch(BaseInfo)
{
    Nullable!ReadFlowCall previous_flow_call;
    Nullable!ReadFlowCall next_flow_call;

    size_t reference_position; // 0-based
    size_t offset; // offset on the query, 0-based

    Base reference_base;
    BaseInfo base;

    bool shares_flowcall_with_left_neighbour() @property const
    {
        return !previous_flow_call.isNull && previous_flow_call == base.flow_call;
    }

    bool shares_flowcall_with_right_neighbour() @property const
    {
        return !next_flow_call.isNull && next_flow_call == base.flow_call;
    }
}

auto mismatchEvents(BaseInfo)(BaseInfo[] bases)
{
    static struct Result
    {
        private BaseInfo[] _bases;

        this(BaseInfo[] bases)
        {
            _bases = bases;
        }

        int opApply(scope int delegate(const ref Mismatch!BaseInfo mismatch) dg)
        {
            Mismatch!BaseInfo result = void;

            foreach (size_t offset, base; _bases)
            {
                if (base.cigar_operation.type == 'M' && base != base.reference_base)
                {
                    result.previous_flow_call.nullify();
                    result.next_flow_call.nullify();

                    if (offset > 0)
                        result.previous_flow_call = _bases[offset - 1].flow_call;

                    if (offset + 1 < _bases.length)
                        result.next_flow_call = _bases[offset + 1].flow_call;

                    result.base = base;

                    result.reference_position = base.position;
                    result.offset = offset;

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
