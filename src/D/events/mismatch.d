module events.mismatch;

import bio.bam.iontorrent.flowcall;
import bio.core.base;

import std.typecons;

struct Mismatch(BaseInfo)
{
    Nullable!BaseInfo previous_base_info;
    Nullable!BaseInfo next_base_info;

    size_t reference_position; // 0-based
    size_t offset; // offset on the query, 0-based

    Base reference_base;
    BaseInfo base_info;

    bool shares_flowcall_with_left_neighbour() @property const
    {
        return !previous_base_info.isNull && 
               previous_base_info.flow_call == base_info.flow_call;
    }

    bool shares_flowcall_with_right_neighbour() @property const
    {
        return !next_base_info.isNull && 
               next_base_info.flow_call == base_info.flow_call;
    }

    bool starts_swap() @property const
    {
        if (next_base_info.isNull)
            return false;

        return next_base_info.base == reference_base &&
               next_base_info.reference_base == base_info.base;       
    }

    bool ends_swap() @property const
    {
        if (previous_base_info.isNull)
            return false;

        return previous_base_info.base == reference_base &&
               previous_base_info.reference_base == base_info.base;
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
                    result.previous_base_info.nullify();
                    result.next_base_info.nullify();

                    if (offset > 0)
                        result.previous_base_info = _bases[offset - 1];

                    if (offset + 1 < _bases.length)
                        result.next_base_info = _bases[offset + 1];

                    result.base_info = base;

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
