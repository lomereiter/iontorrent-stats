module events.insertion;

import bio.bam.iontorrent.flowcall;

import std.typecons;
import std.algorithm;
import std.range;

struct Insertion(BaseInfo)
{
    Nullable!ReadFlowCall previous_flow_call;
    Nullable!ReadFlowCall next_flow_call;

    size_t start_offset;
    size_t end_offset;
    size_t reference_position;

    size_t length() @property const {
        return end_offset - start_offset;
    }

    size_t number_of_flowcalls;

    BaseInfo[] bases;

    auto flowcalls() @property const {
        return bases.map!"a.flow_call"().uniq();
    }

    ReadFlowCall first_flow_call() @property const {
        return bases[0].flow_call;
    }

    ReadFlowCall last_flow_call() @property const {
        return bases[$ - 1].flow_call;
    }
}

auto insertionEvents(BaseInfo)(BaseInfo[] bases)
{
    static struct Result
    {
        private BaseInfo[] _bases;

        this(BaseInfo[] bases)
        {
            _bases = bases;
        }

        int opApply(scope int delegate(const ref Insertion!BaseInfo insertion) dg)
        {
            Insertion!BaseInfo result = void;
            
            foreach (size_t offset, baseinfo; _bases)
            {
                if (baseinfo.cigar_operation.type == 'I' &&
                    baseinfo.cigar_operation_offset == 0)
                {
                    result.start_offset = offset;
                    result.end_offset = offset + baseinfo.cigar_operation.length;

                    assert(result.end_offset <= _bases.length);

                    result.bases = _bases[result.start_offset .. result.end_offset];

                    result.number_of_flowcalls = result.bases.map!"a.flow_call".uniq().walkLength();

                    result.reference_position = baseinfo.position;

                    result.previous_flow_call.nullify();
                    result.next_flow_call.nullify();

                    if (result.start_offset > 0)
                        result.previous_flow_call = _bases[result.start_offset - 1].flow_call;
                        
                    if (result.end_offset < _bases.length)
                        result.next_flow_call = _bases[result.end_offset].flow_call;

                    auto res = dg(result);

                    if (res != 0) return res;
                }        
            }

            return 0;
        }
    }

    return Result(bases);
}
