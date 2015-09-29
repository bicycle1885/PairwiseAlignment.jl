# AlignmentResult
# ---------------

typealias PairedSequences{S1,S2} Pair{GappedSequence{S1},GappedSequence{S2}}

type AlignmentResult{T,S1,S2}
    score::T
    seqpair::Nullable{PairedSequences{S1,S2}}
end

AlignmentResult{S1,S2}(::Type{S1}, ::Type{S2}, score) = AlignmentResult(score, Nullable{PairedSequences{S1,S2}}())
function AlignmentResult(score, a′::GappedSequence, b′::GappedSequence)
    if length(a′) != length(b′)
        error("gapped sequences have different length")
    end
    AlignmentResult(score, Nullable(a′ => b′))
end

function Base.show(io::IO, aln::AlignmentResult)
    print(io, "Score: ", aln.score)
    if !isnull(aln.seqpair)
        pair = get(aln.seqpair)
        a′ = string(pair.first)
        b′ = string(pair.second)
        print(io, '\n')
        print(io, "  a: ", a′, '\n')
        print(io, "     ", matching_string(a′, b′), '\n')
        print(io, "  b: ", b′)
    end
end

function Base.getindex(aln::AlignmentResult, i::Integer)
    if isnull(aln.seqpair)
        error("alignment is undefined")
    end
    pair = get(aln.seqpair)
    if i == 1
        return pair.first
    elseif i == 2
        return pair.second
    end
    throw(BoundsError(i))
end

function matching_string(a, b)
    @assert length(a) == length(b)
    str = Char[]
    for i in 1:length(a)
        if a[i] == b[i] && !(a[i] == GapChar || b[i] == GapChar)
            push!(str, '|')
        else
            push!(str, ' ')
        end
    end
    return ASCIIString(str)
end
