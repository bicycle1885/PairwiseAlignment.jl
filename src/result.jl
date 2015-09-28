# AlignmentResult
# ---------------

typealias PairedSequences{S1,S2} Pair{GappedSequence{S1},GappedSequence{S2}}

type AlignmentResult{T,S1,S2}
    score::T
    seqpair::Nullable{PairedSequences{S1,S2}}
end

AlignmentResult{S1,S2}(::Type{S1}, ::Type{S2}, score) = AlignmentResult(score, Nullable{PairedSequences{S1,S2}}())
AlignmentResult(score, a′::GappedSequence, b′::GappedSequence) = AlignmentResult(score, Nullable(a′ => b′))

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
