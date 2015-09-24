# AlignmentResult
# ---------------

const AlignmentType = Pair{ASCIIString,ASCIIString}

type AlignmentResult{T}
    score::T
    alignment::Nullable{AlignmentType}
end

AlignmentResult(score) = AlignmentResult(score, Nullable{AlignmentType}())
AlignmentResult(score, a′, b′) = AlignmentResult(score, Nullable(a′ => b′))

function Base.show(io::IO, aln::AlignmentResult)
    print(io, "Score: ", aln.score)
    if !isnull(aln.alignment)
        a = get(aln.alignment)
        println(io)
        println(io, "  a: ", a.first)
        print(io,   "  b: ", a.second)
    end
    return
end
