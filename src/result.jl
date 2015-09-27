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
        println(io, "     ", matching_string(a.first, a.second))
        print(io,   "  b: ", a.second)
    end
    return
end

function matching_string(a, b)
    @assert length(a) == length(b)
    str = Char[]
    for i in 1:length(a)
        if a[i] == b[i] && !(a[i] == '-' || b[i] == '-')
            push!(str, '|')
        else
            push!(str, ' ')
        end
    end
    return ASCIIString(str)
end
