# Gapped Sequence
# ---------------

immutable GappedSequence{S}
    seq::S
    startpos::Int
    len::Int
    counts::Vector{Int}  # the alternative number of char and gap counts
    function GappedSequence(seq::S, startpos::Integer, counts::Vector)
        @assert 1 ≤ startpos
        @assert rem(length(counts), 2) == 0
        len = sum(counts)
        return new(seq, startpos, len, counts)
    end
end

function GappedSequence{S}(seq::S, counts::Vector)
    return GappedSequence{S}(seq, 1, counts)
end

function GappedSequence{S}(seq::S, startpos::Integer, counts::Vector)
    return GappedSequence{S}(seq, startpos, counts)
end

# sequence without gaps
Base.convert{S}(::Type{GappedSequence}, seq::S) = GappedSequence(seq, [length(seq), 0])

Base.length(gseq::GappedSequence) = gseq.len

const GapChar = '-'

function Base.show(io::IO, gseq::GappedSequence)
    i = gseq.startpos
    for j in 1:div(length(gseq.counts), 2)
        count = gseq.counts[2j-1]
        while count > 0
            print(io, convert(Char, gseq.seq[i]))
            i += 1
            count -= 1
        end
        count = gseq.counts[2j]
        while count > 0
            print(io, GapChar)
            count -= 1
        end
    end
end
