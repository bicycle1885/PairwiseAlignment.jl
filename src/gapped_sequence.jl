# Gapped Sequence
# ---------------

type GappedSequence{S}
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

function GappedSequence{S}(seq::S, startpos::Integer)
    return GappedSequence{S}(seq, startpos, [0, 0])
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

# count the number of chars/gaps
leading_chars(gseq::GappedSequence) = gseq.counts[1]
function leading_gaps(gseq::GappedSequence)
    return leading_chars(gseq) > 0 ? 0 : gseq.counts[2]
end
trailing_gaps(gseq::GappedSequence) = gseq.counts[end]
function trailing_chars(gseq::GappedSequence)
    return trailing_gaps(gseq) > 0 ? 0 : gseq.counts[end-1]
end

function counts(gseq::GappedSequence)
    return copy(gseq.counts)
end

function reversed_counts(gseq::GappedSequence)
    counts = copy(gseq.counts)
    # char/gap counts => gap/char counts
    reverse!(counts)
    if counts[1] == 0
        shift!(counts)
    else
        unshift!(counts, 0)
    end
    if counts[end] == 0
        pop!(counts)
    else
        push!(counts, 0)
    end
    return counts
end

function push_chars!(gseq::GappedSequence, n_chars::Integer)
    if gseq.counts[end] == 0
        gseq.counts[end-1] += n_chars
    else
        push!(gseq.counts, n_chars, 0)
    end
    gseq.len += n_chars
    return gseq
end

function push_gaps!(gseq::GappedSequence, n_gaps::Integer)
    gseq.counts[end] += n_gaps
    gseq.len += n_gaps
    return gseq
end

function Base.append!(gseq::GappedSequence, counts::Vector)
    @assert rem(length(counts), 2) == 0
    d = sum(counts)
    # TODO: check boundaries
    append!(gseq.counts, counts)
    gseq.len += d
    return gseq
end

function gapmap(gseq::GappedSequence)
    gmap = falses(length(gseq))
    i = gseq.startpos
    j = 0
    for k in 1:div(length(gseq.counts), 2)
        count = gseq.counts[2k-1]
        j += count
        count = gseq.counts[2k]
        while count > 0
            gmap[j+=1] = true
            count -= 1
        end
    end
    return gmap
end

const GapChar = '-'

function Base.show(io::IO, gseq::GappedSequence)
    i = gseq.startpos
    for k in 1:div(length(gseq.counts), 2)
        count = gseq.counts[2k-1]
        while count > 0
            print(io, convert(Char, gseq.seq[i]))
            i += 1
            count -= 1
        end
        count = gseq.counts[2k]
        while count > 0
            print(io, GapChar)
            count -= 1
        end
    end
end
