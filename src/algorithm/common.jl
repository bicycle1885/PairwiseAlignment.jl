# some utils shared among algorithms

# k: gap length
function affinegap_score(k, gap_open_penalty, gap_extend_penalty)
    return -(gap_open_penalty + gap_extend_penalty * k)
end


# trace type for pairwise alignment
typealias Trace UInt8

# trace bitmap
const TRACE_NONE      = 0b00000
const TRACE_MATCH     = 0b00001
const TRACE_GAPOPEN_A = 0b00010
const TRACE_GAPEXTD_A = 0b00100
const TRACE_GAPOPEN_B = 0b01000
const TRACE_GAPEXTD_B = 0b10000


# update counts of a gapped sequence
function update_counts!(counts, isgap)
    if isgap
        if counts[end] > 0
            push!(counts, 0, 0)
        end
        counts[end-1] += 1
    else
        counts[end] += 1
    end
    return counts
end

# traceback utils:
#   * match
#   * gap
#   * gapext
#   * gapopen

macro match()
    esc(quote
        update_counts!(counts_a, false)
        update_counts!(counts_b, false)
        i -= 1
        j -= 1
        gapext_a = false
        gapext_b = false
        continue
    end)
end

macro gap(seq)
    if seq == :a
        return esc(quote
            update_counts!(counts_a, true)
            update_counts!(counts_b, false)
            j -= 1
            continue
        end)
    elseif seq == :b
        return esc(quote
            update_counts!(counts_a, false)
            update_counts!(counts_b, true)
            i -= 1
            continue
        end)
    end
    error("unknown seq: $seq")
end

macro gapopen(seq)
    if seq == :a
        return esc(quote
            gapext_a = false
            @gap a
        end)
    elseif seq == :b
        return esc(quote
            gapext_b = false
            @gap b
        end)
    end
    error("unknown seq: $seq")
end

macro gapext(seq)
    if seq == :a
        return esc(quote
            gapext_a = true
            @gap a
        end)
    elseif seq == :b
        return esc(quote
            gapext_b = true
            @gap b
        end)
    end
    error("unknown seq: $seq")
end

