# global-local alignment
function affinegap_semiglobal_align{T}(a, b, submat::AbstractSubstitutionMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
    m = length(a)
    n = length(b)
    go = gap_open_penalty
    ge = gap_extend_penalty
    goe = go + ge
    trace = Matrix{Trace}(m + 1, n + 1)
    H = Vector{T}(m + 1)
    E = Vector{T}(m)
    # run dynamic programming column by column
    @inbounds begin
        H[1] = T(0)
        for i in 1:m
            H[i+1] = affinegap_score(i, go, ge)
            trace[i+1,1] = i == 1 ? TRACE_GAPOPEN_B : TRACE_GAPEXTD_B
        end
        best_score = H[m+1]
        best_score_column = 0
        for j in 1:n
            h_diag = H[1]
            H[1] = T(0)
            trace[1,j+1] = j == 1 ? TRACE_GAPOPEN_A : TRACE_GAPEXTD_A
            # any value goes well since this will be set in the first iteration
            F = T(0)
            for i in 1:m
                # gap in the sequence A
                e = H[i+1] - goe
                gap_a = TRACE_GAPOPEN_A
                if j > 1
                    e′ = E[i] - ge
                    if e′ == e
                        gap_a |= TRACE_GAPEXTD_A
                    elseif e′ > e
                        gap_a  = TRACE_GAPEXTD_A
                        e = e′
                    end
                end
                # gap in the sequence B
                f = H[i] - goe
                gap_b = TRACE_GAPOPEN_B
                if i > 1
                    f′ = F - ge
                    if f′ == f
                        gap_b |= TRACE_GAPEXTD_B
                    elseif f′ > f
                        gap_b  = TRACE_GAPEXTD_B
                        f = f′
                    end
                end
                # match
                h = h_diag + submat[a[i],b[j]]
                # find the best score and its trace
                best = max(e, f, h)
                t = TRACE_NONE
                e == best && (t |= gap_a)
                f == best && (t |= gap_b)
                h == best && (t |= TRACE_MATCH)
                # update
                E[i] = e
                F = f
                h_diag = H[i+1]
                H[i+1] = best
                trace[i+1,j+1] = t
            end
            if H[m+1] > best_score
                best_score = H[m+1]
                best_score_column = j
            end
        end
    end
    return best_score, trace, (m, best_score_column)
end


function semiglobal_traceback(a, b, trace, endpos)
    # gap/character counts (reversed order)
    counts_a = [0, 0]
    counts_b = [0, 0]
    # if gap extension is selected in the previous traceback step, either
    # gap extension or gap open in that direction should be selected.
    gapext_a = false
    gapext_b = false
    i, j = endpos
    while i ≥ 1
        @assert !(gapext_a && gapext_b)
        t = trace[i+1,j+1]
        if gapext_a
            if j ≥ 2 && t & TRACE_GAPEXTD_A > 0
                @gapext a
            elseif t & TRACE_GAPOPEN_A > 0
                @gapopen a
            end
        elseif gapext_b
            if i ≥ 2 && t & TRACE_GAPEXTD_B > 0
                @gapext b
            elseif t & TRACE_GAPOPEN_B > 0
                @gapopen b
            end
        elseif t & TRACE_MATCH > 0
            @match
        elseif j ≥ 2 && t & TRACE_GAPEXTD_A > 0
            @gapext a
        elseif t & TRACE_GAPOPEN_A > 0
            @gapopen a
        elseif i ≥ 2 && t & TRACE_GAPEXTD_B > 0
            @gapext b
        elseif t & TRACE_GAPOPEN_B > 0
            @gapopen b
        end
        # do not come here
        @assert false
    end
    reverse!(counts_a)
    reverse!(counts_b)
    return GappedSequence(a, 1, counts_a), GappedSequence(b, j + 1, counts_b)
end
