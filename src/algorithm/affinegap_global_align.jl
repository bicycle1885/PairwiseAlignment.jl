function affinegap_global_align{T}(a, b, submat::AbstractSubstitutionMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
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
        end
        for j in 1:n
            H′ = H[1]
            H[1] = affinegap_score(j, go, ge)
            F = T(0)  # any value goes well since this will be initialized in the first iteration
            for i in 1:m
                # gap in A
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
                # gap in B
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
                h = H′ + submat[a[i],b[j]]
                best = max(e, f, h)
                t = 0x00
                if e == best
                    t |= gap_a
                end
                if f == best
                    t |= gap_b
                end
                if h == best
                    t |= TRACE_MATCH
                end
                trace[i+1,j+1] = t
                E[i] = e
                F = f
                H′ = H[i+1]
                H[i+1] = best
            end
        end
    end
    # return the best score and the trace
    return H[end], trace
end

function affinegap_global_traceback(a, b, trace, endpos)
    # gap/character counts (reversed order)
    counts_a = [0, 0]
    counts_b = [0, 0]
    # if gap extension is selected in the previous traceback step, either
    # gap extension or gap open in that direction should be selected.
    gapext_a = false
    gapext_b = false
    i, j = endpos
    while i ≥ 1 && j ≥ 1
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
    while j ≥ 1 @gap a end
    while i ≥ 1 @gap b end
    reverse!(counts_a)
    reverse!(counts_b)
    return GappedSequence(a, counts_a), GappedSequence(b, counts_b)
end
