# Run dynamic programming only within a band.
# Formally, a[i] and b[j] can be matched only if -L ≤ j - i ≤ U.
function affinegap_banded_global_align{T}(a, b, L::Int, U::Int, submat::AbstractSubstitutionMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
    m = length(a)
    n = length(b)
    go = gap_open_penalty
    ge = gap_extend_penalty
    goe = go + ge
    L = min(L, m)
    U = min(U, n)
    # band width
    W = L + U + 1
    trace = Matrix{Trace}(W, n + 1)
    H = Vector{T}(W)
    E = Vector{T}(W)
    # In order to save the working space, the matrices are vertically sheared.
    # Namely, the coordinate is transformed as (i, j) → (i-j+U, j), where
    # (i, j) points to the subproblem of the optimal alignment between two
    # prefixes, a[1:i] and b[1:j], within the band.
    # The position (i, j) is accessed with:
    #   * trace[i-j+U+1,j+1]
    #   * H[i-j+U+1]
    #   * E[i-j+U+1]  (j ≥ 1)
    @inbounds begin
        # (i, j) = (0, 0)
        H[0-0+U+1] = T(0)
        trace[0-0+U+1,0+1] = TRACE_NONE
        for i in 1:L
            # j = 0
            H[i-0+U+1] = affinegap_score(i, go, ge)
            trace[i-0+U+1,0+1] = i == 1 ? TRACE_GAPOPEN_B : TRACE_GAPEXTD_B
        end
        # add gap_extend_penalty to avoid overflow for integers
        minimum = typemin(T) + ge
        for j in 1:n
            if j ≤ U
                # i = 0
                H[0-j+U+1] = affinegap_score(j, go, ge)
                trace[0-j+U+1,j+1] = j == 1 ? TRACE_GAPOPEN_A : TRACE_GAPEXTD_A
            end
            F = T(0)
            # vertical bounds along the j-th column
            lo = max(1, j - U)
            hi = min(m, j + L)
            for i in lo:hi
                # gap in the sequence A
                if j == 1 || i == j + L
                    e = minimum
                    gap_a = TRACE_NONE
                else
                    e = H[i-(j-1)+U+1] - goe
                    gap_a = TRACE_GAPOPEN_A
                    e′ = E[i-(j-1)+U+1] - ge
                    if e′ == e
                        gap_a |= TRACE_GAPEXTD_A
                    elseif e′ > e
                        gap_a  = TRACE_GAPEXTD_A
                        e = e′
                    end
                end
                # gap in the sequence B
                if i == 1 || j == i + U
                    f = minimum
                    gap_b = TRACE_NONE
                else
                    f = H[(i-1)-j+U+1] - goe
                    gap_b = TRACE_GAPOPEN_B
                    f′ = F - ge
                    if f′ == f
                        gap_b |= TRACE_GAPEXTD_B
                    elseif f′ > f
                        gap_b  = TRACE_GAPEXTD_B
                        f = f′
                    end
                end
                # match
                h = H[i-j+U+1] + submat[a[i],b[j]]
                # find the best score and its trace
                best = max(e, f, h)
                t = TRACE_NONE
                e == best && (t |= gap_a)
                f == best && (t |= gap_b)
                h == best && (t |= TRACE_MATCH)
                # update
                E[i-j+U+1] = e
                F = f
                H[i-j+U+1] = best
                trace[i-j+U+1,j+1] = t
            end
        end
    end
    return H[m-n+U+1], trace
end

function affinegap_banded_global_traceback(a, b, L, U, trace, endpos)
    L = min(L, length(a))
    U = min(U, length(b))
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
        t = trace[i-j+U+1,j+1]
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

function isinband(i, j, L, U, a, b)
    return 0 ≤ i ≤ length(a) && 0 ≤ j ≤ length(b) && i - L ≤ j ≤ i + U
end
