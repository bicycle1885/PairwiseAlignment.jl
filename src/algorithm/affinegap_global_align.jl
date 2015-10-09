# Gotoh's algorithm (Global)
# --------------------------

function affinegap_global_align{T}(a, b, submat::AbstractSubstitutionMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
    m = length(a)
    n = length(b)
    go = gap_open_penalty
    ge = gap_extend_penalty
    goe = go + ge
    H = Matrix{T}(m + 1, n + 1)
    E = Matrix{T}(m, n)
    F = Matrix{T}(m, n)
    # run dynamic programming column by column
    @inbounds begin
        H[1,1] = T(0)
        for i in 1:m
            H[i+1,1] = affinegap_score(i, go, ge)
        end
        for j in 1:n
            H[1,j+1] = affinegap_score(j, go, ge)
            for i in 1:m
                E[i,j] = if j == 1
                    H[i+1,j] - goe
                else
                    max(
                        E[i,j-1] - ge,
                        H[i+1,j] - goe
                    )
                end
                F[i,j] = if i == 1
                    H[i,j+1] - goe
                else
                    max(
                        F[i-1,j] - ge,
                        H[i,j+1] - goe
                    )
                end
                H[i+1,j+1] = max(
                    E[i,j],
                    F[i,j],
                    H[i,j] + submat[a[i],b[j]]
                )
            end
        end
    end
    return H, E, F
end

function affinegap_global_traceback(a, b, H, E, F, endpos, submat, gap_open_penalty, gap_extend_penalty)
    ge = gap_extend_penalty
    goe = gap_open_penalty + ge
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
        if gapext_a
            if j ≥ 2 && E[i,j] == E[i,j-1] - ge
                @gapext a
            elseif E[i,j] == H[i+1,j] - goe
                @gapopen a
            end
        elseif gapext_b
            if i ≥ 2 && F[i,j] == F[i-1,j] - ge
                @gapext b
            elseif F[i,j] == H[i,j+1] - goe
                @gapopen b
            end
        elseif H[i+1,j+1] == H[i,j] + submat[a[i],b[j]]
            @match
        elseif H[i+1,j+1] == E[i,j]
            # gap in a
            if j ≥ 2 && E[i,j] == E[i,j-1] - ge
                @gapext a
            elseif E[i,j] == H[i+1,j] - goe
                @gapopen a
            end
        elseif H[i+1,j+1] == F[i,j]
            # gap in b
            if i ≥ 2 && F[i,j] == F[i-1,j] - ge
                @gapext b
            elseif F[i,j] == H[i,j+1] - goe
                @gapopen b
            end
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
