# Gotoh's algorithm (Global)
# --------------------------

function affinegap_global_align{T}(a, b, subst_matrix::AbstractSubstitutionMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
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
            H[i+1,1] = -go - ge * i
        end
        for j in 1:n
            H[1,j+1] = -go - ge * j
            for i in 1:m
                if j == 1
                    E[i,j] = H[i+1,j] - goe
                else
                    E[i,j] = max(
                        E[i,j-1] - ge,
                        H[i+1,j] - goe
                    )
                end
                if i == 1
                    F[i,j] = H[i,j+1] - goe
                else
                    F[i,j] = max(
                        F[i-1,j] - ge,
                        H[i,j+1] - goe
                    )
                end
                H[i+1,j+1] = max(
                    E[i,j],
                    F[i,j],
                    H[i,j] + subst_matrix[a[i],b[j]]
                )
            end
        end
    end
    return H, E, F
end

function traceback(a, b, H, E, F, subst_matrix, gap_open_penalty, gap_extend_penalty)
    ge = gap_extend_penalty
    goe = gap_open_penalty + ge
    # gap/character counts (reversed order)
    counts_a = [0, 0]
    counts_b = [0, 0]
    i = length(a)
    j = length(b)
    while i ≥ 1 || j ≥ 1
        if i ≥ 1 && j ≥ 1 && H[i+1,j+1] == H[i,j] + subst_matrix[a[i],b[j]]
            # ↖
            # ...a[i]...
            #    |
            # ...b[j]...
            gap_a = false
            gap_b = false
            i -= 1
            j -= 1
        elseif i == 0 || (j ≥ 1 && H[i+1,j+1] == E[i,j] && ((j ≥ 2 && E[i,j] == E[i,j-1] - ge) || E[i,j] == H[i+1,j] - goe))
            # ←
            # ... -    ...
            #     |
            # ... b[j] ...
            gap_a = true
            gap_b = false
            j -= 1
        elseif j == 0 || (i ≥ 1 && H[i+1,j+1] == F[i,j] && ((i ≥ 2 && F[i,j] == F[i-1,j] - ge) || F[i,j] == H[i,j+1] - goe))
            # ↑
            # ... a[i] ...
            #     |
            # ... -    ...
            gap_a = false
            gap_b = true
            i -= 1
        else
            @assert false
        end
        # update counts
        update_counts!(counts_a, gap_a)
        update_counts!(counts_b, gap_b)
    end
    reverse!(counts_a)
    reverse!(counts_b)
    return GappedSequence(a, counts_a), GappedSequence(b, counts_b)
end
