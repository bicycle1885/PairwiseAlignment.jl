# Gotoh's algorithm (Local)
# -------------------------

function affinegap_local_align{T}(a, b, subst_matrix::AbstractSubstitutionMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
    m = length(a)
    n = length(b)
    ge = gap_extend_penalty
    goe = gap_open_penalty + ge
    H = Matrix{T}(m + 1, n + 1)
    E = Matrix{T}(m, n)
    F = Matrix{T}(m, n)
    # run dynamic programming column by column
    @inbounds begin
        H[1,1] = T(0)
        for i in 1:m
            H[i+1,1] = T(0)
        end
        best_score = typemin(T)
        best_endpos = (0, 0)
        for j in 1:n
            H[1,j+1] = T(0)
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
                    T(0),
                    E[i,j],
                    F[i,j],
                    H[i,j] + subst_matrix[a[i],b[j]]
                )
                if H[i+1,j+1] > best_score
                    best_score = H[i+1,j+1]
                    best_endpos = (i, j)
                end
            end
        end
    end
    return H, E, F, best_endpos
end

function traceback(a, b, H, E, F, best_endpos, subst_matrix, gap_open_penalty, gap_extend_penalty)
    ge = gap_extend_penalty
    goe = gap_open_penalty + ge
    # gap/character counts (reversed order)
    counts_a = [0, 0]
    counts_b = [0, 0]
    i, j = best_endpos
    while H[i,j] > 0
        if i ≥ 1 && j ≥ 1 && H[i+1,j+1] == H[i,j] + subst_matrix[a[i],b[j]]
            # ↖
            gap_a = false
            gap_b = false
            i -= 1
            j -= 1
        elseif i == 0 || (j ≥ 1 && H[i+1,j+1] == E[i,j] && ((j ≥ 2 && E[i,j] == E[i,j-1] - ge) || E[i,j] == H[i+1,j] - goe))
            # ←
            gap_a = true
            gap_b = false
            j -= 1
        elseif j == 0 || (i ≥ 1 && H[i+1,j+1] == F[i,j] && ((i ≥ 2 && F[i,j] == F[i-1,j] - ge) || F[i,j] == H[i,j+1] - goe))
            # ↑
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
    # update counts
    update_counts!(counts_a, false)
    update_counts!(counts_b, false)
    reverse!(counts_a)
    reverse!(counts_b)
    return GappedSequence(a, i, counts_a), GappedSequence(b, j, counts_b)
end
