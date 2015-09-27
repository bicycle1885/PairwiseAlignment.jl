# Gotoh's algorithm (Local)
# -------------------------

function affinegap_local_align{T}(a, b, subst_matrix::AbstractMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
    m = length(a)
    n = length(b)
    ge = gap_extend_penalty
    goe = gap_open_penalty + ge
    H = Matrix{T}(m + 1, n + 1)
    E = Matrix{T}(m, n)
    F = Matrix{T}(m, n)
    # run dynamic programming column by column
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
    return H, E, F, best_endpos
end

function traceback{T}(a, b, H, E, F, best_endpos, subst_matrix::AbstractMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
    ge = gap_extend_penalty
    goe = gap_open_penalty + ge
    a′ = Char[]
    b′ = Char[]
    i, j = best_endpos
    while H[i,j] > 0
        if i ≥ 1 && j ≥ 1 && H[i+1,j+1] == H[i,j] + subst_matrix[a[i],b[j]]
            # ↖
            push!(a′, a[i])
            push!(b′, b[j])
            i -= 1
            j -= 1
        elseif i == 0 || (j ≥ 1 && H[i+1,j+1] == E[i,j] && ((j ≥ 2 && E[i,j] == E[i,j-1] - ge) || E[i,j] == H[i+1,j] - goe))
            # ←
            push!(a′, '-')
            push!(b′, b[j])
            j -= 1
        elseif j == 0 || (i ≥ 1 && H[i+1,j+1] == F[i,j] && ((i ≥ 2 && F[i,j] == F[i-1,j] - ge) || F[i,j] == H[i,j+1] - goe))
            # ↑
            push!(a′, a[i])
            push!(b′, '-')
            i -= 1
        else
            @assert false
        end
    end
    push!(a′, a[i])
    push!(b′, b[j])
    reverse!(a′)
    reverse!(b′)
    return ASCIIString(a′), ASCIIString(b′)
end
