# Gotoh's algorithm
# -----------------

function global_affinegap_align{T}(a, b, affinegap::AffineGap{T})
    m = length(a)
    n = length(b)
    subst_matrix = affinegap.subst_matrix
    gap_open_penalty = affinegap.gap_open_penalty
    gap_extend_penalty = affinegap.gap_extend_penalty
    H = Matrix{T}(m + 1, n + 1)
    E = Matrix{T}(m, n)
    F = Matrix{T}(m, n)
    # run dynamic programming column by column
    H[1,1] = T(0)
    for i in 1:m
        H[i+1,1] = -gap_open_penalty - gap_extend_penalty * i
    end
    for j in 1:n
        H[1,j+1] = -gap_open_penalty - gap_extend_penalty * j
        for i in 1:m
            if j == 1
                E[i,j] = H[i+1,j] - gap_open_penalty
            else
                E[i,j] = max(
                    E[i,j-1] - gap_extend_penalty,
                    H[i+1,j] - gap_open_penalty
                )
            end
            if i == 1
                F[i,j] = H[i,j+1] - gap_open_penalty
            else
                F[i,j] = max(
                    F[i-1,j] - gap_extend_penalty,
                    H[i,j+1] - gap_open_penalty
                )
            end
            H[i+1,j+1] = max(
                E[i,j],
                F[i,j],
                H[i,j] + subst_matrix[a[i],b[j]]
            )
        end
    end
    return H, E, F
end

function traceback(a, b, H, E, F, affinegap)
    subst_matrix = affinegap.subst_matrix
    go = affinegap.gap_open_penalty
    ge = affinegap.gap_extend_penalty
    a′ = Char[]
    b′ = Char[]
    i = length(a)
    j = length(b)
    while i ≥ 1 || j ≥ 1
        if i ≥ 1 && j ≥ 1 && H[i+1,j+1] - subst_matrix[a[i],b[j]] == H[i,j]
            # ↖
            push!(a′, a[i])
            push!(b′, b[j])
            i -= 1
            j -= 1
        elseif i == 0 || (j ≥ 1 && H[i+1,j+1] == E[i,j] && ((j ≥ 2 && E[i,j] + ge == E[i,j-1]) || E[i,j] + go == H[i+1,j]))
            # ←
            push!(a′, '-')
            push!(b′, b[j])
            j -= 1
        elseif j == 0 || (i ≥ 1 && H[i+1,j+1] == F[i,j] && ((i ≥ 2 && F[i,j] + ge == F[i-1,j]) || F[i,j] + go == H[i,j+1]))
            # ↑
            push!(a′, a[i])
            push!(b′, '-')
            i -= 1
        else
            @assert false
        end
    end
    reverse!(a′)
    reverse!(b′)
    return ASCIIString(a′), ASCIIString(b′)
end
