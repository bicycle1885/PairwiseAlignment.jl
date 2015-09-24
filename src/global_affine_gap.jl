# Gotoh's algorithm (Global)
# --------------------------

function global_affinegap_align{T}(a, b, affinegap::AffineGap{T})
    m = length(a)
    n = length(b)
    subst_matrix = affinegap.subst_matrix
    ge = affinegap.gap_extend_penalty
    go = affinegap.gap_open_penalty
    goe = go + ge
    H = Matrix{T}(m + 1, n + 1)
    E = Matrix{T}(m, n)
    F = Matrix{T}(m, n)
    # run dynamic programming column by column
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
                    H[i+1,j] - go
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
    return H, E, F
end

function traceback(a, b, H, E, F, affinegap)
    subst_matrix = affinegap.subst_matrix
    ge = affinegap.gap_extend_penalty
    goe = affinegap.gap_open_penalty + ge
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
        elseif i == 0 || (j ≥ 1 && H[i+1,j+1] == E[i,j] && ((j ≥ 2 && E[i,j] + ge == E[i,j-1]) || E[i,j] + goe == H[i+1,j]))
            # ←
            push!(a′, '-')
            push!(b′, b[j])
            j -= 1
        elseif j == 0 || (i ≥ 1 && H[i+1,j+1] == F[i,j] && ((i ≥ 2 && F[i,j] + ge == F[i-1,j]) || F[i,j] + goe == H[i,j+1]))
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
