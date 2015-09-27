# Run dynamic programming only within a band.
# Formally, a[i] and b[j] can be matched only if -L ≤ j - i ≤ U.
function affinegap_banded_global_align{T}(a, b, L::Int, U::Int, subst_matrix::AbstractMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
    m = length(a)
    n = length(b)
    go = gap_open_penalty
    ge = gap_extend_penalty
    goe = go + ge
    L = min(L, m)
    U = min(U, n)
    # band width
    W = L + U + 1
    # matrices for dynamic programming
    H = Matrix{T}(W, n + 1)
    E = Matrix{T}(W, n)
    F = Matrix{T}(W, n)
    # In order to save the working space, the matrices are vertically sheared.
    # Namely, the coordinate is transformed as (i, j) → (i-j+U, j), where
    # H[i+1,j+1] stores the best score of two prefixes: a[1:i] and b[1:j].
    # To avoid confusion, the matrices are specified as following formats:
    #   * H[i-j+U+1,j+1]
    #   * E[i-j+U+1,j  ]  (j ≥ 1)
    #   * F[i-j+U+1,j  ]  (j ≥ 1)
    # Hope these expressions are optimized away.
    H[0-0+U+1,0+1] = T(0)
    for i in 1:L
        H[i-0+U+1,0+1] = -(go + ge * i)
    end
    for j in 1:n
        if j ≤ U
            H[0-j+U+1,j+1] = -(go + ge * j)
        end
        # vertical bounds along the j-th column
        lo = max(1, j - U)
        hi = min(m, j + L)
        for i in lo:hi
            if j == 1 || i == j + L
                # add gap_extend_penalty to avoid overflow
                E[i-j+U+1,j] = typemin(T) + ge
            else
                E[i-j+U+1,j] = max(
                    H[i-(j-1)+U+1,(j-1)+1] - goe,
                    E[i-(j-1)+U+1,(j-1)  ] - ge
                )
            end
            if i == 1 || j == i + U
                # add gap_extend_penalty to avoid overflow
                F[i-j+U+1,j] = typemin(T) + ge
            else
                F[i-j+U+1,j] = max(
                    H[(i-1)-j+U+1,j+1] - goe,
                    F[(i-1)-j+U+1,j  ] - ge
                )
            end
            H[i-j+U+1,j+1] = max(
                E[i-j+U+1,j],
                F[i-j+U+1,j],
                H[(i-1)-(j-1)+U+1,(j-1)+1] + subst_matrix[a[i],b[j]]
            )
        end
    end
    # the best alignment score is H[m-n+U+1,n+1]
    return H, E, F, (m-n+U+1, n+1)
end

function traceback(a, b, H, E, F, L, U, subst_matrix, gap_open_penalty, gap_extend_penalty)
    m = length(a)
    n = length(b)
    ge = gap_extend_penalty
    goe = gap_open_penalty
    L = min(L, m)
    U = min(U, n)
    i = m
    j = n
    while i ≥ 1 || j ≥ 1
        i′ = i - j + U
        if i ≥ 1 && j ≥ 1 && H[(i-1)-(j-1)+U+1,j] == H[i′+1,j+1] - subst_matrix[a[i],b[j]]
            # ↖
            push!(a′, a[i])
            push!(b′, b[j])
            i -= 1
            j -= 1
        elseif i == 0
            # ←
            push!(a′, '-')
            push!(b′, b[j])
            j -= 1
        elseif j == 0
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
