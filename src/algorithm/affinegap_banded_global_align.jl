# Run dynamic programming only within a band.
# Formally, a[i] and b[j] can be matched only if -L ≤ j - i ≤ U.
function affinegap_banded_global_align{T}(a, b, L::Int, U::Int, subst_matrix::AbstractSubstitutionMatrix{T}, gap_open_penalty::T, gap_extend_penalty::T)
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
    # To avoid confusion, the indices are specified as following formats:
    #   * H[i-j+U+1,j+1]
    #   * E[i-j+U+1,j  ]  (j ≥ 1)
    #   * F[i-j+U+1,j  ]  (j ≥ 1)
    # Hope these expressions are optimized away.
    @inbounds begin
        H[0-0+U+1,0+1] = T(0)
        for i in 1:L
            H[i-0+U+1,0+1] = affinegap_score(i, go, ge)
        end
        # add gap_extend_penalty to avoid overflow for integers
        minimum = typemin(T) + ge
        for j in 1:n
            if j ≤ U
                H[0-j+U+1,j+1] = affinegap_score(j, go, ge)
            end
            # vertical bounds along the j-th column
            lo = max(1, j - U)
            hi = min(m, j + L)
            for i in lo:hi
                E[i-j+U+1,j] = if j == 1 || i == j + L
                    minimum
                else
                    max(
                        H[i-(j-1)+U+1,(j-1)+1] - goe,
                        E[i-(j-1)+U+1,(j-1)  ] - ge
                    )
                end
                F[i-j+U+1,j] = if i == 1 || j == i + U
                    minimum
                else
                    max(
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
    end
    # the best alignment score is H[m-n+U+1,n+1]
    return H, E, F, (m-n+U+1, n+1)
end

function traceback(a, b, H, E, F, L, U, subst_matrix, gap_open_penalty, gap_extend_penalty)
    m = length(a)
    n = length(b)
    ge = gap_extend_penalty
    goe = gap_open_penalty + ge
    L = min(L, m)
    U = min(U, n)
    # gap/character counts (reversed order)
    counts_a = [0, 0]
    counts_b = [0, 0]
    # if gap extension is selected in the previous traceback step, either
    # gap extension or gap open in that direction should be selected.
    gapext_a = false
    gapext_b = false
    i = m
    j = n
    while i ≥ 1 && j ≥ 1
        @assert !(gapext_a && gapext_b)
        if gapext_a
            if j ≥ 2 && E[i-j+U+1,j] == E[i-(j-1)+U+1,j-1] - ge
                @gapext a
            elseif E[i-j+U+1,j] == H[(i+1)-j+U+1,j] - goe
                @gapopen a
            end
        elseif gapext_b
            if i ≥ 2 && F[i-j+U+1,j] == F[(i-1)-j+U+1,j] - ge
                @gapext b
            elseif F[i-j+U+1,j] == H[(i-1)-j+U+1,j+1] - goe
                @gapopen b
            end
        elseif H[i-j+U+1,j+1] == H[(i-1)-(j-1)+U+1,(j-1)+1] + subst_matrix[a[i],b[j]]
            @match
        elseif H[i-j+U+1,j+1] == E[i-j+U+1,j] && j > i - L
            if E[i-j+U+1,j] == E[i-(j-1)+U+1,j-1] - ge
                @gapext a
            elseif E[i-j+U+1,j] == H[i-(j-1)+U+1,(j-1)+1] - goe
                @gapopen a
            end
        elseif H[i-j+U+1,j+1] == F[i-j+U+1,j] && j < i + U
            if F[i-j+U+1,j] == F[(i-1)-j+U+1,j] - ge
                @gapext b
            elseif F[i-j+U+1,j] == H[(i-1)-j+U+1,j+1] - goe
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

function isinband(i, j, L, U, a, b)
    return 0 ≤ i ≤ length(a) && 0 ≤ j ≤ length(b) && i - L ≤ j ≤ i + U
end
