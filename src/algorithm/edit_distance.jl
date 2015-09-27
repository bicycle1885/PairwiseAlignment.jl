# Edit distance
# -------------

function edit_distance{T}(a, b, subst_matrix::AbstractSubstitutionMatrix{T}, insertion_cost::T, deletion_cost::T)
    m = length(a)
    n = length(b)
    D = Matrix{T}(m + 1, n + 1)
    D[1,1] = T(0)
    for i in 1:m
        D[i+1,1] = i * deletion_cost
    end
    for j in 1:n
        D[1,j+1] = j * insertion_cost
        for i in 1:m
            D[i+1,j+1] = min(
                D[i,  j+1] + deletion_cost,
                D[i+1,j  ] + insertion_cost,
                D[i  ,j  ] + subst_matrix[a[i],b[j]]
            )
        end
    end
    return D
end

function traceback(a, b, D, subst_matrix, insertion_cost, deletion_cost)
    a′ = Char[]
    b′ = Char[]
    i = length(a)
    j = length(b)
    while i ≥ 1 || j ≥ 1
        if i ≥ 1 && j ≥ 1 && D[i+1,j+1] == D[i,j] + subst_matrix[a[i],b[j]]
            push!(a′, a[i])
            push!(b′, b[j])
            i -= 1
            j -= 1
        elseif i == 0 || (j ≥ 1 && D[i+1,j+1] == D[i+1,j] + insertion_cost)
            push!(a′, '-')
            push!(b′, b[j])
            j -= 1
        elseif j == 0 || (i ≥ 1 && D[i+1,j+1] == D[i,j+1] + deletion_cost)
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
