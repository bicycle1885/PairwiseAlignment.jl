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
    # gap/character counts (reversed order)
    counts_a = [0, 0]
    counts_b = [0, 0]
    i = length(a)
    j = length(b)
    while i ≥ 1 || j ≥ 1
        if i ≥ 1 && j ≥ 1 && D[i+1,j+1] == D[i,j] + subst_matrix[a[i],b[j]]
            gap_a = false
            gap_b = false
            i -= 1
            j -= 1
        elseif i == 0 || (j ≥ 1 && D[i+1,j+1] == D[i+1,j] + insertion_cost)
            gap_a = true
            gap_b = false
            j -= 1
        elseif j == 0 || (i ≥ 1 && D[i+1,j+1] == D[i,j+1] + deletion_cost)
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
