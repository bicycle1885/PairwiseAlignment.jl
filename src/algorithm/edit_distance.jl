function edit_distance{T}(a, b, submat::AbstractSubstitutionMatrix{T}, insertion_cost::T, deletion_cost::T)
    m = length(a)
    n = length(b)
    trace = Matrix{Trace}(m + 1, n + 1)
    D = Vector{T}(m + 1)
    D[1] = T(0)
    for i in 1:m
        D[i+1] = i * deletion_cost
        trace[i+1,1] = TRACE_INSERT
    end
    for j in 1:n
        d_diag = D[1]
        D[1] = j * insertion_cost
        trace[1,j+1] = TRACE_DELETE
        for i in 1:m
            del = D[i]   + deletion_cost
            ins = D[i+1] + insertion_cost
            mat = d_diag + submat[a[i],b[j]]
            # find the best score and its trace
            best = min(del, ins, mat)
            t = TRACE_NONE
            del == best && (t |= TRACE_DELETE)
            ins == best && (t |= TRACE_INSERT)
            mat == best && (t |= TRACE_MATCH)
            d_diag = D[i+1]
            D[i+1] = best
            trace[i+1,j+1] = t
        end
    end
    return D[m+1], trace
end

function edit_traceback(a, b, trace, endpos)
    # gap/character counts (reversed order)
    counts_a = [0, 0]
    counts_b = [0, 0]
    i, j = endpos
    while i ≥ 1 && j ≥ 2
        t = trace[i+1,j+1]
        if t & TRACE_MATCH > 0
            @match
        elseif t & TRACE_DELETE > 0
            @gap a
        elseif t & TRACE_INSERT > 0
            @gap b
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
