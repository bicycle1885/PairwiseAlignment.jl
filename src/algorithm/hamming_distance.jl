function hamming_distance{T}(::Type{T}, a, b)
    @assert length(a) == length(b)
    d = T(0)
    for i in 1:length(a)
        d += ifelse(a[i] == b[i], T(0), T(1))
    end
    return d
end
