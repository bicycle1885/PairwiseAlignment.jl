module PairwiseAlignment

export
    pairalign,
    # scoring systems
    AffineGap,
    CostModel,
    AbstractSubstitutionMatrix,
    # alignment types
    GlobalAlignment,
    LocalAlignment,
    # distances
    EditDistance,
    HammingDistance,
    LevenshteinDistance

# alignment types
# ---------------

# alignment
immutable GlobalAlignment end
immutable LocalAlignment end
immutable SemiGlobalAlignment end

# distances
immutable EditDistance end
immutable HammingDistance end
immutable LevenshteinDistance end


# the gap penalty of length `k` is gap_open_penalty + gap_extend_penalty * k
immutable AffineGap{T<:Real}
    subst_matrix::AbstractMatrix{T}
    gap_open_penalty::T
    gap_extend_penalty::T
end

function AffineGap{T}(subst_matrix::AbstractMatrix{T};
                      gap_open_penalty=nothing, gap_extend_penalty=nothing)
    if gap_open_penalty === nothing || gap_extend_penalty === nothing
        error("both gap_open_penalty and gap_extend_penalty should be set")
    end
    return AffineGap(subst_matrix, gap_open_penalty, gap_extend_penalty)
end


include("cost.jl")
include("result.jl")
include("pairalign.jl")
include("algorithm/affinegap_global_align.jl")
include("algorithm/affinegap_local_align.jl")
include("algorithm/affinegap_banded_global_align.jl")
include("algorithm/edit_distance.jl")
include("algorithm/hamming_distance.jl")

end # module
