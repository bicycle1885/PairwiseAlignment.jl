module PairwiseAlignment

export
    pairalign,
    # types
    GappedSequence,
    AlignmentResult,
    # scoring systems
    AbstractSubstitutionMatrix,
    SubstitutionMatrix,
    AffineGapScoreModel,
    CostModel,
    AbstractSubstitutionMatrix,
    # alignment types
    GlobalAlignment,
    LocalAlignment,
    SemiGlobalAlignment,
    # distances
    EditDistance,
    HammingDistance,
    LevenshteinDistance,
    # operations
    leading_chars,
    leading_gaps,
    trailing_chars,
    trailing_gaps,
    counts,
    reversed_counts,
    push_chars!,
    push_gaps!


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


include("gapped_sequence.jl")
include("model.jl")
include("result.jl")
include("pairalign.jl")
include("algorithm/common.jl")
include("algorithm/affinegap_global_align.jl")
include("algorithm/affinegap_local_align.jl")
include("algorithm/affinegap_semiglobal_align.jl")
include("algorithm/affinegap_banded_global_align.jl")
include("algorithm/edit_distance.jl")
include("algorithm/hamming_distance.jl")

end # module
