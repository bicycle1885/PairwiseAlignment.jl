module PairwiseAlignment

export
    pairalign,
    # scoring systems
    AffineGap,
    # alignment types
    GlobalAlignment,
    LocalAlignment

# alignment types
immutable GlobalAlignment end
immutable LocalAlignment end
immutable SemiGlobalAlignment end

include("align.jl")
include("result.jl")
include("global_affine_gap.jl")
include("local_affine_gap.jl")

end # module
