using PairwiseAlignment
using Base.Test

immutable SimpleSubstMatrix{T} <: AbstractMatrix{T}
    matching_score::T
    mismatching_score::T
end
Base.getindex(m::SimpleSubstMatrix, x, y) = ifelse(x == y, m.matching_score, m.mismatching_score)

let
    affine_gap = AffineGap(
        SimpleSubstMatrix(0, -6),
        gap_open_penalty=5,
        gap_extend_penalty=3
    )

    # equivalent sequences
    a = "A"
    b = "A"
    algn = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test algn.score == 0
    a = "ACGT"
    b = "ACGT"
    algn = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test algn.score == 0

    # substitution(s)
    a = "ACGT"
    b = "ATGT"
    algn = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test algn.score == -6
    a = "ACGT"
    b = "ATTT"
    algn = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test algn.score == -12

    # b is shorter than a
    a = "ACGT"
    b = "ACG"
    algn = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test algn.score == -5
    a = "ACGT"
    b = "AC"
    algn = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test algn.score == -8
end
