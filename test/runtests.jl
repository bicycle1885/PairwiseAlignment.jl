using PairwiseAlignment
using Base.Test

immutable SimpleSubstMatrix{T} <: AbstractMatrix{T}
    matching_score::T
    mismatching_score::T
end
Base.getindex(m::SimpleSubstMatrix, x, y) = ifelse(x == y, m.matching_score, m.mismatching_score)

rmspaces(s) = replace(s, ' ', "")

let
    # Global Alignment
    affine_gap = AffineGap(
        SimpleSubstMatrix(0, -6),
        gap_open_penalty=5,
        gap_extend_penalty=3
    )

    # equivalent sequences
    a = "A"
    b = "A"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test aln.score == 0
    a = "ACGT"
    b = "ACGT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test aln.score == 0

    # substitution(s)
    a = "ACGT"
    b = "ATGT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test aln.score == -6
    a = "ACGT"
    b = "ATTT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test aln.score == -12

    # b is shorter than a
    a = "ACGT"
    b = "ACG"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test aln.score == -8
    a = "ACGT"
    b = "AC"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap)
    @test aln.score == -11
end

let
    # Banded Global Alignment
    affine_gap = AffineGap(
        SimpleSubstMatrix(0, -6),
        gap_open_penalty=5,
        gap_extend_penalty=3
    )

    a = "A"
    b = "A"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap,
                    banded=true, score_only=true)
    @test aln.score == 0
    a = "ACGT"
    b = "ACGT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap,
                    banded=true, score_only=true)
    @test aln.score == 0

    # substitution(s)
    a = "ACGT"
    b = "ATGT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap,
                    banded=true, score_only=true)
    @test aln.score == -6
    a = "ACGT"
    b = "ATTT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap,
                    banded=true, score_only=true)
    @test aln.score == -12

    # b is shorter than a
    a = "ACGT"
    b = "ACG"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap,
                    banded=true, score_only=true,
                    lower=1, upper=1)
    @test aln.score == -8
    a = "ACGT"
    b = "AC"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap,
                    banded=true, score_only=true,
                    lower=2, upper=2)
    @test aln.score == -11
end

let
    # Local Alignment
    affine_gap = AffineGap(
        SimpleSubstMatrix(3, -1),
        gap_open_penalty=2,
        gap_extend_penalty=1
    )

    # equivalent sequences
    a = "A"
    b = "A"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test aln.score == 3
    a = "ACGT"
    b = "ACGT"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test aln.score == 12

    # substitution(s)
    a = "ACGT"
    b = "ATGT"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test aln.score == 3 * 3 - 1 * 1
    a = "ACGT"
    b = "ATTT"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test aln.score == 3 * 2 - 1 * 2

    # b is shorter than a
    a = "ACGT"
    b = "ACG"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test aln.score == 9
    a = "ACGT"
    b = "AC"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test aln.score == 6

    # partial approximate matching
    a =     "ACGTAG"
    b = "AGATACGTAGGGAGA"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test aln.score == 18
    a =     "ACG  GGG" |> rmspaces
    b = "AGATACGTAGGGAGA"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test aln.score == 3 * 6 - (2 + 1 + 1)
end
