using PairwiseAlignment
using Base.Test

immutable SimpleSubstMatrix{T} <: AbstractSubstitutionMatrix{T}
    matching_score::T
    mismatching_score::T
end
Base.getindex(m::SimpleSubstMatrix, x, y) = ifelse(x == y, m.matching_score, m.mismatching_score)

rmspaces(s) = replace(s, ' ', "")

let
    # Global Alignment
    affine_gap = AffineGapScoreModel(
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
    affine_gap = AffineGapScoreModel(
        SimpleSubstMatrix(0, -6),
        gap_open_penalty=5,
        gap_extend_penalty=3
    )

    a = "A"
    b = "A"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true)
    @test aln.score == 0
    a = "ACGT"
    b = "ACGT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true)
    @test aln.score == 0

    # substitution(s)
    a = "ACGT"
    b = "ATGT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true)
    @test aln.score == -6
    a = "ACGT"
    b = "ATTT"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true)
    @test aln.score == -12

    # b is shorter than a
    a = "ACGT"
    b = "ACG"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap,
                    banded=true, lower=1, upper=1)
    @test aln.score == -8
    a = "ACGT"
    b = "AC"
    aln = pairalign(GlobalAlignment(), a, b, affine_gap,
                    banded=true, lower=2, upper=2)
    @test aln.score == -11
end

let
    affine_gap = AffineGapScoreModel(
        SimpleSubstMatrix(0, -6),
        gap_open_penalty=5,
        gap_extend_penalty=3
    )
    # same length
    for (a, b) in [("A", "A"), ("ACGT", "ACGT"), ("ACGT", "ACGG")]
        aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=false, score_only=false)
        aln2 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=false, score_only=true)
        aln3 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true,  score_only=false)
        aln4 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true,  score_only=true)
        @test aln1.score == aln2.score == aln3.score == aln4.score
    end
    for (a, b) in [("A", "AT"), ("A", "ATT"), ("ACGT", "ACGTT"), ("ACGT", "ACGGT"), ("ACGT", "ACCCGT")]
        # b is longer
        aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=false, score_only=false)
        aln2 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=false, score_only=true)
        aln3 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true,  lower=2, upper=2, score_only=false)
        aln4 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true,  lower=2, upper=2, score_only=true)
        @test aln1.score == aln2.score == aln3.score == aln4.score

        # a is longer
        a, b = b, a
        aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=false, score_only=false)
        aln2 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=false, score_only=true)
        aln3 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true,  lower=2, upper=2, score_only=false)
        aln4 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true,  lower=2, upper=2, score_only=true)
        @test aln1.score == aln2.score == aln3.score == aln4.score
    end
end

let
    # Local Alignment
    affine_gap = AffineGapScoreModel(
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

let
    # Hamming Distance
    a = "ACGT"
    b = "ACGT"
    aln = pairalign(HammingDistance(), a, b)
    @test aln.score == 0
    a = "ACGT"
    b = "ACCT"
    aln = pairalign(HammingDistance(), a, b)
    @test aln.score == 1
    a = "ACGT"
    b = "ACG"
    @test_throws ErrorException pairalign(HammingDistance(), a, b)
end

let
    # Levenshtein Distance
    a = "ACGT"
    b = "ACGT"
    aln = pairalign(LevenshteinDistance(), a, b)
    @test aln.score == 0
    a = "ACGT"
    b = "ACCT"
    aln = pairalign(LevenshteinDistance(), a, b)
    @test aln.score == 1
    a = "ACGT"
    b = "ACGGA"
    aln = pairalign(LevenshteinDistance(), a, b)
    @test aln.score == 2
    aln = pairalign(LevenshteinDistance(), b, a)
    @test aln.score == 2
end

immutable SimpleCostMatrix{T} <: AbstractSubstitutionMatrix{T}
    mismatching_cost::T
end

Base.getindex{T}(m::SimpleCostMatrix{T}, x, y) = ifelse(x == y, 0, m.mismatching_cost)

let
    # Edit Distance
    costmodel = CostModel(
        SimpleCostMatrix(3),
        insertion_cost=4,
        deletion_cost=5
    )

    # equivalent sequences
    a = "ACGT"
    b = "ACGT"
    aln = pairalign(EditDistance(), a, b, costmodel)
    @test aln.score == 0

    # substitution
    a = "ACGT"
    b = "ACCT"
    aln = pairalign(EditDistance(), a, b, costmodel)
    @test aln.score == 3

    # indel
    a = "ACG"
    b = "ACGT"
    aln = pairalign(EditDistance(), a, b, costmodel)
    @test aln.score == 4
    aln = pairalign(EditDistance(), b, a, costmodel)
    @test aln.score == 5
end
