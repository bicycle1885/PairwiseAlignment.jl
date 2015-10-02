using PairwiseAlignment
using Base.Test

let
    seq = "ACGT"
    gseq = GappedSequence(seq)
    @test isa(gseq, GappedSequence{ASCIIString})
    @test length(gseq) == 4

    gseq = GappedSequence(seq, 2)
    @test length(gseq) == 0
    @test string(gseq) == ""
    @test n_chars(gseq) == 0
    @test n_gaps(gseq) == 0
    @test leading_chars(gseq) == leading_gaps(gseq) == 0
    @test trailing_chars(gseq) == trailing_gaps(gseq) == 0
    @test gapmap(gseq) == BitVector([])
    push_chars!(gseq, 1)
    @test length(gseq) == 1
    @test string(gseq) == "C"
    @test n_chars(gseq) == 1
    @test n_gaps(gseq) == 0
    @test leading_chars(gseq) == trailing_chars(gseq) == 1
    @test leading_gaps(gseq) == trailing_gaps(gseq) == 0
    @test gapmap(gseq) == BitVector([false])
    push_gaps!(gseq, 1)
    @test length(gseq) == 2
    @test string(gseq) == "C-"
    @test n_chars(gseq) == 1
    @test n_gaps(gseq) == 1
    @test leading_chars(gseq) == 1
    @test leading_gaps(gseq) == 0
    @test trailing_chars(gseq) == 0
    @test trailing_gaps(gseq) == 1
    @test gapmap(gseq) == BitVector([false, true])

    counts = [2, 2, 2, 0]
    gseq = GappedSequence(seq, counts)
    @test length(gseq) == 6
    @test string(gseq) == "AC--GT"
    @test leading_chars(gseq) == trailing_chars(gseq) == 2
    @test leading_gaps(gseq) == trailing_gaps(gseq) == 0
    @test PairwiseAlignment.counts(gseq) == counts
    @test gapmap(gseq) == BitVector([false, false, true, true, false, false])

    counts = [2, 2]
    gseq = GappedSequence(seq, counts)
    @test length(gseq) == 4
    @test string(gseq) == "AC--"
    @test n_chars(gseq) == 2
    @test n_gaps(gseq) == 2
    @test leading_chars(gseq) == 2
    @test leading_gaps(gseq) == 0
    @test trailing_chars(gseq) == 0
    @test trailing_gaps(gseq) == 2
    @test PairwiseAlignment.counts(gseq) == counts

    counts = [0, 2, 2, 2, 2, 2]
    gseq = GappedSequence(seq, counts)
    @test length(gseq) == 10
    @test string(gseq) == "--AC--GT--"
    @test n_chars(gseq) == 4
    @test n_gaps(gseq) == 6
    @test leading_chars(gseq) == trailing_chars(gseq) == 0
    @test leading_gaps(gseq) == leading_gaps(gseq) == 2
    @test PairwiseAlignment.counts(gseq) == counts

    counts = [1, 1, 1, 1]
    gseq = GappedSequence(seq, 3, counts)
    @test length(gseq) == 4
    @test string(gseq) == "G-T-"

    counts = [2, 0]
    gseq = GappedSequence(seq, counts)
    @test length(gseq) == 2
    @test string(gseq) == "AC"
    push_chars!(gseq, 1)
    @test length(gseq) == 3
    @test string(gseq) == "ACG"
    append!(gseq, [1, 2])
    @test length(gseq) == 6
    @test string(gseq) == "ACGT--"
end

let
    # a substitution matrix
    subst_matrix = SubstitutionMatrix(eye(5, 5))
    @test subst_matrix[0,0] == 1.0
    @test subst_matrix[0,1] == 0.0
    subst_matrix[0,1] = -1.0
    @test subst_matrix[0,1] == -1.0
    # models
    affine_gap = AffineGapScoreModel(subst_matrix, gap_open_penalty=2, gap_extend_penalty=1)
    @test isa(affine_gap, AffineGapScoreModel{Float64})
    cost = CostModel(subst_matrix, insertion_cost=1, deletion_cost=2)
    @test isa(cost, CostModel{Float64})
end

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
        gap_extend_penalty=2
    )

    a = ""
    b = "A"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -7

    a = ""
    b = "AC"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -9

    a = "A"
    b = "A"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == 0

    a = "ACGT"
    b = "ACGT"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == 0

    a = "ACGT"
    b = "ATGT"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -6

    a = "ACGT"
    b = "TCGA"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -12

    a = "ACGT"
    b = "ATTT"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -12

    a = "ACGT"
    b = "ACG"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -7

    a = "ACGT"
    b = "AC"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -9

    a = "ACGT"
    b = "A  T" |> rmspaces
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -9

    a = "GC   T" |> rmspaces
    b = "GCACGT"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -11
    @test string(aln1[1]) == "GC---T"
    @test string(aln1[2]) == "GCACGT"

    a = "CGTATCCAACGG     TTGTGGGA " |> rmspaces
    #     ||||| || ||      |||| ||
    b = " GTATCAAATGGACCCACTGTGTGAG" |> rmspaces
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap)
    @test aln1.score == aln2.score == -7 + -6 + -6 + -15 + -6 + -6 + -7
end

let
    # Banded Global Alignment
    affine_gap = AffineGapScoreModel(
        SimpleSubstMatrix(0, -6),
        gap_open_penalty=5,
        gap_extend_penalty=2
    )

    a = ""
    b = "A"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=1, upper=1)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=1, upper=1)
    @test aln1.score == aln2.score == -7

    a = ""
    b = "AC"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=2, upper=2)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=2, upper=2)
    @test aln1.score == aln2.score == -9

    a = "A"
    b = "A"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=1, upper=1)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=1, upper=1)
    @test aln1.score == aln2.score == 0

    a = "ACGT"
    b = "ACGT"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=1, upper=1)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=1, upper=1)
    @test aln1.score == aln2.score == 0

    a = "ACGT"
    b = "ATGT"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=1, upper=1)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=1, upper=1)
    @test aln1.score == aln2.score == -6

    a = "ACGT"
    b = "TCGA"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=1, upper=1)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=1, upper=1)
    @test aln1.score == aln2.score == -12

    a = "ACGT"
    b = "ATTT"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=1, upper=1)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=1, upper=1)
    @test aln1.score == aln2.score == -12

    a = "ACGT"
    b = "ACG"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=1, upper=1)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=1, upper=1)
    @test aln1.score == aln2.score == -7

    a = "ACGT"
    b = "AC"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=2, upper=2)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=2, upper=2)
    @test aln1.score == aln2.score == -9

    a = "ACGT"
    b = "A  T" |> rmspaces
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=2, upper=2)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=2, upper=2)
    @test aln1.score == aln2.score == -9

    a = "GC   T" |> rmspaces
    b = "GCACGT"
    aln1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=3, upper=3)
    aln2 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=3, upper=3)
    @test aln1.score == aln2.score == -11
    @test string(aln1[1]) == "GC---T"
    @test string(aln1[2]) == "GCACGT"

    a = "TATGGCCCG TAG   TTT" |> rmspaces
    b = "TA    CCGATAGGGGTTA" |> rmspaces
    aln1_1 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=1, upper=1)
    aln2_1 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=1, upper=1)
    aln1_3 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=3, upper=3)
    aln2_3 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=3, upper=3)
    aln1_5 = pairalign(GlobalAlignment(), a, b, affine_gap, banded=true, lower=5, upper=5)
    aln2_5 = pairalign(GlobalAlignment(), b, a, affine_gap, banded=true, lower=5, upper=5)
    #@show aln1_1 aln1_3 aln1_5
    @test aln1_1.score == aln2_1.score < aln1_3.score == aln1_3.score < aln1_5.score == aln2_5.score == -13 + -7 + -11 + -6
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

    a =   "TTAGC   TAGT" |> rmspaces
    b = "GATTAGCACGTAGTA"
    aln = pairalign(LocalAlignment(), a, b, affine_gap)
    @test string(aln[1]) == "TTAGC---TAGT"
    @test string(aln[2]) == "TTAGCACGTAGT"
    @test aln.score == 3 * 5 + (-2 + -1 * 3) + 3 * 4

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
    # Semi-global Alignment
    affinegap = AffineGapScoreModel(
        SimpleSubstMatrix(0, -6),
        gap_open_penalty=5,
        gap_extend_penalty=2
    )

    a = "A"
    b = "A"
    aln = pairalign(SemiGlobalAlignment(), a, b, affinegap)
    @test aln.score == 0

    a =        "ACGT"
    b = "GATHAGCACGTAGTA"
    aln = pairalign(SemiGlobalAlignment(), a, b, affinegap)
    @test aln.score == 0

    a =   "TTAGC   TAGT" |> rmspaces
    b = "GATTAGCACGTAGTA"
    aln = pairalign(SemiGlobalAlignment(), a, b, affinegap)
    @test string(aln[1]) == "TTAGC---TAGT"
    @test string(aln[2]) == "TTAGCACGTAGT"
    @test aln.score == -11
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
