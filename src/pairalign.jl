# Interfaces
# ----------

function pairalign{S1,S2}(::GlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false,
                          banded::Bool=false, lower::Int=0, upper::Int=0,
                          #linear_space::Bool=score_only,
                          )
    submat = score.submat
    gap_open_penalty = score.gap_open_penalty
    gap_extend_penalty = score.gap_extend_penalty
    if banded
        L = lower
        U = upper
        # check whether the starting and ending positions of the DP matrix are included in the band.
        if !isinband(0, 0, L, U, a, b)
            error("the starting position is not included in the band")
        elseif !isinband(length(a), length(b), L, U, a, b)
            error("the ending position is not included in the band")
        end
        if score_only
            score, _ = affinegap_banded_global_align(a, b, L, U, submat, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(S1, S2, score)
        else
            score, trace = affinegap_banded_global_align(a, b, L, U, submat, gap_open_penalty, gap_extend_penalty)
            a′, b′ = affinegap_banded_global_traceback(a, b, L, U, trace, (length(a), length(b)))
            return AlignmentResult(score, a′, b′)
        end
    else
        if score_only
            score, _ = affinegap_global_align(a, b, submat, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(S1, S2, score)
        else
            score, trace = affinegap_global_align(a, b, submat, gap_open_penalty, gap_extend_penalty)
            a′, b′ = affinegap_global_traceback(a, b, trace, (length(a), length(b)))
            return AlignmentResult(score, a′, b′)
        end
    end
    error("not implemented")
end

function pairalign{S1,S2}(::LocalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false,
                          #linear_space::Bool=score_only
                          )
    submat = score.submat
    gap_open_penalty = score.gap_open_penalty
    gap_extend_penalty = score.gap_extend_penalty
    if score_only
        score, _, _ = affinegap_local_align(a, b, submat, gap_open_penalty, gap_extend_penalty)
        return AlignmentResult(S1, S2, score)
    else
        score, trace, best_endpos = affinegap_local_align(a, b, submat, gap_open_penalty, gap_extend_penalty)
        a′, b′ = affine_local_traceback(a, b, trace, best_endpos)
        return AlignmentResult(score, a′, b′)
    end
    error("not implemented")
end

function pairalign{S1,S2}(::SemiGlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false)
    submat = score.submat
    gap_open_penalty = score.gap_open_penalty
    gap_extend_penalty = score.gap_extend_penalty
    if length(a) > length(b)
        error("the first sequence should be shorter or equal to the second sequence")
    end
    if score_only
        score, _, _ = affinegap_semiglobal_align(a, b, submat, gap_open_penalty, gap_extend_penalty)
        return AlignmentResult(S1, S2, score)
    else
        score, trace, endpos = affinegap_semiglobal_align(a, b, submat, gap_open_penalty, gap_extend_penalty)
        a′, b′ = semiglobal_traceback(a, b, trace, endpos)
        return AlignmentResult(score, a′, b′)
    end
    error("not implemented")
end

function pairalign{S1,S2}(::EditDistance, a::S1, b::S2, cost::CostModel;
                          distance_only::Bool=false)
    submat = cost.submat
    insertion_cost = cost.insertion_cost
    deletion_cost = cost.deletion_cost
    if distance_only
        cost, _ = edit_distance(a, b, submat, insertion_cost, deletion_cost)
        return AlignmentResult(S1, S2, cost)
    else
        cost, trace = edit_distance(a, b, submat, insertion_cost, deletion_cost)
        a′, b′ = edit_traceback(a, b, trace, (length(a), length(b)))
        return AlignmentResult(cost, a′, b′)
    end
end

function pairalign{S1,S2}(::LevenshteinDistance, a::S1, b::S2;
                   distance_only::Bool=false)
    unitcost = CostModel(
        UnitSubstitutionCost{Int}(),
        insertion_cost=1,
        deletion_cost=1
    )
    return pairalign(EditDistance(), a, b, unitcost, distance_only=distance_only)
end

function pairalign{S1,S2}(::HammingDistance, a::S1, b::S2;
                          distance_only::Bool=false)
    if length(a) != length(b)
        error("two sequences should be same length")
    end
    if distance_only
        return AlignmentResult(S1, S2, hamming_distance(Int, a, b))
    else
        aln = hamming_distance(Int, a, b)
        # no gaps
        a′ = GappedSequence(a)
        b′ = GappedSequence(b)
        return AlignmentResult(aln, a′, b′)
    end
end
