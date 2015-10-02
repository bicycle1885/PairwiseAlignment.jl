# Interfaces
# ----------

function pairalign{S1,S2}(::GlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false,
                          banded::Bool=false, lower::Int=0, upper::Int=0,
                          #linear_space::Bool=score_only,
                          )
    subst_matrix = score.subst_matrix
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
            H, _, _, bestpos = affinegap_banded_global_align(a, b, L, U, subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(S1, S2, H[bestpos...])
        else
            H, E, F, bestpos = affinegap_banded_global_align(a, b, L, U, subst_matrix, gap_open_penalty, gap_extend_penalty)
            a′, b′ = traceback(a, b, H, E, F, L, U, subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(H[bestpos...], a′, b′)
        end
    else
        if score_only
            H, _, _ = affinegap_global_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(S1, S2, H[end,end])
        else
            H, E, F = affinegap_global_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
            a′, b′ = affinegap_global_traceback(a, b, H, E, F, (length(a), length(b)), subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(H[end,end], a′, b′)
        end
    end
    error("not implemented")
end

function pairalign{S1,S2}(::LocalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false,
                          #linear_space::Bool=score_only
                          )
    subst_matrix = score.subst_matrix
    gap_open_penalty = score.gap_open_penalty
    gap_extend_penalty = score.gap_extend_penalty
    if score_only
        H, _, _, best_endpos = affinegap_local_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
        return AlignmentResult(S1, S2, H[best_endpos[1]+1,best_endpos[2]+1])
    else
        H, E, F, best_endpos = affinegap_local_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
        a′, b′ = traceback(a, b, H, E, F, best_endpos, subst_matrix, gap_open_penalty, gap_extend_penalty)
        return AlignmentResult(H[best_endpos[1]+1,best_endpos[2]+1], a′, b′)
    end
    error("not implemented")
end

function pairalign{S1,S2}(::SemiGlobalAlignment, a::S1, b::S2, score::AffineGapScoreModel;
                          score_only::Bool=false)
    subst_matrix = score.subst_matrix
    gap_open_penalty = score.gap_open_penalty
    gap_extend_penalty = score.gap_extend_penalty
    if length(a) > length(b)
        error("the first sequence should be shorter or equal to the second sequence")
    end
    if score_only
        H, _, _, best_endpos = affinegap_semiglobal_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
        return AlignmentResult(S1, S2, H[best_endpos[1]+1,best_endpos[2]+1])
    else
        H, E, F, best_endpos = affinegap_semiglobal_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
        a′, b′ = semiglobal_traceback(a, b, H, E, F, best_endpos, subst_matrix, gap_open_penalty, gap_extend_penalty)
        return AlignmentResult(H[best_endpos[1]+1,best_endpos[2]+1], a′, b′)
    end
    error("not implemented")
end

function pairalign{S1,S2}(::EditDistance, a::S1, b::S2, cost::CostModel;
                          distance_only::Bool=false)
    subst_matrix = cost.subst_matrix
    insertion_cost = cost.insertion_cost
    deletion_cost = cost.deletion_cost
    if distance_only
        D = edit_distance(a, b, subst_matrix, insertion_cost, deletion_cost)
        return AlignmentResult(S1, S2, D[end,end])
    else
        D = edit_distance(a, b, subst_matrix, insertion_cost, deletion_cost)
        a′, b′ = traceback(a, b, D, subst_matrix, insertion_cost, deletion_cost)
        return AlignmentResult(D[end,end], a′, b′)
    end
end

function pairalign(::LevenshteinDistance, a, b;
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
