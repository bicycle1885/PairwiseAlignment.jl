# Interfaces
# ----------

function pairalign(::GlobalAlignment, a, b, score::AffineGap;
                   score_only::Bool=false,
                   #linear_space::Bool=score_only,
                   banded::Bool=false, lower::Int=0, upper::Int=0)
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
            return AlignmentResult(H[bestpos...])
        else
            H, E, F, bestpos = affinegap_banded_global_align(a, b, L, U, subst_matrix, gap_open_penalty, gap_extend_penalty)
            a′, b′ = traceback(a, b, H, E, F, L, U, subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(H[bestpos...], a′, b′)
        end
    else
        if score_only
            H, _, _ = affinegap_global_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(H[end,end])
        else
            H, E, F = affinegap_global_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
            a′, b′ = traceback(a, b, H, E, F, subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(H[end,end], a′, b′)
        end
    end
    error("not implemented")
end

function pairalign(::LocalAlignment, a, b, score::AffineGap;
                   score_only::Bool=false,
                   #linear_space::Bool=score_only
                   )
    subst_matrix = score.subst_matrix
    gap_open_penalty = score.gap_open_penalty
    gap_extend_penalty = score.gap_extend_penalty
    if score_only
        H, _, _, best_endpos = affinegap_local_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
        return AlignmentResult(H[best_endpos[1]+1,best_endpos[2]+1])
    else
        H, E, F, best_endpos = affinegap_local_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
        a′, b′ = traceback(a, b, H, E, F, best_endpos, subst_matrix, gap_open_penalty, gap_extend_penalty)
        return AlignmentResult(H[best_endpos[1]+1,best_endpos[2]+1], a′, b′)
    end
    error("not implemented")
end

function pairalign(::EditDistance, a, b, cost::AbstractCostModel;
                   distance_only::Bool=false)
    subst_matrix = cost.subst_matrix
    insertion_cost = cost.insertion_cost
    deletion_cost = cost.deletion_cost
    if distance_only
        D = edit_distance(a, b, subst_matrix, insertion_cost, deletion_cost)
        return AlignmentResult(D[end,end])
    else
        D = edit_distance(a, b, subst_matrix, insertion_cost, deletion_cost)
        a′, b′ = traceback(a, b, D, subst_matrix, insertion_cost, deletion_cost)
        return AlignmentResult(D[end,end], a′, b′)
    end
end

function pairalign(::LevenshteinDistance, a, b;
                   distance_only::Bool=false)
    unit_cost = UnitCostModel{Int}()
    insertion_cost = 1
    deletion_cost  = 1
    cost = CostModel(unit_cost, insertion_cost, deletion_cost)
    return pairalign(EditDistance(), a, b, cost,
                     distance_only=distance_only)
end

function pairalign(::HammingDistance, a, b;
                   distance_only::Bool=false)
    if length(a) != length(b)
        error("two sequences should be same length")
    end
    if distance_only
        return AlignmentResult(hamming_distance(Int, a, b))
    else
        return AlignmentResult(hamming_distance(Int, a, b), a, b)
    end
end
