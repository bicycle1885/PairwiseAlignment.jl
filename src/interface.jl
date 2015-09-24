# Interfaces
# ----------

function pairalign(::GlobalAlignment, a, b, score::AffineGap;
                   score_only::Bool=false, k_band::Int=-1, linear_space::Bool=score_only)
    subst_matrix = score.subst_matrix
    gap_open_penalty = score.gap_open_penalty
    gap_extend_penalty = score.gap_extend_penalty
    if k_band < 0
        if score_only
            H, _, _ = affinegap_global_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(H[end,end])
        else
            H, E, F = affinegap_global_align(a, b, subst_matrix, gap_open_penalty, gap_extend_penalty)
            a′, b′ = traceback(a, b, H, E, F, subst_matrix, gap_open_penalty, gap_extend_penalty)
            return AlignmentResult(H[end,end], a′, b′)
        end
    else
    end
    error("not implemented")
end

function pairalign(::LocalAlignment, a, b, score::AffineGap;
                   score_only::Bool=false, linear_space::Bool=score_only)
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
