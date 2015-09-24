function pairalign(::GlobalAlignment, a, b, score::AffineGap;
                   score_only::Bool=false, k_band::Int=-1, linear_space::Bool=score_only)
    if k_band < 0
        if score_only
            H, _, _ = affinegap_global_align(a, b, score)
            return AlignmentResult(H[end,end])
        else
            H, E, F = affinegap_global_align(a, b, score)
            a′, b′ = traceback(a, b, H, E, F, score)
            return AlignmentResult(H[end,end], a′, b′)
        end
    else
    end
    error("not implemented")
end

function pairalign(::LocalAlignment, a, b, score::AffineGap;
                   score_only::Bool=false, linear_space::Bool=score_only)
    if score_only
        H, _, _, best_endpos = affinegap_local_align(a, b, score)
        return AlignmentResult(H[best_endpos[1]+1,best_endpos[2]+1])
    else
        H, E, F, best_endpos = affinegap_local_align(a, b, score)
        a′, b′ = traceback(a, b, H, E, F, best_endpos, score)
        return AlignmentResult(H[best_endpos[1]+1,best_endpos[2]+1], a′, b′)
    end
    error("not implemented")
end
