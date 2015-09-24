# the gap penalty of length `k` is gap_open_penalty + gap_extend_penalty * k
immutable AffineGap{T<:Real,M<:AbstractMatrix}
    subst_matrix::M  # I hope I could write as `M{T}`
    gap_open_penalty::T
    gap_extend_penalty::T
end

function AffineGap{T}(subst_matrix::AbstractMatrix{T};
                      gap_open_penalty=nothing, gap_extend_penalty=nothing)
    if gap_open_penalty === nothing || gap_extend_penalty === nothing
        error("both gap_open_penalty and gap_extend_penalty should be set")
    end
    return AffineGap(subst_matrix, gap_open_penalty, gap_extend_penalty)
end


function pairalign(::GlobalAlignment, a, b, score::AffineGap;
                   score_only::Bool=false, k_band::Int=-1, linear_space::Bool=score_only)
    if k_band < 0
        if score_only
            H, _, _ = global_affinegap_align(a, b, score)
            return AlignmentResult(H[end,end])
        else
            H, E, F = global_affinegap_align(a, b, score)
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
        H, _, _, best_endpos = local_affinegap_align(a, b, score)
        return AlignmentResult(H[best_endpos[1]+1,best_endpos[2]+1])
    else
        H, E, F, best_endpos = local_affinegap_align(a, b, score)
        a′, b′ = traceback(a, b, H, E, F, best_endpos, score)
        return AlignmentResult(H[best_endpos[1]+1,best_endpos[2]+1], a′, b′)
    end
    error("not implemented")
end
