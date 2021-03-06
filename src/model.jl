# Substitution Matrices
# ---------------------

"""
Supertype of substitution matrix.

The required method:

    * `Base.getindex(submat, x, y) = <substitution score/cost from x to y>`
"""
abstract AbstractSubstitutionMatrix{T<:Real}

"""
Substitution matrix.
"""
type SubstitutionMatrix{T} <: AbstractSubstitutionMatrix{T}
    submat::Matrix{T}
end

Base.getindex(m::SubstitutionMatrix, x, y) = m.submat[convert(UInt8, x)+1,convert(UInt8, y)+1]
Base.setindex!(m::SubstitutionMatrix, v, x, y) = m.submat[convert(UInt8, x)+1,convert(UInt8, y)+1] = v


# Score Models
# ------------

"""
Supertype of score model.

Every score model is a problem of finding the maximum score and its alignment.
"""
abstract AbstractScoreModel{T<:Real}

"""
Affine gap scoring model.

The gap penalty of length `k` is `gap_open_penalty + gap_extend_penalty * k`.

Fields:

    * `submat`: a substitution matrix
    * `gap_open_penalty`: a penalty of opening a new gap
    * `gap_extend_penalty`: a penalty of extending a gap
"""
type AffineGapScoreModel{T} <: AbstractScoreModel{T}
    submat::AbstractSubstitutionMatrix{T}
    gap_open_penalty::T
    gap_extend_penalty::T
end

function AffineGapScoreModel{T}(submat::AbstractSubstitutionMatrix{T},
                                gap_open_penalty, gap_extend_penalty)
    return AffineGapScoreModel{T}(submat, T(gap_open_penalty), T(gap_extend_penalty))
end

function AffineGapScoreModel{T}(submat::AbstractSubstitutionMatrix{T};
                                gap_open_penalty=nothing, gap_extend_penalty=nothing)
    if gap_open_penalty === nothing || gap_extend_penalty === nothing
        error("both gap_open_penalty and gap_extend_penalty should be set")
    end
    return AffineGapScoreModel(submat, gap_open_penalty, gap_extend_penalty)
end

function AffineGapScoreModel{T}(submat::AbstractMatrix{T};
                                gap_open_penalty=nothing, gap_extend_penalty=nothing)
    if gap_open_penalty === nothing || gap_extend_penalty === nothing
        error("both gap_open_penalty and gap_extend_penalty should be set")
    end
    return AffineGapScoreModel(SubstitutionMatrix(submat), gap_open_penalty, gap_extend_penalty)
end


# Cost Models
# -----------

immutable UnitSubstitutionCost{T} <: AbstractSubstitutionMatrix{T} end
Base.getindex{T}(::UnitSubstitutionCost{T}, x, y) = ifelse(x == y, T(0), T(1))

"""
Supertype of cost model.

Every cost model is a problem of finding the minimum cost and its alignment.
"""
abstract AbstractCostModel{T}

"""
Cost model.

Fields:

    * `submat`: a substitution matrix
    * `insertion_cost`: a cost of inserting a character into the first sequence
    * `deletion_cost`: a cost of deleting a character from the first sequence
"""
type CostModel{T} <: AbstractCostModel{T}
    submat::AbstractSubstitutionMatrix{T}
    insertion_cost::T
    deletion_cost::T
end

function CostModel{T}(submat::AbstractSubstitutionMatrix{T},
                      insertion_cost, deletion_cost)
    return CostModel{T}(submat, insertion_cost, deletion_cost)
end

function CostModel{T}(submat::AbstractSubstitutionMatrix{T};
                      insertion_cost=nothing, deletion_cost=nothing)
    if insertion_cost === nothing || deletion_cost === nothing
        error("both insertion_cost and deletion_cost should be set")
    end
    return CostModel(submat, insertion_cost, deletion_cost)
end

function CostModel{T}(submat::AbstractMatrix{T};
                      insertion_cost=nothing, deletion_cost=nothing)
    if insertion_cost === nothing || deletion_cost === nothing
        error("both insertion_cost and deletion_cost should be set")
    end
    return CostModel(SubstitutionMatrix(submat), insertion_cost, deletion_cost)
end
