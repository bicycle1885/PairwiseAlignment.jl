# abstract type of substitution matrix
abstract AbstractSubstitutionMatrix{T}

immutable UnitCostModel{T} <: AbstractSubstitutionMatrix{T} end
Base.getindex{T}(m::UnitCostModel{T}, x, y) = ifelse(x == y, T(0), T(1))

# abstract type of cost model
abstract AbstractCostModel{T}

immutable CostModel{T} <: AbstractCostModel{T}
    subst_matrix::AbstractSubstitutionMatrix{T}
    insertion_cost::T
    deletion_cost::T
end

function CostModel(subst_matrix; insertion_cost=nothing, deletion_cost=nothing)
    if insertion_cost === nothing || deletion_cost === nothing
        error("both insertion_cost and deletion_cost should be set")
    end
    return CostModel(subst_matrix, insertion_cost, deletion_cost)
end
