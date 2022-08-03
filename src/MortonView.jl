struct MortonView{T, A <: AbstractArray{T}} <: AbstractVector{T}
    parent::A
end
MortonView(vector::A) where {T, A <: AbstractVector{T}} = MortonView{T, A}(vector)
function MortonView(parent::A) where {T, N, A <: AbstractArray{T,N}}
    all(ispow2.(size(parent))) ||
        throw(DimensionMismatch("The parent of a MortonView must have size dimensions that are a power of 2."))
    firstsize = size(parent, 1)
    not_halved = true
    for d in 2:N
        if size(parent, d) == firstsize
            continue
        elseif size(parent, d) == firstsize÷2 && not_halved
            firstsize = firstsize÷2
            not_halved = false
        else
            throw(DimensionMismatch("The parent of a MortonView must be a hypercube with equal size dimensions."))
        end
    end
    MortonView{T, A}(parent)
end
Base.size(mv::MortonView) = (length(mv.parent),)
function Base.getindex(mv::MortonView, i::Int)
    return getindex(mv.parent, MortonIndex(i))
end
function Base.setindex!(mv::MortonView, v, i::Int)
    setindex!(mv.parent, v, MortonIndex(i))
end
Base.IndexStyle(::Type{MortonView}) = IndexLinear()