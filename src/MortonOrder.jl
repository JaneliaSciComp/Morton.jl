"""
    MortonOrder{N,A}
    MortonOrder(x::Union{CartesianIndices,AbstractArray,Tuple})

An iterator in Morton order over cartesian indices.

# Parameters
* N - Number of dimensions
* A - Axes Type
"""
struct MortonOrder{N,A}
    cartesian::CartesianIndices{N,A}
end
MortonOrder(x::Tuple) = MortonOrder(CartesianIndices(x))
MortonOrder(x::AbstractArray) = MortonOrder(CartesianIndices(x))
function Base.iterate(mo::MortonOrder{N}) where N
    if isempty(mo.cartesian)
        return nothing
    else
        m = MortonIndex(1)
        return CartesianIndex{N}(m), m
    end
end
function Base.iterate(mo::MortonOrder{N}, mi::MortonIndex) where N
    last_ci = last(mo.cartesian)
    last_mi = MortonIndex(last_ci)

    mi  = MortonIndex(mi.m + 1)
    ci = CartesianIndex{N}(mi)
    # Optimize this by finding the large power of 2 that will fit in dimensions as initial bound
    while ci ∉ mo.cartesian && mi.m < last_mi.m
        mi  = MortonIndex(mi.m + 1)
        ci = CartesianIndex{N}(mi)
    end
    if mi.m > last_mi.m
        return nothing
    end
    return ci, mi
end
Base.length(mo::MortonOrder) = length(mo.cartesian)
Base.eltype(mo::MortonOrder) = eltype(mo.cartesian)
Base.IteratorSize(::Type{MortonOrder{N,A}}) where {N,A} = Base.HasLength()
Base.IteratorEltype(::Type{MortonOrder{N,A}}) where {N,A} = Base.IteratorEltype(CartesianIndices{N,A})
function collect_by_sorting(mo::MortonOrder)
    mis = MortonIndex.(mo.cartesian)[:]
    p = sortperm(mis)
    return mo.cartesian[p]
end

# Allow indexing 
function Base.getindex(A::AbstractArray, ::Type{MortonOrder})
    mo = MortonOrder(A)
    return A[collect(mo)]
end