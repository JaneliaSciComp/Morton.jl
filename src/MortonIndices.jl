"""
    MortonIndices{N,A,I}

# Parameters
* N - Number of dimensions
* A - Axes Type
* I - Type of integer holding the Morton index
"""
struct MortonIndices{N,A,I} <: AbstractArray{MortonIndex{I},N}
    cartesian::CartesianIndices{N,A}
end
MortonIndices(cis::CartesianIndices) = MortonIndices{Int}(cis)
MortonIndices{I}(cis::CartesianIndices{N,A}) where {N,A,I} = MortonIndices{N,A,I}(cis)
MortonIndices(x::Tuple) = MortonIndices{Int}(x)
MortonIndices{I}(x::NTuple{N,T}) where {N,T <: Integer,I} = MortonIndices{I}(CartesianIndices(x))

function Base.getindex(mis::MortonIndices{N,A,II}, I::Int...) where {N,A,II}
    ci = getindex(mis.cartesian, I...)
    MortonIndex(ci)
end
function iterate(mis::MortonIndices)
    t = iterate(mis.cartesian)
    MortonIndex(first(t)), last(t)
end
function iterate(mis::MortonIndices, state)
    t = iterate(mis.cartesian, state)
    MortonIndex(first(t)), last(t)
end
function Base.in(i::MortonIndex, mis::MortonIndices{N}) where N
    in(CartesianIndex{N}(i), mis.cartesian)
end
Base.LinearIndices(mis::MortonIndices) = LinearIndices(mis.cartesian)
Base.:(==)(a::MortonIndices, b::MortonIndices) = a.cartesian == b.cartesian
Base.first(mis::MortonIndices) = MortonIndex(first(mis.cartesian))
Base.last(mis::MortonIndices) = MortonIndex(last(mis.cartesian))
Base.size(mis::MortonIndices) = size(mis.cartesian)
Base.length(mis::MortonIndices) = length(mis.cartesian)
Base.ndims(mis::MortonIndices) = ndims(mis.cartesian)
Base.reverse(mis::MortonIndices) = MortonIndices(reverse(mis.cartesian))