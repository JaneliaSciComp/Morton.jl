struct MortonIndex{I}
    m::I
end
Base.CartesianIndex{2}(mi::MortonIndex) = morton2cartesianindex(mi.m)
Base.CartesianIndex{3}(mi::MortonIndex) = morton3cartesianindex(mi.m)
MortonIndex(ci::CartesianIndex{2}) = MortonIndex(cartesian2morton(ci))
MortonIndex(ci::CartesianIndex{3}) = MortonIndex(cartesian3morton(ci))

function Base.convert(::Type{CartesianIndex{N}}, mi::MortonIndex)::CartesianIndex{N} where N
    CartesianIndex{N}(mi)
end

function Base.convert(::Type{MortonIndex}, ci::CartesianIndex)::MortonIndex
    MortonIndex(ci)
end

function Base.getindex(A::AbstractArray{T,N}, mi::MortonIndex) where {T,N}
    ci = CartesianIndex{N}(mi)
    getindex(A, ci)
end

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