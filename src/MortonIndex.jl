"""
    MortonIndex{I <: Integer}

Represents a 1-based Morton index.
* `I` is a type parameter storing the Morton index itself.
"""
struct MortonIndex{I <: Integer}
    m::I
end
MortonIndex(ci::CartesianIndex{1}) = MortonIndex(ci[1])
MortonIndex(ci::CartesianIndex{2}) = MortonIndex(cartesian2morton(ci))
MortonIndex(ci::CartesianIndex{3}) = MortonIndex(cartesian3morton(ci))
MortonIndex(ci::CartesianIndex{N}) where N = MortonIndex(cart)
Base.isless(a::MortonIndex, b::MortonIndex) = a.m < b.m
Base.getindex(mi::MortonIndex) = mi.m

Base.CartesianIndex{1}(mi::MortonIndex) = CartesianIndex{1}(mi.m)
Base.CartesianIndex{2}(mi::MortonIndex) = morton2cartesianindex(mi.m)
Base.CartesianIndex{3}(mi::MortonIndex) = morton3cartesianindex(mi.m)

function Base.convert(::Type{CartesianIndex{N}}, mi::MortonIndex)::CartesianIndex{N} where N
    CartesianIndex{N}(mi)
end

function Base.convert(::Type{MortonIndex}, ci::Base.CartesianIndex)::MortonIndex
    MortonIndex(ci)
end

function Base.getindex(A::AbstractArray{T,N}, mi::MortonIndex) where {T,N}
    ci = CartesianIndex{N}(mi)
    getindex(A, ci)
end
function Base.getindex(A::AbstractArray{T,N}, MI::AbstractVector{<: MortonIndex}) where {T, N}
    CI = CartesianIndex{N}.(MI)
    Base.getindex(A, CI)
end
function Base.setindex!(A::AbstractArray{T,N}, v, mi::MortonIndex) where {T,N}
    ci = CartesianIndex{N}(mi)
    setindex!(A, v, ci)
end


(::Type{I})(mi::MortonIndex{I}) where I = mi.m
(::Type{N})(mi::MortonIndex{I}) where {N <: Number, I} = N(mi.m)
Base.convert(::Type{I}, mi::MortonIndex{I}) where I = mi.m
Base.convert(::Type{N}, mi::MortonIndex{I}) where {N <: Number, I} = convert(N,mi.m)

