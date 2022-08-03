module Morton

export cartesian2morton, cartesian3morton
export morton2cartesian, morton3cartesian
export morton2cartesianindex, morton3cartesianindex
export tree2morton, tree3morton
export morton2tree, morton3tree
export tree2cartesian, tree3cartesian
export cartesian2tree, cartesian3tree
export MortonIndex, MortonIndices, MortonOrder, MortonView

### https://fgiesen.wordpress.com/2009/12/13/decoding-morton-codes/

# "Insert" a 0 bit after each of the 16 low bits of x
function _Part1By1(x::Integer)
    x &= 0x0000ffff                   # x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x | (x <<  8)) & 0x00ff00ff  # x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x | (x <<  4)) & 0x0f0f0f0f  # x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x | (x <<  2)) & 0x33333333  # x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x | (x <<  1)) & 0x55555555  # x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x+1
end

# "Insert" two 0 bits after each of the 10 low bits of x
function _Part1By2(x::Integer)
    x &= 0x000003ff                   # x = ---- ---- ---- ---- ---- --98 7654 3210
    x = (x | (x << 16)) & 0xff0000ff  # x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x <<  8)) & 0x0300f00f  # x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x <<  4)) & 0x030c30c3  # x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x <<  2)) & 0x09249249  # x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x+1
end

@inline function _Part1ByN(x::Integer, ::Val{N}) where N
    @assert(N > 0, "N must be a positive integer.")
    d = digits(Bool, x, base=2)
    # numdigits = (length(d)-1)*(N+1) + 1
    # powers = 0:(N+1):numdigits
    sum(d .<< ((0:length(d)-1)*(N+1))) + 1
end
#@inline mask(bits) = mapreduce(x->1 << (x-1), +, bits)
@inline function mask(bits, I=Int32)
    m = zero(I)
    for b in bits
        m += one(I) << (b-1)
    end
    return m
end
@inline function _Part1ByN_fast(x::I, ::Val{N}) where {I <: Integer,N}
    #nbits = sizeof(I)*8
    #ndigits = nbits ÷ N
    x &= mask(1:10)
    x = (x | (x << 16)) & (mask(1:8) | mask(25:32))
    x = (x | (x << 8))  & (mask(1:12:32) | mask(2:12:32) | mask(3:12:32) | mask(4:12:32))
    x = (x | (x << 4))  & (mask(1:6:32) | mask(2:6:32))
    x = (x | (x << 2))  & mask(1:3:32)
    x + 1
end
@inline function _Part1By2_fast(x::Integer)
    #nbits = sizeof(I)*8
    #ndigits = nbits ÷ N
    x &= mask(1:10)
    x = (x | (x << 16)) & (mask(1:8) | mask(25:32))
    x = (x | (x << 8))  & (mask(1:12:32) | mask(2:12:32) | mask(3:12:32) | mask(4:12:32))
    x = (x | (x << 4))  & (mask(1:6:32) | mask(2:6:32))
    x = (x | (x << 2))  & mask(1:3:32)
    x + 1
end
function _generate_Part1ByN(I, N)
    nbits = sizeof(I)*8
    S = I(N + 1)
    ndigits = nbits ÷ S
    num_shifts = nbits - leading_zeros(I(nbits)) - 1
    v = Vector{Any}()
    push!(v,:(
        x &= $(mask(1:ndigits, unsigned(I)))
    ))
    for i in 2:num_shifts
        p = 2^(num_shifts-i)
        mp = I(0x0)
        for m in I(0x1):I(p)
            mp |= mask(m:(S*p):nbits, unsigned(I))
        end
        push!(v,:(
            x = (x | x << $(2^(N-1)*p)) & $mp
        ))
    end
    push!(v, :(($I(x)+one($I))))
    return Expr(:block, v...)
end
    #=
    quote
        x &= $(mask(1:ndigits))
        x = (x | (x << 16)) & $(mask(1:8) | mask(25:32))
        x = (x | (x << 8))  & $(mask(1:12:32) | mask(2:12:32) | mask(3:12:32) | mask(4:12:32))
        x = (x | (x << 4))  & $(mask(1:6:32) | mask(2:6:32))
        x = (x | (x << 2))  & $(mask(1:3:32))
        x + 1
    end
    =#

@generated function _Part1By1_gen(x::I)::I where I <: Integer
    _generate_Part1ByN(I, 1)
end

@generated function _Part1By2_gen(x::I)::I where I <: Integer
    _generate_Part1ByN(I, 2)
end

#_Part1ByN(x::Integer, ::Val{1}) = _Part1By1(x)
#_Part1ByN(x::Integer, ::Val{2}) = _Part1By2(x)

"""
    cartesian2morton(c::Vector) -> m::Integer

Given a 2-D Cartesian coordinate, return the corresponding Morton number.

# Examples
```jldoctest
julia> cartesian2morton([5,2])
19
```
"""
function cartesian2morton(c::Union{AbstractVector{T},CartesianIndex{2}}) where T<:Integer
    (_Part1By1(c[2]-1) << 1) + _Part1By1(c[1]-1) - 2
end

"""
    cartesian3morton(c::AbstractVector) -> m::Integer

Given a 3-D Cartesian coordinate, return the corresponding Morton number.

# Examples
```jldoctest
julia> cartesian3morton([5,2,1])
67
```
"""
function cartesian3morton(c::Union{AbstractVector{T},CartesianIndex{3}}) where T<:Integer
    (_Part1By2(c[3]-1) << 2) + (_Part1By2(c[2]-1) << 1) + _Part1By2(c[1]-1) - 6
end


function _Compact1By1(x::Integer)
    x &= 0x55555555                   # x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x = (x | (x >>  1)) & 0x33333333  # x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x | (x >>  2)) & 0x0f0f0f0f  # x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x | (x >>  4)) & 0x00ff00ff  # x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x | (x >>  8)) & 0x0000ffff  # x = ---- ---- ---- ---- fedc ba98 7654 3210
    x+1
end

function _Compact1By2(x::Integer)
    x &= 0x09249249                   # x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x = (x | (x >>  2)) & 0x030c30c3  # x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x >>  4)) & 0x0300f00f  # x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x >>  8)) & 0xff0000ff  # x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x >> 16)) & 0x000003ff  # x = ---- ---- ---- ---- ---- --98 7654 3210
    x+1
end

function _Compact1By2_64(x::Integer)
    x &= 0x1249249249249249           # x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
    x = (x | (x >>  2)) & 0x030c30c3  # x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x >>  4)) & 0x0300f00f  # x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x >>  8)) & 0xff0000ff  # x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x >> 16)) & 0x000003ff  # x = ---- ---- ---- ---- ---- --98 7654 3210
    x+1
end


function _Compact1By3(x::Integer)
    x &= 0x000000111111111             # x = ---9 ---8 ---7 ---6 ---5 ---4 ---3 ---2 ---1 ---0
    x = (x | (x >>  4)) & 0x030c30c3   # x = ---- --98 ---- --76 ---- --54 ---- --32 ---- --10
    x = (x | (x >>  8)) & 0x0300f00f   # x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x >>  16)) & 0xff0000ff  # x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x >>  32)) & 0x000003ff  # x = ---- ---- ---- ---- ---- --98 7654 3210
    x+1
end

function _Compact1ByN(x::Integer, ::Val{N}) where N
    @assert(N > 0, "N must be positive")
    d = digits(Bool, x, base=2)
    bitidx = 1:(N+1):length(d)
    sum(d[bitidx] .<< (0:length(bitidx)-1)) + 1
end
_Compact1ByN(x::Integer, ::Val{1}) = _Compact1By1(x)
_Compact1ByN(x::Integer, ::Val{2}) = _Compact1By2(x)

function _generate_Compact1ByN(I, N)
    nbits = sizeof(I)*8
    S = I(N + 1)
    ndigits = nbits ÷ S
    num_shifts = nbits - leading_zeros(I(nbits)) - 1
    v = Vector{Any}()
    push!(v,:(
        x &= $(mask(1:S:(ndigits*S), unsigned(I)))
    ))
    for i in 2:num_shifts
        p = 2^(i-1)
        mp = I(0x0)
        for m in I(0x1):I(p)
            mp |= mask(m:(S*p):nbits, unsigned(I))
        end
        push!(v,:(
            x = (x | x >> $(p ÷ 2^(2-N))) & $mp
        ))
    end
    push!(v, :(($I(x)+one($I))))
    return Expr(:block, v...)
end

@generated function _Compact1ByN_gen(x::I, ::Val{N}) where {I <: Integer, N}
    _generate_Compact1ByN(I, N)
end

"""
    morton2cartesian(m::Integer) -> [x,y]

Given a Morton number, return the corresponding 2-D Cartesian coordinates.

# Examples
```jldoctest
julia> morton2cartesian(19)
2-element Array{Int64,1}:
 5
 2
```
"""
function morton2cartesian(m::Integer)
    m -= 1
    [_Compact1By1(m>>0), _Compact1By1(m>>1)]
end

function morton2cartesianindex(m::Integer)
    m -= 1
    CartesianIndex(_Compact1By1(m>>0), _Compact1By1(m>>1))
end



"""
    morton3cartesian(m::Integer) -> [x,y,z]

Given a Morton number, return the corresponding 3-D Cartesian coordinates.

# Examples
```jldoctest
julia> morton3cartesian(67)
3-element Array{Int64,1}:
 5
 2
 1
```
"""
function morton3cartesian(m::Integer)
    m -= 1
    [_Compact1By2(m>>0), _Compact1By2(m>>1), _Compact1By2(m>>2)]
end

function morton3cartesianindex(m::Integer)
    m -= 1
    CartesianIndex(_Compact1By2(m>>0), _Compact1By2(m>>1), _Compact1By2(m>>2))
end

function _treeNmorton(t::AbstractVector{T}, ndim::Integer) where T<:Integer
    n=m=0
    it=length(t)
    ndim2=2^ndim
    while it>0
      m += (t[it]-1)*ndim2^n
      n += 1
      it-=1
    end
    m+1
end

"""
    tree2morton(t::AbstractVector) -> m::Integer

Given a quadtree coordinate, return the corresponding Morton number.

# Examples
```jldoctest
julia> tree2morton([2,1,3])
19
```
"""
tree2morton(t::AbstractVector{T}) where T<:Integer = _treeNmorton(t,2)

"""
    tree3morton(t::AbstractVector) -> m::Integer

Given a octree coordinate, return the corresponding Morton number.

# Examples
```jldoctest
julia> tree3morton([2,1,3])
67
```
"""
tree3morton(t::AbstractVector{T}) where T<:Integer = _treeNmorton(t,3)


function _mortonNtree(m::T, ndim::Integer) where T<:Integer
    t=T[]
    ndim2=2^ndim
    while true
        d,r = [divrem(m-1,ndim2)...]+[1,1]
        pushfirst!(t,r)
        d==1 && break
        m=d
    end
    t
end

"""
    morton2tree(m::Integer) -> t::AbstractVector

Given a Morton number, return the corresponding quadtree coordinate.

# Examples
```jldoctest
julia> morton2tree(19)
3-element Array{Any,1}:
 2
 1
 3
```
"""
morton2tree(m::Integer) = _mortonNtree(m,2)

"""
    morton3tree(m::Integer) -> t::AbstractVector

Given a Morton number, return the corresponding octree coordinate.

# Examples
```jldoctest
julia> morton3tree(67)
3-element Array{Any,1}:
 2
 1
 3
```
"""
morton3tree(m::Integer) = _mortonNtree(m,3)


function _treeNcartesian(t::AbstractVector{T}, ndim::Integer) where T<:Integer
    c = [((t[1]-1)>>b)&1+1 for b in 0:ndim-1]
    if length(t)>1
        cn = _treeNcartesian(t[2:end], ndim)
        m = 2^(length(t)-1)
        return [m*(c[i]-1)+cn[i] for i in 1:ndim]
    else
        return c
    end
end

"""
   tree2cartesian(t::AbstractVector) -> c::AbstractVector

Given quadtree coordinate, return the corresponding 2-D Cartesian coordinate.

# Examples
```jldoctest
julia> tree2cartesian([2,1,3])
2-element Array{Int64,1}:
 5
 2
```
"""
tree2cartesian(t::AbstractVector{T}) where T<:Integer = _treeNcartesian(t, 2)

"""
   tree3cartesian(t::AbstractVector) -> c::AbstractVector

Given octree coordinate, return the corresponding 3-D Cartesian coordinate.

# Examples
```jldoctest
julia> tree3cartesian([2,1,3])
3-element Array{Int64,1}:
 5
 2
 1
```
"""
tree3cartesian(t::AbstractVector{T}) where T<:Integer = _treeNcartesian(t, 3)


function _cartesianNtree(c::AbstractVector{T}, half, ndim::Integer) where T<:Integer
    t = 1
    for d=1:ndim
        t += 2^(d-1)*(c[d]>half)
    end
    if half>1
        return [t, _cartesianNtree(map(x->(x-1)%half+1,c), half>>1, ndim)...]
    else
        return [t]
    end
end

"""
   cartesian2tree(c::AbstractVector) -> t::AbstractVector

Given a 2-D Cartesian coordinate, return the corresponding quadtree coordinate.

# Examples
```jldoctest
julia> cartesian2tree([5,2])
3-element Array{Int64,1}:
 2
 1
 3
```
"""
cartesian2tree(c::AbstractVector{T}) where T<:Integer =
      _cartesianNtree(c, max(2,nextpow(2, widen(maximum(c))))>>1, 2)

"""
   cartesian3tree(c::AbstractVector) -> t::AbstractVector

Given a 3-D Cartesian coordinate, return the corresponding octree coordinate.

# Examples
```jldoctest
julia> cartesian3tree([5,2,1])
3-element Array{Int64,1}:
 2
 1
 3
```
"""
cartesian3tree(c::AbstractVector{T}) where T<:Integer =
      _cartesianNtree(c, max(2,nextpow(2, widen(maximum(c))))>>1, 3)

include("MortonIndex.jl")
include("MortonIndices.jl")
include("MortonOrder.jl")
include("MortonView.jl")

end # module
