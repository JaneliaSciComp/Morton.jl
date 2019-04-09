module Morton

export cartesian2morton, cartesian3morton
export morton2cartesian, morton3cartesian
export tree2morton, tree3morton
export morton2tree, morton3tree
export tree2cartesian, tree3cartesian
export cartesian2tree, cartesian3tree

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

"""
    cartesian2morton(c::Vector) -> m::Integer

Given a 2-D Cartesian coordinate, return the corresponding Morton number.

# Examples
```jldoctest
julia> cartesian2morton([5,2])
19
```
"""
function cartesian2morton(c::AbstractVector{T}) where T<:Integer
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
function cartesian3morton(c::AbstractVector{T}) where T<:Integer
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

end # module
