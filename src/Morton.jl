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
function cartesian2morton(c::Vector{T}) where T<:Integer
  (_Part1By1(c[2]-1) << 1) + _Part1By1(c[1]-1) - 2
end

"""
    cartesian3morton(c::Vector) -> m::Integer

Given a 3-D Cartesian coordinate, return the corresponding Morton number.

# Examples
```jldoctest
julia> cartesian3morton([5,2,1])
67
```
"""
function cartesian3morton(c::Vector{T}) where T<:Integer
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


function _treeNmorton(t::Vector{T}, ndim::Integer) where T<:Integer
  n=m=0
  it=length(t)
  while it>0
    m += (t[it]-1)*(2^ndim)^n
    n += 1
    it-=1
  end
  m+1
end

"""
    tree2morton(t::Vector) -> m::Integer

Given a quadtree coordinate, return the corresponding Morton number.

# Examples
```jldoctest
julia> tree2morton([2,1,3])
19
```
"""
tree2morton(t::Vector{T}) where T<:Integer = _treeNmorton(t,2)

"""
    tree3morton(t::Vector) -> m::Integer

Given a octree coordinate, return the corresponding Morton number.

# Examples
```jldoctest
julia> tree3morton([2,1,3])
67
```
"""
tree3morton(t::Vector{T}) where T<:Integer = _treeNmorton(t,3)


function _mortonNtree(m::T, ndim::Integer) where T<:Integer
  t=T[]
  while true
    d,r = [divrem(m-1,ndim)...]+[1,1]
    unshift!(t,r)
    d==1 && break
    m=d
  end
  t
end

"""
    morton2tree(m::Integer) -> t::Vector

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
morton2tree(m::Integer) = _mortonNtree(m,4)

"""
    morton3tree(m::Integer) -> t::Vector

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
morton3tree(m::Integer) = _mortonNtree(m,8)


"""
   tree2cartesian(t::Vector) -> [x,y]

Given quadtree coordinate, return the corresponding 2-D Cartesian coordinate.

# Examples
```jldoctest
julia> tree2cartesian([2,1,3])
2-element Array{Int64,1}:
 5
 2
"""
function tree2cartesian(t::Vector{T}) where T<:Integer
  x = isodd(t[1]) ? 1 : 2
  y = t[1]<3 ? 1 : 2
  if length(t)>1
    xn,yn = tree2cartesian(t[2:end])
    m = 2^(length(t)-1)
    return [m*(x-1)+xn, m*(y-1)+yn]
  else
    return [x,y]
  end
end

"""
   tree3cartesian(t::Vector) -> [x,y,z]

Given octree coordinate, return the corresponding 3-D Cartesian coordinate.

# Examples
```jldoctest
julia> tree3cartesian([2,1,3])
3-element Array{Int64,1}:
 5
 2
 1
"""
function tree3cartesian(t::Vector{T}) where T<:Integer
  x = isodd(t[1]) ? 1 : 2
  y = (t[1]-1)&3<2 ? 1 : 2
  z = t[1]<5 ? 1 : 2
  if length(t)>1
    xn,yn,zn = tree3cartesian(t[2:end])
    m = 2^(length(t)-1)
    return [m*(x-1)+xn, m*(y-1)+yn, m*(z-1)+zn]
  else
    return [x,y,z]
  end
end


"""
   cartesian2tree(c::Vector) -> t::Vector

Given a 2-D Cartesian coordinate, return the corresponding quadtree coordinate.

# Examples
```jldoctest
julia> cartesian2tree([5,2])
3-element Array{Int64,1}:
 2
 1
 3
"""
cartesian2tree(c::Vector{T}) where T<:Integer =
      cartesian2tree(c, max(2,nextpow2(widen(maximum(c))))>>1)

function cartesian2tree(c::Vector{T}, half) where T<:Integer
  t = (c[1]>half) + 2*(c[2]>half) + 1
  if half>1
    return [t, cartesian2tree(map(x->(x-1)%half+1,c), half>>1)...]
  else
    return [t]
  end
end

"""
   cartesian3tree(c::Vector) -> t::Vector

Given a 3-D Cartesian coordinate, return the corresponding octree coordinate.

# Examples
```jldoctest
julia> cartesian3tree([5,2,1])
3-element Array{Int64,1}:
 2
 1
 3
"""
cartesian3tree(c::Vector{T}) where T<:Integer =
      cartesian3tree(c, max(2,nextpow2(widen(maximum(c))))>>1)

function cartesian3tree(c::Vector{T}, half) where T<:Integer
  t = (c[1]>half) + 2*(c[2]>half) + 4*(c[3]>half) + 1
  if half>1
    return [t, cartesian3tree(map(x->(x-1)%half+1,c), half>>1)...]
  else
    return [t]
  end
end

end # module
