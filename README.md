# Morton

[![Build Status](https://travis-ci.org/bjarthur/Morton.jl.svg?branch=master)](https://travis-ci.org/bjarthur/Morton.jl) [![Coverage Status](https://coveralls.io/repos/bjarthur/Morton.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/bjarthur/Morton.jl?branch=master) [![codecov.io](http://codecov.io/github/bjarthur/Morton.jl/coverage.svg?branch=master)](http://codecov.io/github/bjarthur/Morton.jl?branch=master)

This package provides functions to convert between Morton number (a.k.a.
[Z-order](https://en.wikipedia.org/wiki/Z-order_curve)), [Cartesian
coordinates](https://en.wikipedia.org/wiki/Cartesian_coordinate_system),
and [quadtree](https://en.wikipedia.org/wiki/Quadtree) and
[octree](https://en.wikipedia.org/wiki/Octree) coordinates.

Say for example you have a 4x4 matrix.  The sixteen cells could be addressed in
each of the following three ways.

Morton order:

| | | | |
|---|---|---|---|
1|2|5|6
3|4|7|8
9|10|13|14
11|12|15|16

Cartesian coordinates:

| | | | |
|---|---|---|---|
1,1|2,1|3,1|4,1
1,2|2,2|3,2|4,2
1,3|2,3|3,3|4,3
1,4|2,4|3,4|4,4

Quadtree coordinates:

| | | | |
|---|---|---|---|
1,1|1,2|2,1|2,2
1,3|1,4|2,3|2,4
3,1|3,2|4,1|4,2
3,3|3,4|4,3|4,4

To convert from Morton to Cartesian, use the `morton2cartesian` function:

```
julia> Pkg.add("Morton")
julia> using Morton

julia> morton2cartesian(13)
2-element Array{Int64,1}:
3
3
```

Similarly, one can convert from Morton to quadtree, or Cartesian to quadtree:

```
julia> morton2tree(13)
2-element Array{Int64,1}:
4
1

julia> cartesian2tree([3,3])
2-element Array{Int64,1}:
4
1
```

Of course each of the functions can be reversed:

```
julia> cartesian2morton([3,3])
13

julia> tree2morton([4,1])
13

julia> tree2cartesian([4,1])
2-element Array{Int64,1}:
3
3
```

Corresponding functions also exist for three dimensional matrices (i.e.
octrees).  Simply replace the 2 with a 3:  `morton3cartesian`, `morton3tree`,
etc.


# Related packages

[RegionTrees](https://github.com/rdeits/RegionTrees.jl)


# Author

[Ben Arthur](http://www.janelia.org/people/research-resources-staff/ben-arthur), arthurb@hhmi.org  
[Scientific Computing](http://www.janelia.org/research-resources/computing-resources)  
[Janelia Research Campus](http://www.janelia.org)  
[Howard Hughes Medical Institute](http://www.hhmi.org)
