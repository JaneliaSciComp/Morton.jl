# Morton

[![Build Status](https://travis-ci.org/JaneliaSciComp/Morton.jl.svg?branch=master)](https://travis-ci.org/JaneliaSciComp/Morton.jl)  [![codecov.io](http://codecov.io/github/JaneliaSciComp/Morton.jl/coverage.svg?branch=master)](http://codecov.io/github/JaneliaSciComp/Morton.jl?branch=master)

This package provides functions to convert between Morton number (a.k.a.
[Z-order](https://en.wikipedia.org/wiki/Z-order_curve)), [Cartesian
coordinates](https://en.wikipedia.org/wiki/Cartesian_coordinate_system),
and [quadtree](https://en.wikipedia.org/wiki/Quadtree) and
[octree](https://en.wikipedia.org/wiki/Octree) coordinates.

Say for example you have a 4x4 matrix.  The sixteen cells could be addressed in
each of the following three ways.

Morton order (column major):

| | | | |
|---|---|---|---|
1|3|9|11
2|4|10|12
5|7|13|15
6|8|14|16

Cartesian indices (column major):

| | | | |
|---|---|---|---|
1,1|1,2|1,3|1,4
2,1|2,2|2,3|2,4
3,1|3,2|3,3|3,4
4,1|4,2|4,3|4,4

Quadtree coordinates (column major):

| | | | |
|---|---|---|---|
1,1|3,1|3,1|3,3
2,1|4,1|3,2|3,4
2,1|2,3|4,1|4,3
2,2|2,4|4,2|4,4

To convert from Morton to Cartesian, use the `morton2cartesian` function:

```
julia> Pkg.add("Morton")
julia> using Morton

julia> morton2cartesian(13)
2-element Array{Int64,1}:
3
3

julia> morton2cartesianindex(21)
CartesianIndex(7, 1)

julia> CartesianIndex{2}(MortonIndex(29))
CartesianIndex(7, 3)
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

julia> cartesian2morton(CartesianIndex(3,3))
13

julia> MortonIndex(CartesianIndex(3,3))
MortonIndex{Int64}(13)

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

There are also un-exported N-dimensional functions to convert between tree and
Morton, and tree and Cartesian (e.g. `Morton._treeNmorton`).  Please let me
know if you have a clever way to convert directly between Morton and Cartesian
in arbitrary dimensions!


# Related packages

[RegionTrees](https://github.com/rdeits/RegionTrees.jl)


# Author

[Ben Arthur](http://www.janelia.org/people/research-resources-staff/ben-arthur), arthurb@hhmi.org  
[Scientific Computing](http://www.janelia.org/research-resources/computing-resources)  
[Janelia Research Campus](http://www.janelia.org)  
[Howard Hughes Medical Institute](http://www.hhmi.org)
