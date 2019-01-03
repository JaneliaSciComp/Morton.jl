using Morton
using Test
using Random

@test cartesian2morton([5,2]) == 19
@test cartesian3morton([5,2,1]) == 67
@test morton2cartesian(19) == [5,2]
@test morton3cartesian(67) == [5,2,1]
@test tree2morton([2,1,3]) == 19
@test tree3morton([2,1,3]) == 67
@test morton2tree(19) == [2,1,3]
@test morton3tree(67) == [2,1,3]
@test tree2cartesian([2,1,3]) == [5,2]
@test tree3cartesian([2,1,3]) == [5,2,1]
@test cartesian2tree([5,2]) == [2,1,3]
@test cartesian3tree([5,2,1]) ==  [2,1,3]


Random.seed!(1234)

m = rand(0x01:0xff,10)
@test map(cartesian2morton ∘ morton2cartesian, m) == m
@test map(cartesian3morton ∘ morton3cartesian, m) == m
@test map(tree2morton ∘ morton2tree, m) == m
@test map(tree3morton ∘ morton3tree, m) == m

for i=1:10
    c = rand(0x01:0xff,2)
    @test (tree2cartesian ∘ cartesian2tree)(c) == c
    @test (morton2cartesian ∘ cartesian2morton)(c) == c
    c = rand(0x01:0xff,3)
    @test (tree3cartesian ∘ cartesian3tree)(c) == c
    @test (morton3cartesian ∘ cartesian3morton)(c) == c
end

for i=1:10
    t = rand(1:4,10)
    while t[1]==1; popfirst!(t); end
    @test (morton2tree ∘ tree2morton)(t) == t
    @test (cartesian2tree ∘ tree2cartesian)(t) == t
    t = rand(1:8,10)
    while t[1]==1; popfirst!(t); end
    @test (morton3tree ∘ tree3morton)(t) == t
    @test (cartesian3tree ∘ tree3cartesian)(t) == t
end
