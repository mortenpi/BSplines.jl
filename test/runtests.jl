using BSplines
using Test

@testset "Knot sets" begin
    t = LinearKnotSet(7, 0, 1, 10)
    @test order(t) == 7
    @test numintervals(t) == 10
    @test length(t) == 23
    @test first(t) == 0
    @test last(t) == 1
    for i = 1:7
        @test t[i] == 0
        @test t[end-i+1] == 1
    end
    tt = range(0, stop=1, length=11)
    for i = 7:17
        @test t[i] == tt[i-6]
    end

    ttt = collect(t)
    @test length(t) == length(ttt)
    for i ∈ eachindex(t)
        @test t[i] == ttt[i]
    end

    @test eltype(t) <: Real
end

@testset "Quadrature" begin
    t = LinearKnotSet(1, 0, 1, 2)
    x,w = BSplines.lgwt(t)
    @test all(w .== 1/4)
    @test x == [-1,1,-1,1]/(4*√3) + [1,1,3,3]/4
end
