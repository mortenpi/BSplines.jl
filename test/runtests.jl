using BSplines
using Test
using LinearAlgebra
using SparseArrays

function vecdist(a::AbstractVector, b::AbstractVector,
                 ϵ = eps(eltype(a)))
    δ = √(sum(abs2, a-b))
    δ, δ/√(sum(abs2, a .+ ϵ))
end

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

@testset "Basis" begin
    t = LinearKnotSet(3, 0, 1, 2)
    basis = BSplines.Basis(t)
    @testset "Eval on subintervals" begin
        x₁ = range(-1,stop=-0.5,length=10)
        Bᵢ₁ = basis(x₁)
        @test norm(Bᵢ₁) == 0

        function testbasis(x)
            Bᵢ = basis(x)
            B̃ = spzeros(Float64, length(x), 4)
            B̃[:,1] = (x .>= 0) .* (x .< 0.5) .* ((2x .- 1).^2)
            B̃[:,2] = (x .>= 0) .* (x .< 0.5) .* (2/3*(1 .- (3x .- 1).^2)) +
                (x .>= 0.5) .* (x .< 1)  .* (2*(x .- 1).^2)
            B̃[:,3] = (x .>= 0) .* (x .< 0.5) .* (2*x.^2) +
                (x .>= 0.5) .* (x .< 1)  .* (2/3*(1 .- (3x .- 2).^2))
            B̃[:,4] = (x .>= 0.5) .* (x .< 1) .* ((2x .- 1).^2)
            for j = 1:4
                δ,δr = vecdist(Bᵢ[end][:,j],B̃[:,j])
                @test δ < 1e-15
                @test δr < 1e-15
            end
        end
        testbasis(range(0,stop=1,length=50))
        testbasis(range(-0.5,stop=0.6,length=40))
        testbasis(range(0.5,stop=1.6,length=40))
    end

    @testset "Overlap matrix" begin
        B = basis(I)
        # # The reference is generated thusly:
        # Br = zeros{Any}(undef,4,4)
        # using SymPy
        # @syms x
        # Bs = [piecewise(((2x-1)^2, (Gt(x,0)) ∧ Lt(2x,1)),(0,true)),
        #       piecewise((2//3*(1-(3x-1)^2), (Gt(x,0)) ∧ Lt(2x,1)), (2*(x-1)^2, (Gt(2x,1)) ∧ Lt(x,1)),(0,true)),
        #       piecewise((2*x^2, (Gt(x,0)) ∧ Lt(2x,1)), (2//3*(1-(3x-2)^2), (Gt(2x,1)) ∧ Lt(x,1)),(0,true)),
        #       piecewise(((2x-1)^2, (Gt(2x,1)) ∧ Lt(x,1)),(0,true))]
        # for i = 1:4
        #     for j = 1:4
        #         Br[i,j] = integrate(Bs[i]*Bs[j],x)
        #     end
        # end
        Br = [  1/10  7/120  1/120      0
                7/120    1/6   1/10  1/120
                1/120   1/10    1/6  7/120
                0  1/120  7/120   1/10]
        @test norm(B-Br) < 1e-15
    end
end
