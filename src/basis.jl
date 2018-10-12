using LinearAlgebra
using SparseArrays
using BandedMatrices

import Base: show

struct Basis
    t::AbstractKnotSet
    x::AbstractVector
    w::AbstractVector
    B::Vector{AbstractMatrix}
end

# See http://pages.cs.wisc.edu/~deboor/pgs/bsplvb.f
function evaluate!(Bᵢ, t::AbstractKnotSet, x::AbstractVector)
    #=
    \[B_{i,1,t}(x)=\begin{cases}
    1, & t_i\leq x < t_{i+1}\\
    0, & \textrm{else}
    \end{cases}\]
    =#
    for i = 1:length(t)-1
        Bᵢ[1][:,i] = (t[i] .<= x .< t[i+1])
    end

    #=
    \[B_{i,k,t}(x)=\frac{x-t_i}{t_{i+k-1}-t_i}B_{i,k-1,t}(x)-
    \frac{t_{i+k}-x}{t_{i+k}-t_{i+1}}B_{i+1,k-1,t}(x)\]
    =#
    for kk = 2:order(t)
        for i = 1:length(t) - kk
            d₁ = t[i+kk-1]-t[i]
            d₂ = t[i+kk]-t[i+1]
            f₁ = d₁ == 0 ? 0 : (x .- t[i])/d₁
            f₂ = d₂ == 0 ? 0 : (t[i+kk] .- x)/d₂
            Bᵢ[kk][:,i] = f₁.*Bᵢ[kk-1][:,i] + f₂.*Bᵢ[kk-1][:,i+1]
        end
    end
end

function Basis(t::AbstractKnotSet)
    x,w = lgwt(t)
    B = [spzeros(eltype(x), length(x), length(t)-kk)
         for kk = 1:order(t)]
    evaluate!(B, t, x)
    Basis(t, x, w, B)
end

function (basis::Basis)(x::AbstractVector)
    Bᵢ = [spzeros(eltype(x), length(x), length(basis.t)-kk)
          for kk = 1:order(basis.t)]
    evaluate!(Bᵢ, basis.t, x)
    Bᵢ
end

function (basis::Basis)(::UniformScaling)
    m = size(basis.B[end],2)
    k = order(basis.t)
    B = BandedMatrix{eltype(basis.t)}(undef, m,m, 0, k-1)
    for j = 1:m
        for i = max(1,j-k):j
            B[i,j] = dot(basis.B[end][:,i], (weights(basis) .* basis.B[end][:,j]))
        end
    end
    Symmetric(B)
end

locs(basis::Basis) = basis.x
weights(basis::Basis) = basis.w

function show(io::IO, basis::Basis)
    write(io, "BSpline basis with $(basis.t)")
end

@recipe function plot(basis::Basis, n=100)
    x = range(first(basis.t), stop=last(basis.t),
              length=n*numintervals(basis.t))
    x, Matrix(basis(x)[end])
end

export locs, weights
