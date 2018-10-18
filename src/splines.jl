mutable struct Spline
    α::AbstractVector
    k::Integer
end

#=
\[f(x)=\sum_i\alpha_iB_{i,k,t}(x)\equiv B(x)\vec{\alpha} = \langle B|\vec{\alpha}\rangle.\]
=#

(S::Spline)(Bᵢ::AbstractMatrix) = Bᵢ*S.α
function (S::Spline)(B::Vector{M}) where {M<:AbstractMatrix}
    kB = length(B)
    S.k > kB && error("No support for splines of order $(S.k) in basis of order $(kB)")
    S(B[S.k])
end
(S::Spline)(B::Basis) =  S(B.B)

Spline(f::Function, basis::Basis) = Spline(basis.B[end]\f.(locs(basis)), order(basis.t))

export Spline
