module BSplines

using Polynomials: Polynomials, Poly, polyder

using RecipesBase

include("knot_sets.jl")
include("quadrature.jl")
include("basis.jl")
include("splines.jl")

struct PiecewisePoly2{T,U <: Real}
    ps :: Vector{Poly{T}}
    ts :: Vector{U}
    function PiecewisePoly2(ts::AbstractVector{U}, ps::AbstractVector{Poly{T}}) where {T, U}
        issorted(ts) || error("ts not sorted")
        zs = findall(x -> x != 0.0, diff(ts))
        @assert issorted(zs)
        _ps = ps[zs]
        _ts = ts[zs]
        push!(_ts, ts[end])
        @assert length(_ps) + 1 == length(_ts)
        new{T,U}(_ps, _ts)
    end
end
function (pp::PiecewisePoly2)(x::Real)
    (x < pp.ts[1] || x > pp.ts[end]) && return 0.0
    i = findfirst(i -> x <= pp.ts[i+1], 1:length(pp.ps))
    (pp.ps[i])(x)
end

function Polynomials.polyder(pp::PiecewisePoly2, k::Integer=1)
    PiecewisePoly2(pp.ts, polyder.(pp.ps, k))
end

"""
    abstract type AbstractBSpline

....

An `AbstractBSpline` should have the following interface:

* `order(s)`: order ``k`` of the B-spline
* `knot(s, i)`: value of the `i`-th knot

Additional supported methods that have default implementations:

* `length(s) = order(s) + 1`
* `getindex(s, i::Integer)`: return the value of the `i`-th knot
* `lastindex()`
"""
abstract type AbstractBSpline end

Base.length(s::AbstractBSpline) = order(s) + 1
Base.lastindex(s::AbstractBSpline) = order(s) + 1
Base.getindex(s::AbstractBSpline, i::Integer) = knot(s, i)
Base.getindex(s::AbstractBSpline, is::AbstractVector) = [knot(s, i) for i = is]
Base.lastindex(s::AbstractBSpline) = length(s)

function knot end

struct BSpline2 <: AbstractBSpline
    order :: Int
    # length(knots) == order + 1, sorted
    knots :: Vector{Float64}

    function BSpline2(order::Integer, knots::Vector{<:Real})
        order >= 1 || throw(DomainError(order, "B-spline order must be ≥ 1"))
        length(knots) == order + 1 || throw(DomainError("Bad number of knots ($(length(knots))) for B-spline of order $(order)"))
        new(order, sort(knots))
    end
end
BSpline2(knots::Vector{<:Real}) = BSpline2(length(knots)-1, knots)
order(s::BSpline2) = s.order
knot(s::BSpline2, i::Integer) = s.knots[i]

function _polyrep(s::BSpline2)
    k = order(s)
    if k == 1
        return [Poly(1.0)]
    else
        p = fill(Poly(0.0), k)
        s1, s2 = BSpline2(k - 1, s[1:end-1]), BSpline2(k - 1, s[2:end])
        Δt₁, Δt₂ = s1[end] - s1[1], s2[end] - s2[1]
        α₁, α₂ = Poly([-s1[1]/Δt₁, 1/Δt₁]), Poly([s2[end]/Δt₂, -1/Δt₂])
        p[1:end-1] .+= [α₁ * p for p in _polyrep(s1)]
        p[2:end] .+= [α₂ * p for p in _polyrep(s2)]
        return p
    end
end

polyrep(s::BSpline2) = PiecewisePoly2(s.knots, _polyrep(s))

end
