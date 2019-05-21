"""
    abstract type AbstractKnotSet{T <: Real}

Abstract parametric supertype of a B-spline knot set. The type parameter `T` is the type
of the knot values.

An `AbstractKnotSet` should have the following interface

* `order(t::AbstractKnotSet) -> Int`: returns the order ``k`` of the B-splines to be
  generated from this knot set
* `Base.length(t::AbstractKnotSet) -> Int`: returns the number of knots (``\\geq k + 1``)
* `Base.getindex(t::AbstractKnotSet, i::Integer) -> T`: returns the value of the `i`-th knot
"""
abstract type AbstractKnotSet{T <: Real} end

"""
    order(::AbstractKnotSet) -> Int

Returns the order ``k`` of the knot set.
"""
function order end

"""
    numintervals(::AbstractKnotSet) -> Int

Returns the number of ``(k + 1)``-element intervals in the knot set.
"""
numintervals(t::AbstractKnotSet) = length(t) - order(t)

knotsetname(t::AbstractKnotSet) = string(typeof(t))

Base.lastindex(t::AbstractKnotSet) = length(t)
Base.eachindex(t::AbstractKnotSet) = 1:length(t)
Base.first(t::AbstractKnotSet) = getindex(t, 1)
Base.last(t::AbstractKnotSet) = getindex(t, lastindex(t))


function Base.iterate(t::AbstractKnotSet, i=1)
    i > length(t) && return nothing
    t[i], i + 1
end

Base.eltype(t::AbstractKnotSet{T}) where T = T
#similar(t::AbstractKnotSet{T}) where T = Vector{T}(undef, length(t))

function Base.show(io::IO, t::AbstractKnotSet)
    write(io, "$(knotsetname(t)) of order k=$(order(t)) on [$(first(t)),$(last(t))] ($(numintervals(t)) intervals)")
end

# Knot set with arbitrary knots
"""
    struct KnotSet{T} <: AbstractKnotSet{T}

Represents a knot set that has all the knots explicitly specified.
"""
struct KnotSet{T} <: AbstractKnotSet{T}
    order :: Int
    ts :: Vector{T}

    function KnotSet(order::Integer, ts::AbstractVector{<:Real})
        length(ts) >= order + 1 || throw(DomainError(order, "Too few knots ($(length(ts))) for KnotSet of order $order"))
        new{eltype(ts)}(order, sort(ts))
    end
end
order(t::KnotSet) = t.order
Base.length(t::KnotSet) = length(t.ts)
Base.getindex(t::KnotSet, i::Integer) = t.ts[i]


struct LinearKnotSet{T,R <: AbstractRange{T}} <: AbstractKnotSet{T}
    k :: Integer
    ts :: R

    function LinearKnotSet(k::Integer, a::T, b::T, N::Integer) where {T <: Real}
        r = range(a, stop=b, length=N)
        new{eltype(r),typeof(r)}(k, r)
    end
end
order(t::LinearKnotSet) = t.k
Base.length(t::LinearKnotSet) = length(t.ts) + 2*t.k - 2
function Base.getindex(t::LinearKnotSet, i::Integer)
    if 1 <= i <= t.k
        first(t.ts)
    elseif t.k < i < t.k + length(t.ts)
        t.ts[i - t.k + 1]
    elseif t.k + length(t.ts) <= i <= 2*t.k + length(t.ts) - 2
        last(t.ts)
    else
        throw(BoundsError(t, i))
    end
end
knotsetname(t::LinearKnotSet) = "LinearKnotSet{$(eltype(t))}"

# struct ExpKnotSet{T} <: AbstractKnotSet{T}
#     k::Integer
#     exponents::AbstractRange{T}
#     base::T
#     t::AbstractVector{T}
#     include0::Bool
# end
# function ExpKnotSet(k::Integer, a::T, b::T, N::Integer;
#                     base::T=T(10), include0::Bool=true) where T
#     exponents = range(a, stop=b, length=include0 ? N : N+1)
#     t = base .^ exponents
#     ExpKnotSet(k, exponents, eltype(t)(base), include0 ? vcat(0,t) : t, include0)
# end
#
# function Base.getproperty(ks::ExpKnotSet, property::Symbol)
#     if property === :ts
#         getfield(ks, :t)
#     else
#         getfield(ks, property)
#     end
# end
#
#
# function show(io::IO, t::ExpKnotSet)
#     write(io, "$(typeof(t)) of order k=$(order(t)) on [")
#     t.include0 && write(io, "0,")
#     write(io, "$(t.base^first(t.exponents))..$(t.base^last(t.exponents))] ($(numintervals(t)) intervals)")
# end

# function arcsin_knot_set(k::Integer, a::Integer, b::Integer, N::Integer)
#     N2 = N/2
#     arcsin.(range(-N2, stop=N2, lengthN+1)/N2)*(b-a)/π+(a+b)/2
# end

# arcsin_half_knot_set(a, b, N) = arcsin(linspace(0,N,N+1)/N)*(b-a)*2.0/π+a

# function exp_linear_knot_set(a, b, N)
#     ap = a != 0 ? a : 1e-1
#     @assert t[2][1]+t[2][2] == N
#     [logspace(log10(ap),log10(t[1]), t[2][1]); linspace(t[1],b,t[2][2]+1)[2:end]]
# end

@recipe function plot(t::AbstractKnotSet)
    markershape --> :circle
    y = similar(t)
    y′ = 1
    t′ = -Inf
    for i in eachindex(t)
        if t[i] == t′
            y′ += 1
        else
            y′ = 1
        end
        t′ = t[i]
        y[i] = y′
    end
    collect(t),y
end

export KnotSet, LinearKnotSet, #=ExpKnotSet,=# order, numintervals
