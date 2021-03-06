import Base: first, last, length,
    getindex, lastindex, eachindex, iterate,
    eltype, similar,
    show
abstract type AbstractKnotSet{T} end

order(t::AbstractKnotSet) = t.k
numintervals(t::AbstractKnotSet) = length(t.t)-1

first(t::AbstractKnotSet) = first(t.t)
last(t::AbstractKnotSet) = last(t.t)
length(t::AbstractKnotSet) = length(t.t) + 2t.k - 2

function getindex(t::AbstractKnotSet, i::Integer)
    if i < order(t)
        first(t)
    elseif i < order(t) + numintervals(t)
        t.t[i-order(t)+1]
    else
        last(t)
    end
end
lastindex(t::AbstractKnotSet) = length(t)
eachindex(t::AbstractKnotSet) = 1:length(t)
function iterate(t::AbstractKnotSet, state=(t[1],0))
    element, i = state
    if i >= length(t)
        return nothing
    end
    element, (t[i+2], i+1)
end

eltype(t::AbstractKnotSet{T}) where T = T
similar(t::AbstractKnotSet{T}) where T = Vector{T}(undef, length(t))

function show(io::IO, t::AbstractKnotSet)
    write(io, "$(typeof(t)) of order k=$(order(t)) on [$(first(t)),$(last(t))] ($(numintervals(t)) intervals)")
end

struct LinearKnotSet{T} <: AbstractKnotSet{T}
    k::Integer
    t::AbstractRange{T}
end
LinearKnotSet(k::Integer, a::Integer, b::Integer, N::Integer) =
    LinearKnotSet(k, range(a, stop=b, length=N+1))

struct ExpKnotSet{T} <: AbstractKnotSet{T}
    k::Integer
    exponents::AbstractRange{T}
    base::T
    t::AbstractVector{T}
    include0::Bool
end
function ExpKnotSet(k::Integer, a::T, b::T, N::Integer;
                    base::T=T(10), include0::Bool=true) where T
    exponents = range(a, stop=b, length=include0 ? N : N+1)
    t = base .^ exponents
    ExpKnotSet(k, exponents, eltype(t)(base), include0 ? vcat(0,t) : t, include0)
end

function show(io::IO, t::ExpKnotSet)
    write(io, "$(typeof(t)) of order k=$(order(t)) on [")
    t.include0 && write(io, "0,")
    write(io, "$(t.base^first(t.exponents))..$(t.base^last(t.exponents))] ($(numintervals(t)) intervals)")
end

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

export LinearKnotSet, ExpKnotSet, order, numintervals
