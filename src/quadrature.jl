using FastGaussQuadrature

#=
Gaußian quadrature:

\[\int\limits_a^b dx\;f(x)\approx
\frac{b-a}{2}\sum_{i=1}^n w_i
f\left(\frac{b-a}{2}x_i+\frac{a+b}{2}\right)\]
=#
function lgwt!(x::AbstractVector{T}, w::AbstractVector{T},
               xs::AbstractVector{T}, ws::AbstractVector{T},
               i,
               a::T=zero(T), b::T=one(T)) where {T<:Real}
    n = length(x)
    interval = (i-1)*n+1 : i*n
    xs[interval] = 0.5((b-a)*x .+ a .+ b)
    ws[interval] = 0.5(b-a)*w
    nothing
end

function lgwt(t::AbstractKnotSet{T}) where T
    k = order(t)
    N = numintervals(t)
    @warn "Non-standard k"
    x, w = gausslegendre(k+20)
    xo = zeros(T, N*length(x))
    wo = zeros(T, N*length(x))

    for i = 1:N
        lgwt!(x, w,
              xo, wo, i,
              t[k+i-1], t[k+i])
    end

    xo,wo
end
