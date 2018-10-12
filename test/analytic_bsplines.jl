# The reference is generated thusly:
using SymPy
using ProgressMeter
function gen_basis(t)
    @syms x
    d(a,b) = b == 0 ? 0 : a/b
    interval(x,a,b) = Gt(x,a)∧Lt(x,b)
    interval(x,a::Rational,b::Rational) = Gt(x*denominator(a),numerator(a))∧Lt(x*denominator(b),numerator(b))
    B = [[piecewise((1,interval(x,t[j],t[j+1])),(0,true))
          for j = 1:length(t)-1]]
    for k = 2:order(t)
        Bₖ = [d((x-t[i]),(t[i+k-1]-t[i]))*B[end][i] +
              d((t[i+k]-x),(t[i+k]-t[i+1]))*B[end][i+1]
              for i = 1:length(t)-k]
        push!(B,Bₖ)
    end
    B
end
function scalar_op(Bₖ, t, f=x->1)
    B = Bₖ[end]
    n = order(t)+1
    O = Matrix{Any}(undef, n,n)
    @syms x
    l(v::Rational) = denominator(v) == 1 ? numerator(v) : 1.0*v
    l(v) = v
    @showprogress for i = 1:n
        for j = 1:n
            O[i,j] = integrate(B[i]*f(x)*B[j], (x,l(first(t)),l(last(t))))
        end
    end
    O
end
