struct Basis
    k::Integer
    t::AbstractVector
    x::AbstractVector
    B::Vector{Matrix}
end

# See http://pages.cs.wisc.edu/~deboor/pgs/bsplvb.f
function Basis(k::Integer, t::AbstractVector, x::AbstractVector)
    basis = [zeros(eltype(x), (length(x),length(t)-kk)) for kk = 1:k]

    #=
    \[B_{i,1,t}(x)=\begin{cases}
    1, & t_i\leq x < t_{i+1}\\
    0, & \textrm{else}
    \end{cases}\]
    =#
    for i = 1:length(t)-1
        basis[1][:,i] = (t[i] .<= x .< t[i+1])
    end

    #=
    \[B_{i,k,t}(x)=\frac{x-t_i}{t_{i+k-1}-t_i}B_{i,k-1,t}(x)-
    \frac{t_{i+k}-x}{t_{i+k}-t_{i+1}}B_{i+1,k-1,t}(x)\]
    =#
    for kk = 2:k
        for i = 1:length(t) - kk
            d₁ = t[i+kk-1]-t[i]
            d₂ = t[i+kk]-t[i+1]
            f₁ = d₁ == 0 ? 0 : (x-t[i])/d₁
            f₂ = d₂ == 0 ? 0 : (t[i+kk]-x)/d₂
            basis[kk][:,i] = f₁.*basis[kk-1][:,i] + f₂.*basis[kk-1][:,i+1]
        end
    end

    Basis(k, t, x, basis)
end

#=
\[f(x)=\sum_i\alpha_iB_{i,k,t}(x)\equiv B(x)\vec{\alpha} = \langle B|\vec{\alpha}\rangle.\]
=#

function (B::Basis)(S::Spline)
    if S.k > B.k
        error("No support for splines of order $(S.k) in basis of order $(B.k)")
    end
    B.B[S.k]*S.α
end
