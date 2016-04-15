linear_knot_set(a, b, N) = linspace(a, b, N+1)

function arcsin_knot_set(a, b, N)
    N2 = N/2
    arcsin(linspace(-N2, N2, N+1)/N2)*(b-a)/π+(a+b)/2
end

arcsin_half_knot_set(a, b, N) = arcsin(linspace(0,N,N+1)/N)*(b-a)*2.0/π+a

function exp_knot_set(a, b, N)
    ap = a != 0 ? a : 1e-1
    logspace(log10(ap),log10(b),N)
end

function exp_linear_knot_set(a, b, N)
    ap = a != 0 ? a : 1e-1
    @assert t[2][1]+t[2][2] == N
    [logspace(log10(ap),log10(t[1]), t[2][1]); linspace(t[1],b,t[2][2]+1)[2:end]]
end

export linear_knot_set, arcsin_knot_set, arcsin_half_knot_set, exp_knot_set, exp_linear_knot_set
