module BSplines

type BSpline
    α::AbstractVector
    k::Integer
end

include("basis.jl")
include("knot_sets.jl")

export BSpline

end
