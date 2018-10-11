module BSplines

mutable struct Spline
    α::AbstractVector
    k::Integer
end

include("basis.jl")
include("knot_sets.jl")

end
