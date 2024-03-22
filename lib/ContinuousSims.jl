module ContinuousSims

using ControlSystems: StateSpace, Continuous, lsim
using LinearAlgebra: I

export nominal_trajectory

function nominal_trajectory(sys::StateSpace{<:Continuous}, u::Function, period::Real, H::Real, x0::Vector{<:Real})
    @boundscheck length(x0) == sys.nx || throw(DimensionMismatch("x0 must have same length as sys.nx"))
    t = 0:period:H
    y, t, x, uout = lsim(sys, u, t, x0=x0)
    [x; uout]
end

end
