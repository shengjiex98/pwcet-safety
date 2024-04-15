module Benchmarks

export benchmarks
export delay_lqr, pole_place, augment

using ControlSystemsBase
using LinearAlgebra

# Resistor-capacitor network
sys_rcn = let
	r_1 = 100000
	r_2 = 500000
	r_3 = 200000
	c_1 = 0.000002
	c_2 = 0.000010
	A = [-1/c_1 * (1/r_1 + 1/r_2)  1/(r_2*c_1)
	     1/(r_2*c_2)               -1/c_2 * (1/r_2 + 1/r_3)]
	B = [1/(r_1*c_1)
	     1/(r_3*c_2)]
	# C = [1  -1]
    C = I
	D = 0

	ss(A, B, C, D)
end

# F1-tenth car
sys_f1t = let 
    v = 6.5
    L = 0.3302
    d = 1.5
    A = [0 v ; 0 0]
    B = [0; v/L]
    C = [1 0]
    D = 0

    ss(A, B, C, D)
end

# DC motor
sys_dcm = let
    A = [-10 1; -0.02 -2]
    B = [0; 2]
    C = [1 0]
    D = 0

    ss(A, B, C, D)
end

# Car suspension system
sys_css = let
    A = [0. 1 0 0; -8 -4 8 4; 0 0 0 1; 80 40 -160 -60]
    B = [0; 80; 20; -1120]
    C = [1 0 0 0]
    D = 0

    ss(A, B, C, D)
end

# Electronic wedge brake
sys_ewb = let
    A = [0 1; 8.3951e3 0]
    B = [0; 4.0451]
    C = [7.9920e3 0]
    D = 0

    ss(A, B, C, D)
end

# Cruise control 1
sys_cc1 = let
    A = -0.05
    B = 0.01
    C = 1
    D = 0

    ss(A, B, C, D)
end

# Cruise control 2
sys_cc2 = let
    A = [0 1 0; 0 0 1; -6.0476 -5.2856 -0.238]
    B = [0; 0; 2.4767]
    C = [1 0 0]
    D = 0

    ss(A, B, C, D)
end

sys_mpc = ss(tf([3, 1],[1, 0.6, 1]))

# TODOs: figure out why EWB doesn't work, for ControlTimingSafety package.
benchmarks = Dict([
    :RCN => sys_rcn,
    :F1T => sys_f1t,
    :DCM => sys_dcm,
    :CSS => sys_css,
    :EWB => sys_ewb,
    :CC1 => sys_cc1,
    :CC2 => sys_cc2,
    :MPC => sys_mpc
])

function delay_lqr(sys::AbstractStateSpace{<:Continuous}, h::Float64; Q=I, R=I)
    sysd_delay = c2da(sys, h)
    lqr(sysd_delay, Q, R)
end

# Pole placement only works with single input single output (SISO) systems. i.e., 
# sys.B & sys.D can only have 1 column and sys.C & sys.D can only have 1 row
function pole_place(sys::AbstractStateSpace{<:Continuous}, h::Float64; p=0.9)
    sysd_delay = c2da(sys, h)
    place(sysd_delay, vcat([0], fill(p, size(sys.A)[1])))
end

function c2da(sys::AbstractStateSpace{<:Continuous}, h::Float64)
    augment(c2d(sys, h))
end

function augment(sysd::AbstractStateSpace{<:Discrete})
    p = size(sysd.A, 1)
    q = size(sysd.B, 2)
    r = size(sysd.C, 1)
    A = [sysd.A sysd.B; zeros(q, p+q)]
    B = [zeros(p, q); I]
    C = [sysd.C zeros(r, q)]
    D = sysd.D
    ss(A, B, C, D, sysd.Ts)
end

end