module Experiments

export SamplerPWCET
export single_run_deviation, generate_samples, generate_filename, generate_samples_mpc, generate_samples_mpc_with_multi_ref
export binomial_prob, lr_test, lr_test_2
export find_intervals

import Random
using Distributions
using Memoization
using Printf
using ControlSystemsBase

using RealTimeScheduling
using ControlTimingSafety
using MATLABControlTest

# struct SamplerPWCET <: Random.Sampler{BitVector}
struct SamplerPWCET <: SamplerWeaklyHard
    q::Real
    H::Integer
end

function Random.rand!(rng::Random.AbstractRNG, a::BitVector, sp::SamplerPWCET)
    for i = 1:sp.H
        a[i] = Random.rand(rng) < sp.q
    end
    a
end

function Random.rand(rng::Random.AbstractRNG, sp::SamplerPWCET)
    a = falses(sp.H)
    Random.rand!(rng, a, sp)
end

function single_run_deviation(a::Automaton, z_0::AbstractVector{Float64}, input::AbstractVector{Int64};
        nominal_trajectory=nothing)
    @boundscheck length(z_0) == a.nz || throw(DimensionMismatch("z_0 must have length a.nz"))
    @boundscheck nominal_trajectory === nothing || 
        size(nominal_trajectory, 1) == a.nz || throw(DimensionMismatch("First dim of nominal_trajectory must be a.nz"))
    @boundscheck nominal_trajectory === nothing || 
        size(nominal_trajectory, 2) == length(input) + 1 || throw(DimensionMismatch("Second dim of nominal_trajectory must be input length+1"))

    # Dimensions: time, state
    z = evol(a, z_0, input)
    # Dimensions: time, state, corners (just 1)
    reachable = reshape(z, size(z)..., 1)

    if nominal_trajectory !== nothing
        nominal_trajectory = reshape(a.C * nominal_trajectory, size(a.C, 1), 1, length(input) + 1)
        maximum(deviation(a, z_0, reachable, nominal_trajectory=nominal_trajectory))
    else
        maximum(deviation(a, z_0, reachable))
    end
end

function generate_samples(a::Automaton, z0::AbstractVector{<:Real}, q::Real, n::Integer;
        H::Integer=100, sorted=true, nominal_trajectory=nothing)
    sp = SamplerPWCET(q, H)
    samples = Vector{Tuple{BitVector,Float64}}(undef, n)
    Threads.@threads for i in 1:n
        σ = rand(sp)
        samples[i] = (σ, single_run_deviation(a, z0, 2 .- σ, nominal_trajectory=nominal_trajectory))
    end
    if sorted
        sort!(samples, by=x -> x[2])
    else
        samples
    end
end

function generate_samples_mpc(sysd::AbstractStateSpace{<:ControlSystemsBase.Discrete},
        x0::AbstractVector{<:Real}, refs::AbstractVector{<:Real},
        q::Real, n::Integer; H::Integer=100, sorted=true)
    sp = SamplerPWCET(q, H)
    samples = Vector{Tuple{BitVector,Float64}}(undef, n)
    for i in 1:n
        σ = rand(sp)
        d = maximum(abs.(run_simulation(sysd, x0, refs; input = σ) - refs))
        samples[i] = (σ, d)
    end
    if sorted
        sort!(samples, by=x -> x[2])
    else
        samples
    end
end

function generate_samples_mpc_with_multi_ref(sysd::AbstractStateSpace{<:ControlSystemsBase.Discrete},
        x0::AbstractVector{<:Real}, refs1::AbstractVector{<:Real}, refs2::AbstractVector{<:Real},
        q::Real, n::Integer; H::Integer=100, sorted=true)
    sp = SamplerPWCET(q, H)
    samples = Vector{Tuple{BitVector,Float64,Float64}}(undef, n)
    for i in 1:n
        σ = rand(sp)
        d1 = maximum(abs.(run_simulation(sysd, x0, refs1; input = σ) - refs1))
        d2 = maximum(abs.(run_simulation(sysd, x0, refs2; input = σ) - refs2))
        samples[i] = (σ, d1, d2)
    end
    if sorted
        sort!(samples, by=x -> x[2])
    else
        samples
    end
end


function generate_filename(batchsize::Integer, q::Real, h::Real, n::Integer; th::Integer=Threads.nthreads())
    @sprintf "b%.1e-q%.9g-h%.9g-n%i-th%i" batchsize q h n th
end

function generate_filename(batchsize::Integer, q::Real, h::Real; th::Integer=Threads.nthreads())
    @sprintf "b%.1e-q%.9g-h%.9g-th%i" batchsize q h th
end

"""
    binomial_prob(n, p, x)

Calculates the probability that `n` i.i.d. samples has `x` 1s, given that each sample has probability `p` of being 1.
"""
function binomial_prob(n::Integer, p::Real, x::Integer)
    @boundscheck 0 <= p <= 1 || throw(ArgumentError("p has to be within 0 and 1"))
    binomial(BigInt(n), BigInt(x)) * p^x * (1 - p)^(n - x)
end

function lr_test(θ::Real, n::Integer, x::Integer)
    @boundscheck 0 <= θ <= 1 || throw(ArgumentError("θ has to be within 0 and 1"))
    # Calculate the θ value from observed data
    observed = x / n
    θ0 = min(observed, θ)
    (θ0^x * (1 - θ0)^(n - x)) / (observed^x * (1 - observed)^(n - x))
end

function lr_test_2(θ::Real, n::Integer, ϵ::Real)
    @boundscheck 0 <= θ <= 1 || throw(ArgumentError("θ has to be within 0 and 1"))
    # Calculate the θ value from observed data
    x = round(Int64, n * θ)
    observed = x / n
    θ0 = [θ - ϵ, θ + ϵ]
    maximum(@. (θ0 / observed)^x * ((1 - θ0) / (1 - observed))^(n - x))
end

"""
    find_intervals(n, p, α; fullresults=false, centered=false)

Find confidence intervals (expressed by 1-α) of the p-th quantile for n samples. If 
fullresults is set to true, include dominated intervals (e.g., [5, 9] is dominated by [5, 7]);
if centered is set to true, only return the interval with the p-th quantile at the center.
true, 
"""
function find_intervals(n::Integer, p::Real, α::Real; centered=false, fullresults=false)
    @boundscheck 0 <= p <= 1 || throw(ArgumentError("p has to be within 0 and 1"))
    @boundscheck 0 <= α <= 1 || throw(ArgumentError("α has to be within 0 and 1"))

    dist = Binomial(n, p)
    # Use Memoization.jl to cache calculated cdf results (saves time for large n values).
    @memoize cdf_cached(i) = cdf(dist, i)

    if centered
        c = round(Int64, n * p)
        prob_mass = (1 - α - pdf(dist, c)) / 2

        i1 = i2 = c
        while cdf_cached(i2) - cdf_cached(c) < prob_mass
            i2 += 1
        end
        while cdf_cached(c-1) - cdf_cached(i1-1) < prob_mass
            i1 -= 1
        end
        return [(i1, i2)]
    end
    
    # i2 values for each i1 in 1:n, initialized to -1
    i2, i2s = 1, fill(-1, n)
    for i1 in 1:n
        # Find first suitable i2
        while cdf_cached(i1-1) + (1-cdf_cached(i2)) > α && i2 < n
            i2 += 1
        end
        if cdf_cached(i1-1) + (1-cdf_cached(i2)) > α
            break
        end
        i2s[i1] = i2
    end
    if fullresults
        return [(i1, i2s[i1]) for i1 in 1:n-1]
    end

    # Remove suboptimal results (same i2 only keep highest i1)
    i1 = 1
    i2 = i2s[1]
    results = []
    while i2 != -1
        while i2s[i1+1] == i2
            i1 += 1
        end
        push!(results, (i1, i2))
        i1 += 1
        i2 = i2s[i1]
    end

    results
end

function constraint_satisfiability(c::WeaklyHardConstraint, q::Real, H::Integer)
    # return Pr(σ ⊢ c | σ ∈ {0, 1}^H and Pr(σ[i] = 1) = q)
end

end
