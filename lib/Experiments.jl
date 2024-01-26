
module Experiments

export SamplerPWCET
export single_run_deviation, generate_samples, binomial_prob, lr_test, lr_test_2
export find_intervals

import Random
using Distributions
using RealTimeScheduling
using ControlTimingSafety
using Base.Threads

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

function single_run_deviation(a::Automaton, z_0::AbstractVector{Float64}, input::AbstractVector{Int64})
    @boundscheck length(z_0) == a.nz || throw(DimensionMismatch("z_0 must have length a.nz"))

    # Dimensions: time, state
    z = evol(a, z_0, input)
    # Dimensions: time, state, min/max
    reachable = cat(z, z, dims=3)

    maximum(deviation(a, z_0, reachable))
end

function generate_samples(a::Automaton, z0::AbstractVector{<:Real}, q::Real, n::Integer; H::Integer=100)
    sp = SamplerPWCET(q, H)
    samples = Vector{Tuple{BitVector,Float64}}(undef, n)
    Threads.@threads for i in 1:n
        σ = rand(sp)
        samples[i] = (σ, single_run_deviation(a, z0, 2 .- σ))
    end
    sort!(samples, by=x -> x[2])
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
    find_intervals(n, p, α)

Find optimal intervals for a given confidence value (expressed by 1-α)
"""
function find_intervals(n::Integer, p::Real, α::Real, fullresults=false)
    @boundscheck 0 <= p <= 1 || throw(ArgumentError("p has to be within 0 and 1"))
    @boundscheck 0 <= α <= 1 || throw(ArgumentError("α has to be within 0 and 1"))

    dist = Binomial(n, p)

    
    i2s = fill(-1, n-1)
    i2  = 2

    # Iterate over i1
    for i1 in 1:n-1
        # Find first suitable i2
        while cdf(dist, i1-1) + (1-cdf(dist, i2-1)) > α && i2 <= n
            i2 += 1
        end
        # Break if i2 goes out of range
        if cdf(dist, i1-1) + (1-cdf(dist, i2-1)) <= α
            i2s[i1] = i2
        else
            # @info "Loop ends at" i1
            break
        end
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

end
