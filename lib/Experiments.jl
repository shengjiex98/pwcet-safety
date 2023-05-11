
module Experiments

export SamplerPWCET
export single_run_deviation, binomial_prob, lr_test, lr_test_2
export sim, missrow, misstotal, missfirst

import Random
using Distributions
using RealTimeScheduling
using ControlTimingSafety

# struct SamplerPWCET <: Random.Sampler{BitVector}
struct SamplerPWCET <: SamplerWeaklyHard
    p::Real
    H::Integer
end

function Random.rand!(rng::Random.AbstractRNG, a::BitVector, sp::SamplerPWCET)
    for i = 1:sp.H
        a[i] = Random.rand(rng) < sp.p
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

function binomial_prob(n::Integer, p::Real, x::Integer)
    @boundscheck 0 <= p <= 1 || throw(ArgumentError("p has to be within 0 and 1"))
    binomial(n, x) * p^x * (1-p)^(n-x)
end

function lr_test(θ::Real, n::Integer, x::Integer)
    @boundscheck 0 <= θ <= 1 || throw(ArgumentError("theta has to be within 0 and 1"))
    # Calculate the θ value from observed data
    observed = x / n
    θ0 = min(observed, θ)
    (θ0^x * (1-θ0)^(n-x)) / (observed^x * (1-observed)^(n-x))
end

function lr_test_2(θ::Real, n::Integer, ϵ::Real)
    @boundscheck 0 <= θ <= 1 || throw(ArgumentError("theta has to be within 0 and 1"))
    # Calculate the θ value from observed data
    x = round(Int64, n * θ)
    observed = x / n
    θ0 = [θ-ϵ, θ+ϵ]
    maximum(@. (θ0/observed)^x * ((1-θ0)/(1-observed))^(n-x))
end

function sim(a::Automaton, z0::AbstractVector{<:Real}, p::Real, n::Integer; H::Integer=100)
    sp = SamplerPWCET(p, H)
    samples = Vector{Tuple{BitVector, Float64}}(undef, n)
    for i in 1:n
        σ = rand(sp)
        samples[i] = (σ, single_run_deviation(a, z0, 2 .- σ))
    end
    sort!(samples, by=x -> x[2])
end

function missrow(σ::BitVector)
    counter = maxcount = 0
    for σ_i in σ
        if σ_i == 0
            counter += 1
            maxcount = max(maxcount, counter)
        else
            counter = 0
        end
    end
    maxcount
end

function misstotal(σ::BitVector)
    reduce((a, b) -> a + (1 - b), σ; init=0)
end

function missfirst(σ::BitVector)
    f = findfirst(σ .== 0)
    f === nothing ? 101 : f
end


# function visualize(a::Automaton, z_0:Vector{<:Real}, bv::BitVector)
#     x = length(bv)
#     y = deviation(a, z_0, )
# end

end
