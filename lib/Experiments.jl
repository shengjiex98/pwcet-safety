
module Experiments

export SamplerPWCET
export single_run_deviation, binomial_prob, lr_test, lr_test_2
export sim, missrow, misstotal, missfirst
export ev_consc_miss_prep

import Random
using Distributions
using RealTimeScheduling
using ControlTimingSafety
using OffsetArrays

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

using OffsetArrays
function ev_consc_miss_prep(n, p)
    T = OffsetArray(zeros(n+1, n+1), 0:n, 0:n)
    C = OffsetArray(zeros(n+1, n+1), 0:n, 0:n)
    
    # Initializing
    for i in 0:n
        # No misses
        T[i, 0] = p^i
        C[i, 0] = T[i, 0]
        
        # All misses
        T[i, i] = (1-p)^i
        C[i, i] = 1
        
        # Convenience shortcuts
        for j in i+1:n
            T[i, j] = 0
            C[i, j] = 1
        end
    end
    
    for i in 2:n
        ## This line is invalid: a sequence can have multiple separated misses
        # T[i, 1]   = i * p^(i-1) * (1-p)
        ## This line can be covered by the loop below
        # T[i, i-1] = 2 * p * (1-p)^(i-1)
        # Assuming the **first** longest sequence starts at k and has length j
        for j in 1:i-1, k in 1:i-j+1
            if k == 1
                T[i, j] += T[j, j] * p * C[i-j-1, j]
            elseif k == (i-j+1)
                # Since the sequence is at the end, earlier consecutive misses 
                # are at most j-1 in length.
                T[i, j] += C[i-j-1, j-1] * p * T[j, j]
            else
                T[i, j] += C[k-2, j-1] * p * T[j, j] * p * C[i-k-j, j]
            end
        end
        for j in 1:i
            C[i, j] = C[i, j-1] + T[i, j]
        end
    end
    
    return T, C
end

# function visualize(a::Automaton, z_0:Vector{<:Real}, bv::BitVector)
#     x = length(bv)
#     y = deviation(a, z_0, )
# end

end
