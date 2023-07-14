
module Experiments

export SamplerPWCET
export single_run_deviation, binomial_prob, lr_test, lr_test_2
export sim, missrow, misstotal, missfirst, calculate_mean_miss
export beta, Fcv, hello, FcvW, inverse_fcv, Psi2_cv, Tau2, z_alpha, confidence_interval

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

"""
    calculate_mean_miss(V::Function, samples::Vector{Tuple{BitVector, Float64}})

    This function calulates the mean value of the Control Variates
"""
function calculate_mean_miss(V::Function, samples::Vector{Tuple{BitVector, Float64}})
    total_miss = 0
    num_samples = length(samples)

    for (σ, _) in samples
        miss_value = 0
        miss_value = V(σ)
        total_miss += miss_value
    end
    
    mean_miss = total_miss / num_samples
    return mean_miss
end

"""
    beta(y::Real, samples::Vector{Tuple{BitVector, Float64}}, V::Function)

    This function calculates the Control Variates estimator
"""
function beta(y::Real, samples::Vector{Tuple{BitVector, Float64}}, V::Function)
    n = length(samples)
    sum_1 = 0;
    sum_2 = 0;
    sum_3 = 0;
    mean = calculate_mean_miss(V, samples)
    for (σ, devation) in samples
        sum_1 += indicator(devation, y) * V(σ)
        sum_2 += indicator(devation, y)
        sum_3 += (V(σ) - mean) ^ 2
    end
    beta = ((1/n) * sum_1 - (1/n) * sum_2 * mean)/((1/n) * sum_3)
    return beta
end

"""
    indicator(Y::Real, y::Real)

    This is a function for the indicator function, returns 1 when Y <=y, returns 0 otherwise
"""
function indicator(Y::Real, y::Real)
    if Y <= y
        return 1
    else
        return 0
    end
end

"""
    Fcv(y::Real, samples::Vector{Tuple{BitVector, Float64}}, V::Function, mean::Real)

    This function calulates the CDF for control variates estimator beta
"""
function Fcv(y::Real, samples::Vector{Tuple{BitVector, Float64}}, V::Function, mean::Real)
    n = length(samples)
    sum = 0
    for (σ, devation) in samples
        sum += indicator(devation, y)
    end
    nmc = (1/n) * sum
    b = beta(y, samples, V)
    fcv = nmc - b * (calculate_mean_miss(V, samples)-mean)
    return fcv
end

"""
    FcvW(y::Real, m::Real, samples::Vector{Tuple{BitVector, Float64}}, v::Function)

    This function calulates the CDF for control variates with estimator W
"""

function FcvW(y::Real, m::Real, samples::Vector{Tuple{BitVector, Float64}}, v::Function)
    n = length(samples)
    sum_1 = 0
    Fcv = 0
    mean = calculate_mean_miss(v, samples)
    for (σ, _) in samples
        sum_1 += (v(σ) - mean) ^ 2
    end
    for (σ, devation) in samples
        Wi = 1/n + (mean - v(σ)) * (mean - m)/sum_1
        Fcv += Wi * indicator(devation, y)
    end
    return Fcv
end

"""
    inverse_fcv(p::Real, m::Real, samples::Vector{Tuple{BitVector, Float64}}, v::Function)

    This function calulates the inverse function of CV CDF
"""

function inverse_fcv(p::Real, m::Real, samples::Vector{Tuple{BitVector, Float64}}, v::Function)
    n = length(samples)
    sum_1 = 0
    sum_2 = 0
    y = 0
    mean = calculate_mean_miss(v, samples)
    for (σ, _) in samples
        sum_1 += (v(σ) - mean) ^ 2
    end
    for (σ, devation) in samples
        Wi = 1/n + (mean - v(σ)) * (mean - m)/sum_1
        sum_2 += Wi
        if sum_2 >= p
            y = devation
            break
        end
    end
    return y
end

function var(samples::Vector{Tuple{BitVector, Float64}},v::Function)
    n = length(samples)
    ret = 0
    sum_1 = 0
    mean = calculate_mean_miss(v, samples)
    for (σ, _) in samples
        sum_1 += (v(σ)- mean) ^ 2
    end
    ret = 1/(n-1) * sum_1
    return ret
end

function cov(quantile::Real,samples::Vector{Tuple{BitVector, Float64}},v::Function)
    n = length(samples)
    ret = 0
    sum_1 = 0
    sum_2 = 0
    mean = calculate_mean_miss(v, samples)
    for (_, deviation) in samples
        sum_2 += indicator(deviation, quantile)
    end
    mean2 = sum_2 / n
    for (σ, deviation) in samples
        sum_1 += (v(σ)- mean) * (deviation - mean2)
    end
    ret = 1/(n-1) * sum_1
    return ret
end

function Psi2_cv(quantile::Real,p::Real,samples::Vector{Tuple{BitVector, Float64}},v::Function)
    Psi = p*(1-p) - ((cov(quantile,samples,v))^2)/var(samples,v)
    return Psi
end

function Eta(Delta::Real, p::Real, m::Real, samples::Vector{Tuple{BitVector, Float64}}, v::Function)
    Eta = (inverse_fcv(p + Delta, m, samples, v) - inverse_fcv(p - Delta, m, samples, v))/(2 * Delta)
    return Eta
end

function Tau2(quantile::Real,Delta::Real,p::Real,m::Real,samples::Vector{Tuple{BitVector, Float64}},v::Function)
    Tau = Psi2_cv(quantile, p, samples, v) * (Eta(Delta, p, m, samples, v))^2
    return Tau
end

function z_alpha(alpha::Real)
    d = Normal(0.0,1.0)
    z_alpha = quantile.(d,[1-alpha/2])
    return z_alpha[1]
end

function confidence_interval(alpha::Real, quantile::Real,Delta::Real,p::Real,m::Real,samples::Vector{Tuple{BitVector, Float64}},v::Function)
    n = length(samples)
    diff = z_alpha(alpha) * sqrt(Tau2(quantile, Delta, p, m, samples, v))/sqrt(n)
    i1 = quantile - diff
    i2 = quantile + diff
    return [i1, i2, diff]
end


# function visualize(a::Automaton, z_0:Vector{<:Real}, bv::BitVector)
#     x = length(bv)
#     y = deviation(a, z_0, )
# end

end
