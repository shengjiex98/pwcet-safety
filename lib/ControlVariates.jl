module ControlVariates

export missrow, misstotal, missfirst
export calculate_mean_cv, calculate_missrow_prob
export β, fcv, fcvW, inverse_fcv, ψ2_cv, τ2, z_α, confidence_interval, bin_list

using OffsetArrays

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

function missrow(σ::BitVector, n::Integer)
    k = σ[1:n]
    missrow(k)
end

function misstotal(σ::BitVector)
    reduce((a, b) -> a + (1 - b), σ; init=0)
end

function misstotal(σ::BitVector, n::Integer)
    k = σ[1:n]
    misstotal(k)
end

function missfirst(σ::BitVector)
    f = findfirst(σ .== 0)
    f === nothing ? length(σ)+1 : f
end

function missfirst(σ::BitVector, n::Integer)
    k = σ[1:n]
    missfirst(k)
end

"""
    calculate_mean_cv(v, samples)

Calulate the mean value of the control variates with variate function v and given samples.
"""
function calculate_mean_cv(v::Function, samples::Vector{Tuple{BitVector,Float64}})
    ## The following section can be simplified by using a map
    # total_miss = 0
    # num_samples = length(samples)
    # for (σ, _) in samples
    #     miss_value = v(σ)
    #     total_miss += miss_value
    # end
    # mean_miss = total_miss / num_samples

    miss_values = map((σ, _) -> v(σ), samples)
    mean_miss = mean(miss_values)

    return mean_miss
end

function calculate_missrow_prob(n, p)
    T = OffsetArray(zeros(n + 1, n + 1), 0:n, 0:n)
    C = OffsetArray(zeros(n + 1, n + 1), 0:n, 0:n)

    # Initializing
    for i in 0:n
        # No misses
        T[i, 0] = p^i
        C[i, 0] = T[i, 0]

        # All misses
        T[i, i] = (1 - p)^i
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
            elseif k == (i - j + 1)
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

"""
    β(y, samples, v)

Calculate the control variates estimator β with given samples, 
control variate function v and desired deviation y.
"""
function β(y::Real, samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    n = length(samples)
    sum_1 = 0
    sum_2 = 0
    sum_3 = 0
    mean = calculate_mean_cv(v, samples)
    for (σ, devation) in samples
        sum_1 += indicator(devation, y) * v(σ)
        sum_2 += indicator(devation, y)
        sum_3 += (v(σ) - mean)^2
    end
    β = ((1 / n) * sum_1 - (1 / n) * sum_2 * mean) / ((1 / n) * sum_3)
    return β
end

"""
    indicator(Y, y)

Indicator function, returns 1 when Y <= y, returns 0 otherwise.
"""
function indicator(Y::Real, y::Real)
    return if Y <= y 1 else 0 end
end

"""
    fcv(y, samples, v)

Calulate the CDF for control variates estimator β.
"""
function fcv(y::Real, samples::Vector{Tuple{BitVector,Float64}}, v::Function, mean::Real)
    n = length(samples)
    sum = 0
    for (σ, devation) in samples
        sum += indicator(devation, y)
    end
    nmc = (1 / n) * sum
    β = β(y, samples, v)
    fcv = nmc - β * (calculate_mean_cv(v, samples) - mean)
    return fcv
end

"""
    fcvW(y::Real, m::Real, samples::Vector{Tuple{BitVector, Float64}}, v::Function)

Calulates the CDF for control variates with estimator W.
"""
function fcvW(y::Real, m::Real, samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    n = length(samples)
    sum_1 = 0
    fcv = 0
    mean = calculate_mean_cv(v, samples)
    for (σ, _) in samples
        sum_1 += (v(σ) - mean)^2
    end
    for (σ, devation) in samples
        Wi = 1 / n + (mean - v(σ)) * (mean - m) / sum_1
        fcv += Wi * indicator(devation, y)
    end
    return fcv
end

"""
    inverse_fcv(p, m, samples, v)

Calulates the p-quantile using the given p, theoratical mean value of control variate function,
the given sample and control variate function v.
"""
function inverse_fcv(p::Real, m::Real, samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    n = length(samples)
    sum_1 = 0
    sum_2 = 0
    y = 0
    mean = calculate_mean_cv(v, samples)
    for (σ, _) in samples
        sum_1 += (v(σ) - mean)^2
    end
    for (σ, devation) in samples
        Wi = 1 / n + (mean - v(σ)) * (mean - m) / sum_1
        sum_2 += Wi
        if sum_2 >= p
            y = devation
            break
        end
    end
    return y
end

"""
    var(samples,v)

Caluate the variance of the control variates in the given sample.
"""
function var(samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    n = length(samples)
    ret = 0
    sum_1 = 0
    mean = calculate_mean_cv(v, samples)
    for (σ, _) in samples
        sum_1 += (v(σ) - mean)^2
    end
    ret = 1 / (n - 1) * sum_1
    return ret
end

"""
    cov(quantile::Real,samples::Vector{Tuple{BitVector, Float64}},v::Function)

Calculate the covariance for quantile and control variates in the given sample.
"""
function cov(quantile::Real, samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    n = length(samples)
    ret = 0
    sum_1 = 0
    sum_2 = 0
    mean = calculate_mean_cv(v, samples)
    for (_, deviation) in samples
        sum_2 += indicator(deviation, quantile)
    end
    mean2 = sum_2 / n
    for (σ, deviation) in samples
        sum_1 += (v(σ) - mean) * (deviation - mean2)
    end
    ret = 1 / (n - 1) * sum_1
    return ret
end

"""
    ψ2_cv(quantile,p,samples,v)

Calculate the valure of ψ square of control variates. This value is used to calculate τ
Input the p-quantile, p, given samples and control variate function.
"""
function ψ2_cv(quantile::Real, p::Real, samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    ψ = p * (1 - p) - (cov(quantile, samples, v)^2) / var(samples, v)
    return ψ
end


"""
    η(δ, p, m, samples, v)

Calculate the value of η using a user-specified bandwidth δ, p, theoratical mean value of 
control variate function, the given sample and control variate function v. 
This value is used to calculate τ
"""
function η(δ::Real, p::Real, m::Real, samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    η = (inverse_fcv(p + δ, m, samples, v) - inverse_fcv(p - δ, m, samples, v)) / (2 * δ)
    return η
end

"""
    τ2(quantile,δ,p,m,samples,v)

Calculate the value of τ square with p-quantile, user-specified bandwidth δ, p, theoratical mean value of 
control variate function, the given sample and control variate function v. τ is used to calculate the confidence
interval.
"""
function τ2(quantile::Real, δ::Real, p::Real, m::Real, samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    τ = ψ2_cv(quantile, p, samples, v) * (η(δ, p, m, samples, v))^2
    return τ
end


"""
    z_α(α)

Calculate the inverse of CDF for Normal Distribution for a desired confidence level α.
This value is used to calculate the confidence interval.
"""
function z_α(α::Real)
    d = Normal(0.0, 1.0)
    z_α = quantile.(d, [1 - α / 2])
    return z_α[1]
end

"""
    confidence_interval(α, quantile, δ, p, m, samples, v)

Calculate the confidence interval given the desired confidence level α, the p-quantile, 
a user-specified bandwidth δ, p, theoratical mean value of control variate function, 
the given sample and control variate function v. Output a 3 element vector with the first
two being the confidence interval and the last element being the magnitude of the interval.
"""
function confidence_interval(α::Real, quantile::Real, δ::Real, p::Real, m::Real, 
        samples::Vector{Tuple{BitVector,Float64}}, v::Function)
    n = length(samples)
    diff = z_α(α) * sqrt(τ2(quantile, δ, p, m, samples, v)) / sqrt(n)
    i1 = quantile - diff
    i2 = quantile + diff
    return [i1, i2, 2 * diff]
end

function bin_list(n::Integer)
    if n == 1
        return [[1], [0]]
    else
        result = Vector{BitVector}(undef, 0)
        previous = bin_list(n - 1)
        for vector in previous
            vector2 = copy(vector)
            push!(vector, 1)
            push!(vector2, 0)
            push!(result, vector)
            push!(result, vector2)
        end
        return result
    end
end

end
