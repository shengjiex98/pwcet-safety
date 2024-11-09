@info "Importing packages"
flush(stderr)

using CSV
using JSON
using Serialization
using Printf
using ControlSystemsBase
using ControlTimingSafety
using LinearAlgebra: I
import Distributions: Pareto, Normal, cdf, quantile
using Statistics
using Printf
using DataFrames

push!(LOAD_PATH, "$(@__DIR__)/../src")
using Experiments
using Benchmarks
using ContinuousSims: nominal_trajectory

function cdf(dist::Vector{<:Real}, value::Real)
    return (dist .<= value) .|> Int64 |> mean
end

@info "Setting parameters"
flush(stderr)

const JOB_ID = 1
const TASK_ID = 30

# >>> Experiment parameters >>>
# BATCHSIZE = 100
const BATCHSIZE = 1000

# Continuous SYStem definition
const SYSNAME = :F1T
const SYS = benchmarks[SYSNAME]
const K = lqr(SYS, I, I)

const H = 100 * 0.02
const PERIOD = 0.0204241

const PATH = "$(@__DIR__)/../data-csv/$JOB_ID/"
# Create PATH in case it does not exist yet
mkpath(PATH)

const H_STEPS = floor(Int64, H / PERIOD)

# Set initial conditions
x0 = fill(1., SYS.nx)
u0 = 0.
z0 = [x0; u0]

@info "Threads count:" Threads.nthreads()
@info "Distribution:" DIST
@info "System dynamics:" SYS K z0
@info "Parameters:" BATCHSIZE H PERIOD
@info "Calculating the nominal_trajectory."
flush(stderr)

z_nom = nominal_trajectory(SYS, (x, t) -> -K * x, PERIOD, H, x0)

# Construct automaton
a = hold_kill(c2d(SYS, PERIOD), delay_lqr(SYS, PERIOD))

df = DataFrame(
    system=String[],
    distribution=String[],
    batchsize=Int64[],
    hit_chance=Float64[],
    period=Float64[],
    utilization=Float64[],
    p99=Float64[],
    p99_lower=Float64[],
    p99_upper=Float64[]
)
@info "Running simulations"
q = 0.1

@info "Iteration parameters:" PERIOD q H_STEPS
flush(stderr)

t = @elapsed data = generate_samples(a, z0, q, BATCHSIZE; H=H_STEPS, nominal_trajectory=z_nom)
@info data
@info "Elapsed time:" t
p99, p99_lower, p99_upper = summarize_data(data, p=0.99, α=0.05)
@info "Quantiles" p99 p99_lower p99_upper
