@info "Importing packages"
flush(stderr)

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

const JOB_ID = parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])
const TASK_ID = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])

# >>> Experiment parameters >>>
# BATCHSIZE = 100
const BATCHSIZE = 300_000

# Continuous SYStem definition
const SYSNAME = :F1T
const SYS = benchmarks[SYSNAME]
const K = lqr(SYS, I, I)

# Time horizon
const H = 100 * 0.02

# Timing file PATH
@info "Reading timing file"
flush(stderr)
const TIMING_FILE = "$(@__DIR__)/slre.json"
const DATA = open(TIMING_FILE) do file
    JSON.parse(file)
end
const E_VALUES = sort(DATA["t"]) / mean(DATA["t"]) * 0.02

const DIST = E_VALUES

# Execution time distribution
# const DIST = Pareto(1.5, 0.01)
# DIST = Normal(0.03, 0.005)
const P_MIN = quantile(DIST, 0.01)
const P_MAX = quantile(DIST, 0.99)

# Chosen quantiles
const I_VALUES = 0.01:0.01:1.0

# Utilization values
const U_VALUES = 0.01:0.01:1.0

const IDX = TASK_ID/100
@assert IDX ∈ I_VALUES "TASK_ID exceeds index range.."
const PERIOD = P_MIN + (P_MAX - P_MIN) * IDX

const PATH = "$(@__DIR__)/../data-csv/$JOB_ID/"
# Create PATH in case it does not exist yet
mkpath(PATH)

# <<< Experiment parameters <<<

function get_q_from_period(period::Real; util::Real=1)
    @assert 0 < util <= 1 "Utilization value has to be in (0, 1]"
    cdf(DIST, period * util)
end

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
flush(stderr)
for u in U_VALUES
    q = get_q_from_period(PERIOD, util=u)
    if q == 0
        continue
    end

    @info "Iteration parameters:" PERIOD q H_STEPS
    flush(stderr)

    t = @elapsed data = generate_samples(a, z0, q, BATCHSIZE; H=H_STEPS, nominal_trajectory=z_nom)
    @info "Elapsed time:" t
    p99, p99_lower, p99_upper = summarize_data(data, p=0.99, α=0.05)
    push!(df, (String(SYSNAME), "slre", BATCHSIZE, q, PERIOD, u, p99, p99_lower, p99_upper))
end

# Save df to CSV file
CSV.write("$PATH/$TASK_ID.csv", df)
