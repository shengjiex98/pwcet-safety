@info "Importing packages"
flush(stderr)

using Serialization
using Printf
using ControlSystemsBase
using ControlTimingSafety
using LinearAlgebra: I
using Distributions: Pareto, Normal, cdf, quantile

push!(LOAD_PATH, "../src")
using Experiments
using Benchmarks
using ContinuousSims: nominal_trajectory

@info "Setting parameters"
flush(stderr)

const JOB_ID = parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])
const TASK_ID = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])

# >>> Experiment parameters >>>
# BATCHSIZE = 100
const BATCHSIZE = 100_000

# Continuous SYStem definition
const SYS = benchmarks[:F1T]
const K = lqr(SYS, I, I)

# Time horizon
const H = 100 * 0.02

# Execution time distribution
const DIST = Pareto(1.5, 0.03)
# DIST = Normal(0.03, 0.005)
const P_MIN = quantile(DIST, 0.01)
const P_MAX = quantile(DIST, 0.99)

# Chosen quantiles
const I_VALUES = 0.01:0.01:1.0

# Utilization values
const U_VALUES = 0.01:0.01:1.0

const UTIL = TASK_ID/100
@assert UTIL ∈ U_VALUES "TASK_ID exceeds utilization range. Exiting."

PATH = "../data/nmc-grid/$JOB_ID-$DIST/"
# <<< Experiment parameters <<<

mkpath(PATH)

function get_period_uniform(percentage:: Real)
    @assert 0 <= util <= 1 "Percentage value has to be in [0, 1]"
    P_MIN + (P_MAX - P_MIN) * percentage
end

function get_q_from_period(period::Real; util::Real=1)
    @assert 0 < util <= 1 "Utilization value has to be in (0, 1]"
    cdf(DIST, period * util)
end

# Set initial conditions
x0 = fill(1., SYS.nx)
u0 = 0.
z0 = [x0; u0]

@info "Threads count:" Threads.nthreads()
@info "Distribution:" DIST
@info "System dynamics:" SYS K z0
@info "Parameters:" BATCHSIZE H
@info "Calculating the nominal_trajectory."
flush(stderr)

z_nom = nominal_trajectory(SYS, (x, t) -> -K * x, period, H, x0)

@info "Running simulations"
flush(stderr)
for i in I_VALUES
    period = get_period_uniform(i)
    q = get_q_from_period(period, util=UTIL)
    H_steps = floor(Int64, H / period)

    @info "Iteration parameters:" period q H_steps
    flush(stderr)

    # Construct automaton
    a = hold_kill(c2d(SYS, period), delay_lqr(SYS, period))

    filename = generate_filename(BATCHSIZE, q, period)
    if isfile("$PATH/$filename.jls")
        @info "$filename.jls exists, exiting."
        return
    end
    t = @elapsed data = generate_samples(a, z0, q, BATCHSIZE; H=H_steps, nominal_trajectory=z_nom)
    @info "Elapsed time:" t
    serialize("$PATH/$filename.jls", data)
    @info "Saved at $PATH/$filename.jls"
end
