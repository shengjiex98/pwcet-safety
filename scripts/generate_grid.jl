@info "Importing packages"
flush(stderr)

using Serialization
using Printf
using ControlSystemsBase
using ControlTimingSafety
using LinearAlgebra: I
using Distributions: Pareto, Normal, cdf, quantile
using Printf

push!(LOAD_PATH, "$(@__DIR__)/../src")
using Experiments
using Benchmarks
using ContinuousSims: nominal_trajectory

@info "Setting parameters"
flush(stderr)

const JOB_ID = parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])
const TASK_ID = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])

# >>> Experiment parameters >>>
# BATCHSIZE = 100
const BATCHSIZE = 300_000

# Continuous SYStem definition
const SYS = benchmarks[:F1T]
const K = lqr(SYS, I, I)

# Time horizon
const H = 100 * 0.02

# Execution time distribution
const DIST = Pareto(1.5, 0.01)
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

PATH = "$(@__DIR__)/../data/nmc-grid/$JOB_ID-$DIST/$PERIOD/"
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

@info "Running simulations"
flush(stderr)
for u in U_VALUES
    q = get_q_from_period(PERIOD, util=u)
    if q == 0
        continue
    end

    @info "Iteration parameters:" PERIOD q H_STEPS
    flush(stderr)

    # filename = generate_filename(BATCHSIZE, q, PERIOD) * "-u$u"
    filename = @sprintf "b%.1e-q%.6f-h%.6g-u%.2f-th%i" BATCHSIZE q PERIOD u Threads.nthreads()
    if isfile("$PATH/$filename.jls")
        @info "$filename.jls exists, exiting."
        continue
    end
    t = @elapsed data = generate_samples(a, z0, q, BATCHSIZE; H=H_STEPS, nominal_trajectory=z_nom)
    @info "Elapsed time:" t
    serialize("$PATH/$filename.jls", data)
    @info "Saved at $PATH/$filename.jls"
end
