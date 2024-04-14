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

JOB_ID = parse(Int64, ENV["SLURM_ARRAY_JOB_ID"])
TASK_ID = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])

# >>> Experiment parameters >>>
# BATCHSIZE = 100
BATCHSIZE = 1_000_000

DIST = Pareto(1.5, 0.001)
# DIST = Normal(0.03, 0.005)

# Time horizon
H = 100 * 0.02
# Values of quantiles
Q_VALUES = 0.01:0.01:0.99
SYS = :EWB
PATH = "../data/nmc-dist/$(String(SYS))/$DIST/"
# <<< Experiment parameters <<<

mkpath(PATH)

if TASK_ID > length(Q_VALUES)
    println("TASK_ID exceeds available parameters. Exiting.")
    exit()
end

q = Q_VALUES[TASK_ID]
period = quantile(DIST, q)
H_steps = floor(Int64, H / period)

# Set initial conditions
sys = benchmarks[SYS]
K = lqr(sys, I, I)
x0 = fill(1., sys.nx)
u0 = 0.
z0 = [x0; u0]

@info "Threads count:" Threads.nthreads()
@info "Distribution and quantiles:" DIST Q_VALUES
@info "System dynamics:" sys K z0
@info "Parameters:" BATCHSIZE H q period

z_nom = nominal_trajectory(sys, (x, t) -> -K * x, period, H, x0)

# Construct automaton
a = hold_kill(c2d(sys, period), delay_lqr(sys, period))

filename = generate_filename(BATCHSIZE, q, period)
if isfile("$PATH/$filename.jls")
    @info "$filename.jls exists, exiting."
    return
end
t = @elapsed data = generate_samples(a, z0, q, BATCHSIZE; H=H_steps, nominal_trajectory=z_nom)
@info t
serialize("$PATH/$filename.jls", data)
@info "Saved at $PATH$filename.jls"
