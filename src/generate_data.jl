"""
ARGS: MODE, b, q, h[, n]
MODE ∈ {"batch", "normal"}
b is BATCHSIZE
q is hit chance
n (optional) number of batches for batch mode
"""

using Serialization
using Printf
using ControlSystemsBase
using ControlTimingSafety
using LinearAlgebra: I
using Distributions: Pareto, cdf, quantile

push!(LOAD_PATH, "../lib")
using Experiments
using Benchmarks
using ContinuousSims: nominal_trajectory

# BATCHSIZE = 100
BATCHSIZE = 1_000_000
MODE = "normal"
# MODE, N = "batch", 100
H = 100 * 0.02
Q_VALUES = 0.01:0.01:0.99
DIST = Pareto(1.5, 0.01)
PATH = "../data/nmc-dist/"

TASK_ID = parse(Int64, ENV["SLURM_ARRAY_TASK_ID"])
if TASK_ID > length(Q_VALUES)
    println("TASK_ID exceeds available parameters. Exiting.")
    return
end

@info "Threads count:" Threads.nthreads()

q = Q_VALUES[TASK_ID]
period = quantile(DIST, q)
H_steps = floor(Int64, H / period)

# Set initial conditions
sys = benchmarks[:F1T]
K = lqr(sys, I, I)
x0 = fill(1., sys.nx)
u0 = 0.
z0 = [x0; u0]

z_nom = nominal_trajectory(sys, (x, t) -> -K * x, period, H, x0)

# Construct automaton
a = hold_kill(c2d(sys, period), delay_lqr(sys, period))

if MODE == "batch"
    nbatches = parse(Int64, ARGS[5])
    filename = generate_filename(BATCHSIZE, q, period, n)
    if isfile("$PATH/$filename.jls")
        @info "$filename.jls exists, exiting."
        return
    end
    @info "Parameters" MODE BATCHSIZE q period nbatches
    t = @elapsed batches = map(_ -> generate_samples(a, z0, q, BATCHSIZE; H=H_steps, nominal_trajectory=z_nom), 1:nbatches)
    @info t
    serialize("$PATH/$filename.jls", batches)
    # write("$PATH/$filename.txt", "$t")
elseif MODE == "normal"
    filename = generate_filename(BATCHSIZE, q, period)
    if isfile("$PATH/$filename.jls")
        @info "$filename.jls exists, exiting."
        return
    end
    @info "Parameters" MODE BATCHSIZE q period
    t = @elapsed data = generate_samples(a, z0, q, BATCHSIZE; H=H_steps, nominal_trajectory=z_nom)
    @info t
    serialize("$PATH/$filename.jls", data)
    # write("$PATH/$filename.txt", "$t")
else
    @error "First argument must be either 'batch' or 'normal'"
end

## === Generate a single batch for each parameter ===

# BATCHSIZEs = [100_000, 1_000_000, 10_000_000]
# qs = [0.9, 0.99, 0.999]

# for b in BATCHSIZEs, q in qs
#     path = "../data/batches"
#     filename = "b$(b/1_000_000)m-q$q"
#     t = time()
#     batches = generate_samples(a, z0, q, b; H=H)
#     t = time() - t
#     write("$path/$filename.txt", t)
#     serialize("$path/$filename.jls", batches)
# end


# # MODE ∈ {"batch", "normal"}
# MODE = ARGS[1]
# # Number of samples to generate in a batch
# BATCHSIZE = parse(Int64, ARGS[2])
# # Hit chance
# q = parse(Float64, ARGS[3])
# # Period
# period = parse(Float64, ARGS[4])
# H_steps = floor(Int64, H / period)
