"""
ARGS: MODE, b, q, h[, n]
MODE ∈ {"batch", "normal"}
b is batchsize
q is hit chance
n (optional) number of batches for batch mode
"""

using Serialization
using Printf
using ControlSystemsBase
using ControlTimingSafety
using LinearAlgebra: I

push!(LOAD_PATH, "../lib")
using Experiments
using Benchmarks
using ContinuousSims: nominal_trajectory

println("Threads count: $(Threads.nthreads())")

# Sequence length
H = 100 * 0.02

# MODE ∈ {"batch", "normal"}
MODE = ARGS[1]
# Number of samples to generate in a batch
batchsize = parse(Int64, ARGS[2])
# Hit chance
q = parse(Float64, ARGS[3])
# Period
period = parse(Float64, ARGS[4])
H_steps = ceil(Int64, H / period)

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
    path = "../data/batches-samenom"
    nbatches = parse(Int64, ARGS[5])
    filename = generate_filename(batchsize, q, period, n)
    @info "Parameters" batchsize q nbatches
    t = @elapsed batches = map(_ -> generate_samples(a, z0, q, batchsize; H=H_steps, nominal_trajectory=z_nom), 1:nbatches)
    @info t
    serialize("$path/$filename.jls", batches)
    write("$path/$filename.txt", "$t")
elseif MODE == "normal"
    path = "../data/nmc-samenom"
    filename = generate_filename(batchsize, q, period)
    @info "Parameters" batchsize q
    t = @elapsed data = generate_samples(a, z0, q, batchsize; H=H_steps, nominal_trajectory=z_nom)
    @info t
    serialize("$path/$filename.jls", data)
    write("$path/$filename.txt", "$t")
else
    @error "First argument must be either 'batch' or 'normal'"
end

## === Generate a single batch for each parameter ===

# batchsizes = [100_000, 1_000_000, 10_000_000]
# qs = [0.9, 0.99, 0.999]

# for b in batchsizes, q in qs
#     path = "../data/batches"
#     filename = "b$(b/1_000_000)m-q$q"
#     t = time()
#     batches = generate_samples(a, z0, q, b; H=H)
#     t = time() - t
#     write("$path/$filename.txt", t)
#     serialize("$path/$filename.jls", batches)
# end
