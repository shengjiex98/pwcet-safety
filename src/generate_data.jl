"""
ARGS: MODE, b, q[, n]
MODE ∈ {"batch", "normal"}
b is batchsize
q is hit chance
n (optional) number of batches for batch mode
"""

using Serialization
using ControlSystemsBase
using ControlTimingSafety

push!(LOAD_PATH, "../lib")
using Experiments
using Benchmarks
using ControlVariates

println("Threads count: $(Threads.nthreads())")

# Construct automaton
period = 0.02
sys = benchmarks[:F1T]
a = hold_kill(c2d(sys, period), delay_lqr(sys, period))

# Set initial conditions
x0 = 1.
u0 = 0.
z0 = [fill(x0, size(sys.A, 1)); u0]

# Hit chance
# q = 0.999
# Sequence length
H = 100

MODE = ARGS[1]
batchsize = parse(Int64, ARGS[2])
q = parse(Float64, ARGS[3])

if MODE == "batch"
    path = "../data/batches"
    nbatches = parse(Int64, ARGS[4])
    filename = "b$(batchsize/1_000)k-q$q-n$nbatches-th$(Threads.nthreads())"
    @info "Parameters" batchsize q nbatches
    t = @elapsed batches = map(_ -> generate_samples(a, z0, q, batchsize; H=H), 1:nbatches)
    @info t
    serialize("$path/$filename.jls", batches)
    write("$path/$filename.txt", "$t")
elseif MODE == "normal"
    path = "../data/nmc"
    filename = "b$(batchsize/1_000_000)m-q$q-th$(Threads.nthreads())"
    @info "Parameters" batchsize q
    t = @elapsed data = generate_samples(a, z0, q, batchsize; H=H)
    @info t
    serialize("$path/$filename.jls", data)
    write("$path/$filename.txt", "$t")
else
    @error "First argument must be either \"batch\" or \"normal\""
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
