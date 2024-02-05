using Serialization

push!(LOAD_PATH, "../lib")
using Experiments
using Benchmarks
using ControlVariates

using ControlSystemsBase
using ControlTimingSafety

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
q = 0.99
# Sequence length
H = 100

## === Generate multiple batches ===

# batchsize = 10_000
# nbatches = 100
# filename = "100x10k.jls"
# @time batches = map(_ -> generate_samples(a, z0, q, batchsize; H=H), 1:nbatches)
# serialize("../data/batches/$filename", batches)

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


# === Generate a single batch from command line arguments ===

b = parse(Int64, ARGS[1])
q = parse(Float64, ARGS[2])

path = "../data/batches"
filename = "b$(b/1_000_000)m-q$q-th$(Threads.nthreads())"
t = time()
@info b q
t = @elapsed batches = generate_samples(a, z0, q, b; H=H)
@info t
write("$path/$filename.txt", "$t")
serialize("$path/$filename.jls", batches)
