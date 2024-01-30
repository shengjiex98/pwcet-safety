using Serialization

push!(LOAD_PATH, "../lib")
using Experiments
using Benchmarks
using ControlVariates

using ControlSystemsBase
using ControlTimingSafety

@info Threads.nthreads()

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

# Generate multiple batches
# batchsize = 10_000
# nbatches = 100
# filename = "100x10k.jls"
# @time batches = map(_ -> generate_samples(a, z0, q, batchsize; H=H), 1:nbatches)
# serialize("../data/batches/$filename", batches)

# Generate a single batch
batchsize = 10_000_000
filename = "10m_q=$q.jls"
@time batches = generate_samples(a, z0, q, batchsize; H=H)
serialize("../data/batches/$filename", batches)
