@info "Importing packages"
flush(stderr)

using Serialization
using Printf
using CSV
using DataFrames
using ControlSystemsBase
using ControlTimingSafety
using LinearAlgebra: I
using Distributions: Pareto, Normal, cdf, quantile

using MATLABControlTest

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

# Continuous system definition
const SYS = ss(tf([3, 1],[1, 0.6, 1]))

# Reference and cycle reading
const DATA = CSV.read("output-jumping1000-1e3-O1.csv", DataFrame)

# Normalizing factor: 100_000 cycles per 0.2 second
const E_VALUES = sort(DATA[:, :t]) / 100_000 * 0.2

# Time horizon
const H = 1000 * 0.1

# Chosen quantiles
const Q_VALUES = 0.01:0.01:0.99

# Save directory
const PATH = "../data/mpc/$JOB_ID/"
# <<< Experiment parameters <<<

mkpath(PATH)

if TASK_ID > length(Q_VALUES)
    println("TASK_ID exceeds available parameters. Exiting.")
    exit()
end

function get_ref(t::Real)
	t_i = floor(Int64, t / 0.1) + 1
	@boundscheck 1 ≤ t_i ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
	DATA[t_i, :r]
end

function get_period(q::Real)
    e_i = ceil(Int64, q * length(E_VALUES))
    @boundscheck 1 ≤ e_i ≤ length(E_VALUES) || throw(ArgumentError("t=$t out of bound"))
    E_VALUES[e_i]
end

q = Q_VALUES[TASK_ID]
period = get_period(q)
H_steps = floor(Int64, H / period)
ref_values = map(get_ref, 0:period:H)

sysd = c2d(SYS, period)
x0 = zeros(sysd.nx)

@info "Threads count:" Threads.nthreads()
@info "System dynamics:" SYS sysd
@info "Parameters:" BATCHSIZE H q period
flush(stderr)

filename = generate_filename(BATCHSIZE, q, period)
if isfile("$PATH/$filename.jls")
    @info "$filename.jls exists, exiting."
    exit()
end
t = @elapsed data = generate_samples_mpc(sysd, x0, ref_values, q, BATCHSIZE, H=H_steps)
@info t
serialize("$PATH/$filename.jls", data)
@info "Saved at $PATH/$filename.jls"
