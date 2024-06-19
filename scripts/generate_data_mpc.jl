@info "Importing packages"
flush(stderr)

using JSON
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
const BATCHSIZE = 30_000

# Continuous system definition
const SYS = ss(tf([3, 1],[1, 0.6, 1]))

# Number of json file
const FILE_NUM = 1:30

# Reference and cycle reading
# const DATA = CSV.read("output-jumping1000-1e3-O1.csv", DataFrame)

# Normalizing factor: 100_000 cycles per 0.2 second
# const E_VALUES = sort(DATA[:, :t]) / 100_000 * 0.2

DATA = [Dict() for _ in FILE_NUM]
E_VALUES =[Float64[] for _ in FILE_NUM]
for i in FILE_NUM
    file_path = "./O2/output-O2-$i.json"
    DATA[i] = open(file_path) do file
        JSON.parse(file)
    end
    E_VALUES[i] = sort(DATA[i]["t"]) / 100_000 * 0.2
end

# Time horizon
const H = 1000 * 0.1

# Chosen quantiles
const Q_VALUES = 0.01:0.01:0.99

# Save directory
const PATH = "../data/mpc/$JOB_ID/"
# <<< Experiment parameters <<<

for i in FILE_NUM
    dir_path = joinpath(PATH, string(i))
    mkpath(dir_path)
end

if TASK_ID > length(Q_VALUES)
    println("TASK_ID exceeds available parameters. Exiting.")
    exit()
end

function get_ref(t::Real, i::Integer)
    t_i = floor(Int64, t / 0.1) + 1
    @boundscheck 1 ≤ t_i ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
    DATA[i]["r"][t_i]
end

function get_period(q::Real, i::Integer)
    e_i = ceil(Int64, q * length(E_VALUES[i]))
    @boundscheck 1 ≤ e_i ≤ length(E_VALUES[i]) || throw(ArgumentError("t=$t out of 		bound"))
    E_VALUES[i][e_i]
end

function get_y(t::Real, i::Integer)
    t_i = floor(Int64, t / 0.1) + 1
    @boundscheck 1 ≤ t_i ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
    DATA[i]["y"][t_i]
end

for i in FILE_NUM
    q = Q_VALUES[TASK_ID]
    period = get_period(q, i)
    ref_values = map(t -> get_ref(t, i), 0:period:H)
    ref_values_y = map(t -> get_y(t, i), 0:period:H)
    H_steps = length(ref_values)

    sysd = c2d(SYS, period)
    x0 = zeros(sysd.nx)

    @info "Threads count:" Threads.nthreads()
    @info "System dynamics:" SYS sysd
    @info "Parameters:" BATCHSIZE H q period
    flush(stderr)

    filename = generate_filename(BATCHSIZE, q, period)
    if isfile("$PATH/$i/$filename.jls")
        @info "$filename.jls exists, exiting."
        exit()
    end
    t = @elapsed data = generate_samples_mpc_with_multi_ref(sysd, x0, ref_values, ref_values_y, q, BATCHSIZE, H=H_steps)
    @info t
    serialize("$PATH/$i/$filename.jls", data)
    @info "Saved at $PATH/$i/$filename.jls"
end
