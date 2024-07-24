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
const BATCHSIZE = 30_000

# Continuous system definition
const SYS = ss(tf([3, 1],[1, 0.6, 1]))

# Time horizon
const H = 1000 * 0.1

# Chosen quantiles
const I_VALUES = 0.01:0.01:1.0

# Utilization values
const U_VALUES = 0.01:0.01:1.0

const COMPARE_MODE = ENV["COMPARE_MODE"]
@assert COMPARE_MODE ∈ ["y", "ref"] "COMPARE_MODE must be either 'y' or 'ref'"

const UTIL = TASK_ID/100
@assert UTIL ∈ U_VALUES "TASK_ID exceeds utilization range. Exiting."

# Timing file PATH
const TIMING_FILE = "$(@__DIR__)/../data/output-O2-samer/output-O2-samer-1.json"

# Save directory
const SAVE_PATH = "$(@__DIR__)/../data/mpc-grid/$JOB_ID/$TASK_ID/"

# <<< Experiment parameters <<<

# Reference and cycle reading
# const DATA = CSV.read("output-jumping1000-1e3-O1.csv", DataFrame)

# Normalizing factor: 100_000 cycles per 0.2 second
# const E_VALUES = sort(DATA[:, :t]) / 100_000 * 0.2

@info "Reading timing file"
flush(stderr)

const DATA = open(TIMING_FILE) do file
    JSON.parse(file)
end
const E_VALUES = sort(DATA["t"]) / 100_000 * 0.2

@info "Ref length" length(E_VALUES)

mkpath(SAVE_PATH)

function get_ref(t::Real)
    t_i = floor(Int64, t / 0.1) + 1
    @boundscheck 1 ≤ t_i ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
    DATA["r"][t_i]
end

# Get period as a quantile from the runtime distribution
function get_period_distribution(q::Real)
    e_i = ceil(Int64, q * length(E_VALUES))
    @boundscheck 1 ≤ e_i ≤ length(E_VALUES) || throw(ArgumentError("t=$t out of bound"))
    E_VALUES[e_i]
end

# Get period as a quantile uniformly between min and max runtimes
function get_period_uniform(percentage::Real)
    p_min = E_VALUES[1]
    p_max = E_VALUES[end]
    p_min + (p_max - p_min) * percentage
end

function get_q_from_period(period::Real; util::Real=1)
    @assert 0 < util <= 1 "Utilization value has to be in (0, 1]"
    e_i = findfirst(x -> x > period * util, E_VALUES)
    if e_i !== nothing; (e_i-1) / length(E_VALUES) else 1 end
end

function get_y(t::Real)
    t_i = floor(Int64, t / 0.1) + 1
    @boundscheck 1 ≤ t_i ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
    DATA["y"][t_i]
end

@info "Running simulations"
for i in I_VALUES
    # # Choosing the period from the distribution
    # period = get_period_quantile(i)
    # q = i

    # Choosing the period uniformly between min and max period
    period = get_period_uniform(i)
    q = get_q_from_period(period, util=UTIL)

    ref_values = map(t -> get_ref(t), 0:period:H)
    y_values = map(t -> get_y(t), 0:period:H)
    H_steps = length(ref_values)

    sysd = c2d(SYS, period)
    x0 = zeros(sysd.nx)

    @info "Threads count:" Threads.nthreads()
    @info "System dynamics:" SYS sysd
    @info "Parameters:" BATCHSIZE q period UTIL H H_steps
    @info "MPC compare mode:" COMPARE_MODE
    flush(stderr)

    filename = generate_filename(BATCHSIZE, q, period)
    if isfile("$SAVE_PATH/$filename.jls")
        @info "$filename.jls exists, exiting."
        # exit()
        continue
    end
    if COMPARE_MODE == "y"
        t = @elapsed data = generate_samples_mpc(sysd, x0, ref_values, q, BATCHSIZE, H=H_steps, compare=y_values)
    elseif COMPARE_MODE == "ref"
        t = @elapsed data = generate_samples_mpc(sysd, x0, ref_values, q, BATCHSIZE, H=H_steps, compare=ref_values)
    end
    @info "Elapsed time:" t
    serialize("$SAVE_PATH/$filename.jls", data)
    @info "Saved at $SAVE_PATH/$filename.jls"
end
