using Revise, MATLABControlTest, ControlSystemsBase, CSV, DataFrames

sys = ss(tf([3, 1],[1, 0.6, 1])); data = CSV.read("mpc_randomref_O1.csv", DataFrame)

period = 0.1; sysd = c2d(sys, period); x0 = zeros(sysd.nx); sim = run_simulation(sysd, x0, data[:,:r])

sim - data[:,:r]

@info maximum(abs.(sim - data[:,:r]))

@info maximum(abs.(data[:,:r]-data[:,:y]))

period = 0.1; sysd = c2d(sys, period); x0 = zeros(sysd.nx); sim = run_simulation(sysd, x0, data[:,:r])

function get_ref(t::Real)
    t_i = ceil(Int64, t / 2) + 1
    @boundscheck 1 ≤ t_i ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
    REF_VALUES[t_i]
end

# H = 0.1 * 1000
# period = 0.15; ref_values = get_ref(0:period:H);
