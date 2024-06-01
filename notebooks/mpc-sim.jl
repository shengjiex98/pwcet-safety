### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ e37c5e34-0191-11ef-0c30-379698507b64
begin
	import Pkg
	Pkg.activate("..")

	using Printf
	using Serialization
	using DelimitedFiles
	using Distributions: Distribution, Normal, Pareto, Uniform, cdf, pdf, quantile
	using Plots
	plotlyjs()

	push!(LOAD_PATH, "../src")
	using Experiments

	push!(LOAD_PATH, "../../")
	
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 39bbf368-4d9b-447a-be29-a0b74f391cf3
using Revise, MATLABControlTest, ControlSystemsBase, CSV, DataFrames

# ╔═╡ 566cc4d1-e984-44e4-ab06-4d05d903335a
begin
	function get_ref(t::Real)
		t_i = floor(Int64, t / 0.1) + 1
		@boundscheck 1 ≤ t_i ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
		DATA[t_i, :r]
	end
	function get_period(q::Real)
    	e_i = ceil(Int64, q * length(E_VALUES))
		@boundscheck 1 ≤ e_i ≤ length(E_VALUES) || throw(ArgumentError("t=$t out of 		bound"))
		E_VALUES[e_i]
	end
	JOB_ID = 39927918
	DATA = CSV.read("../MATLABControlTest.jl/output-jumping1000-1e3-O1.csv", DataFrame)
	MPC_data = CSV.read("../data-proxy/mpc-$JOB_ID.csv", DataFrame)
	E_VALUES = sort(DATA[:, :t])
	q_values = 0.01:0.01:0.99
	H = 1000 * 0.1
	period = []
	for q in q_values
		push!(period, get_period(q)) 
	end
end

# ╔═╡ cca95454-a0a8-4b99-b620-9b14af48c701
let hs = sort(MPC_data[:, 2])
	n = length(hs)
	cdf_values = (1:n) ./ n
	fig2 = plot(hs, cdf_values, xlabel="period", ylabel="cdf", title="CDF for mpc")
end

# ╔═╡ 07e901e0-5a41-4b02-8bd7-c50e67022b70


# ╔═╡ Cell order:
# ╠═e37c5e34-0191-11ef-0c30-379698507b64
# ╠═39bbf368-4d9b-447a-be29-a0b74f391cf3
# ╠═566cc4d1-e984-44e4-ab06-4d05d903335a
# ╠═cca95454-a0a8-4b99-b620-9b14af48c701
# ╠═07e901e0-5a41-4b02-8bd7-c50e67022b70
