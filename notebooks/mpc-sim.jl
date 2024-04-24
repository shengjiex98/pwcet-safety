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
	using StatsBase, MATLABControlTest, ControlSystemsBase, CSV, DataFrames
	
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 89e10260-cbfe-49b8-932d-f5a1798a123a
pwd()

# ╔═╡ ebd19987-d113-4b8d-a0c7-06c5142310a4
if "LD_LIBRARY_PATH" in keys(ENV)
	ENV["LD_LIBRARY_PATH"] = ENV["LD_LIBRARY_PATH"] * ":$(pwd())/../../MATLABControlTest.jl/src"
else
	ENV["LD_LIBRARY_PATH"] = "$(pwd())/../../MATLABControlTest.jl/src"
end

# ╔═╡ e234e781-9e83-4b89-8700-f63951b65ab6
df = CSV.read("../../MATLABControlTest.jl/output-jumping1000-1e3-O1.csv", DataFrame)

# ╔═╡ c2b10163-1380-4466-aa50-07d0ffe10f7f
ref_signal = df[!,:r]

# ╔═╡ 22f33099-3217-4d44-bf77-65873eae293f
function get_ref(t::Real)
	t_idx = floor(Int64, t / 0.1) + 1
	@boundscheck 1 ≤ t_idx ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
	ref_signal[t_idx]
end

# ╔═╡ 8a2ca62b-f4cd-443e-af7b-ff9ac60caeca
scatter(df[!,:t], xlim=[0, 100])

# ╔═╡ d43f9865-1358-4544-9bc3-97006a340ef9
plot(sort(df[!,:t]), (1:nrow(df))/nrow(df))

# ╔═╡ 8a324244-2037-4dae-bd7d-8d5d6dd3359a
let
	scatter(df[!,:r], xlim=[0, 100])
	plot!(df[!,:y])
end

# ╔═╡ ed5d3ef6-1b05-4a17-b386-0b493e9d848e
sample(df[!,:t])

# ╔═╡ Cell order:
# ╠═e37c5e34-0191-11ef-0c30-379698507b64
# ╠═89e10260-cbfe-49b8-932d-f5a1798a123a
# ╠═ebd19987-d113-4b8d-a0c7-06c5142310a4
# ╠═e234e781-9e83-4b89-8700-f63951b65ab6
# ╠═c2b10163-1380-4466-aa50-07d0ffe10f7f
# ╠═22f33099-3217-4d44-bf77-65873eae293f
# ╠═8a2ca62b-f4cd-443e-af7b-ff9ac60caeca
# ╠═d43f9865-1358-4544-9bc3-97006a340ef9
# ╠═8a324244-2037-4dae-bd7d-8d5d6dd3359a
# ╠═ed5d3ef6-1b05-4a17-b386-0b493e9d848e
