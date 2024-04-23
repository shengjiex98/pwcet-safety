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

# ╔═╡ 89e10260-cbfe-49b8-932d-f5a1798a123a
pwd()

# ╔═╡ ebd19987-d113-4b8d-a0c7-06c5142310a4
ENV["LD_LIBRARY_PATH"] = ENV["LD_LIBRARY_PATH"] * ":$(pwd())/../../MATLABControlTest.jl/src"

# ╔═╡ 566cc4d1-e984-44e4-ab06-4d05d903335a
ref_signal = CSV.read("../../MATLABControlTest.jl/output-jumping1000-1e3-O1.csv", DataFrame)[!,:r]

# ╔═╡ 22f33099-3217-4d44-bf77-65873eae293f
function get_ref(t::Real)
	t_idx = floor(Int64, t / 0.1) + 1
	@boundscheck 1 ≤ t_idx ≤ 1000 || throw(ArgumentError("t=$t out of bound"))
	ref_signal[t_idx]
end

# ╔═╡ 8a2ca62b-f4cd-443e-af7b-ff9ac60caeca


# ╔═╡ Cell order:
# ╠═e37c5e34-0191-11ef-0c30-379698507b64
# ╠═89e10260-cbfe-49b8-932d-f5a1798a123a
# ╠═ebd19987-d113-4b8d-a0c7-06c5142310a4
# ╠═39bbf368-4d9b-447a-be29-a0b74f391cf3
# ╠═566cc4d1-e984-44e4-ab06-4d05d903335a
# ╠═22f33099-3217-4d44-bf77-65873eae293f
# ╠═8a2ca62b-f4cd-443e-af7b-ff9ac60caeca
