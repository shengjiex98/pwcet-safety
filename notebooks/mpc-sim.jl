### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

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
	JOB_ID = 41671552
	file_num = 1
	MPC_data = CSV.read("../data-proxy/mpc-flags/ref/mpc-$JOB_ID-$file_num.csv", DataFrame; header=false)
end

# ╔═╡ 07e901e0-5a41-4b02-8bd7-c50e67022b70
function plot_results(
		q_values::Vector{<:Real},
		h_values::Vector{<:Real},
		quantiles::Vector{<:Real};
		filter_fn::Union{Function, Nothing}=nothing,
		confidence::Union{Matrix{<:Real}, Nothing}=nothing,
		color::Union{Symbol, Vector{Symbol}}=:lightblue,
		title::String="",
		draw_surface::Bool=true,
		draw_all_points::Bool=true,
		mode::String="free",
		az::Real=55,
		el::Real=15)

	# Filter data if needed; `selection` is a BitVector representing points to show
	if !isnothing(filter_fn)
		selection = filter_fn.(q_values, h_values, quantiles)
		q_values = q_values[selection]
		h_values = h_values[selection]
		quantiles = quantiles[selection]
		if color isa Vector
			color = color[selection]
		end
		if confidence isa Matrix
			confidence = confidence[selection,:]
		end
	end
	
	# Set common properties
	plt = plot(title=title, legend=:topleft, proj_type=:ortho)

	hovertext=map(q_values, h_values, quantiles) do q, h, dev
		@sprintf "q=%.3f h=%.3f dev=%.3f" q h dev
	end
	
	# Show different graph depending on `mode`
	if mode == "q"
		scatter!(q_values, quantiles,
			# xlim=(qmin, qmax), 
			# ylim=(0, maximum(quantiles)),
			xlabel="quantile", ylabel="deviation", color=color,
			hover=hovertext
		)
		if !isnothing(confidence)
			plot!(q_values, quantiles,
				ribbon=(confidence[:,1], confidence[:,2]))
		end
	elseif mode == "period"
		scatter!(h_values, quantiles,
			# xlim=(hmin, hmax), 
			# ylim=(0, maximum(quantiles)),
			xlabel="period", ylabel="deviation", color=color)
		if !isnothing(confidence)
			plot!(h_values, quantiles,
				ribbon=(confidence[:,1], confidence[:,2]))
		end
	elseif mode == "free"
	    if draw_surface
	        surface!(q_values, h_values, quantiles)
	    end
		scatter!(q_values, h_values, quantiles,
			markersize=2.5,
			# xlim=(qmin, qmax), ylim=(hmin, hmax), 
			# zlim=(0, maximum(quantiles)),
			xlabel="quantile", ylabel="period", zlabel="deviation",
			color=color)
			plot!(camera=(az, el))
	else
		throw(ArgumentError("`mode` has to be `free`, `q`, or `period`. $mode is given"))
	end
	plt
end

# ╔═╡ b18059b1-4eb3-486c-ad09-11cf386fff20
begin
	b = 1_000_000
	p = 0.99
	hs = sort(MPC_data[:, 2])
	n = length(hs)
	qs = (1:n) ./ (n+1)
end

# ╔═╡ 07df95d3-8763-47bd-a6a4-3b5a42e3ccb4
let go
	md"""
	| | |
	|--:|:--|
	| camera azimuth   | $(@bind az Slider(-180:5:180, default=55, show_value=true)) |
	| camera elevation | $(@bind el Slider(0:5:90, default=15, show_value=true)) |
	| mode             | $(@bind mode Radio(["free", "q", "period"], default="free"))
	| min q value      | $(@bind qmin Slider(qs, default=qs[1], show_value=true)) |
	| max q value      | $(@bind qmax Slider(qs, default=qs[end], show_value=true)) |
	| min h value      | $(@bind hmin Slider(hs, default=hs[1], show_value=true)) |
	| max h value      | $(@bind hmax Slider(hs, default=hs[end], show_value=true)) |
	| surface          | $(@bind surf CheckBox(default=false)) |
	| cap 			   | $(@bind cap  Slider(1500:1:18000, default=18000, show_value=true)) |
	| all points       | $(@bind all_points CheckBox(default=true)) |
	"""
end

# ╔═╡ b75abdc8-f875-48b5-9497-cf67b6c305fe
let
	JOB_ID = 41671552
	file_num = 20
	full_matrix_ref = readdlm("../data-proxy/mpc-flags/ref/mpc-$JOB_ID-$file_num.csv", ',')
	colors = let
		res = fill(:lightblue, size(full_matrix_ref, 1))
		mindev = cap
		for (i, row) in enumerate(eachrow(full_matrix_ref))
			if row[4] < mindev
				mindev = row[4]
				res[i] = :red
			end
		end
		res
	end
	filter_fn = (q, h, dev) -> 
		qmin <= q <= qmax && 
		hmin <= h <= hmax &&
		dev <= cap
	fig1 = plot_results(full_matrix_ref[:,1], full_matrix_ref[:,2], full_matrix_ref[:,4],
		confidence=[full_matrix_ref[:,4]-full_matrix_ref[:,3] full_matrix_ref[:,5]-full_matrix_ref[:,4]],
		filter_fn=filter_fn, title="ref", draw_surface=surf,
		mode=mode, az=az, el=el, color=colors)
	
	JOB_ID_2 = 41703162
	full_matrix_y = readdlm("../data-proxy/mpc-flags/y/mpc-$JOB_ID_2-$file_num.csv", ',')
	colors = let
		res = fill(:lightblue, size(full_matrix_y, 1))
		mindev = cap
		for (i, row) in enumerate(eachrow(full_matrix_y))
			if row[4] < mindev
				mindev = row[4]
				res[i] = :red
			end
		end
		res
	end
	filter_fn = (q, h, dev) -> 
		qmin <= q <= qmax && 
		hmin <= h <= hmax &&
		dev <= cap
	fig2 = plot_results(full_matrix_y[:,1], full_matrix_y[:,2], full_matrix_y[:,4],
		confidence=[full_matrix_y[:,4]-full_matrix_y[:,3] full_matrix_y[:,5]-full_matrix_y[:,4]],
		filter_fn=filter_fn, title="y", draw_surface=surf,
		mode=mode, az=az, el=el, color=colors)
	plot(fig1, fig2, plot_title="MPC")
end

# ╔═╡ Cell order:
# ╠═e37c5e34-0191-11ef-0c30-379698507b64
# ╠═39bbf368-4d9b-447a-be29-a0b74f391cf3
# ╠═566cc4d1-e984-44e4-ab06-4d05d903335a
# ╟─07e901e0-5a41-4b02-8bd7-c50e67022b70
# ╠═b18059b1-4eb3-486c-ad09-11cf386fff20
# ╟─07df95d3-8763-47bd-a6a4-3b5a42e3ccb4
# ╠═b75abdc8-f875-48b5-9497-cf67b6c305fe
