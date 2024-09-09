### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ b216a1aa-dc98-11ee-0312-21d71fee5020
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
	
	using PlutoUI
	TableOfContents()
end

# ╔═╡ 0becd1ea-ade9-43d4-9039-fb335de62c6d
md"""
# Visualization for Experiment Results
"""

# ╔═╡ 305d3f52-0208-4fe7-b362-1f9b9424b705
md"""
## Load Packages
"""

# ╔═╡ 9bd08151-ce7f-4d57-9dd7-7af79ff42d3f
md"""
## Defining Plotting Function
"""

# ╔═╡ ce5e566b-b1e3-4e79-9156-4237a170546e
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
	plt = plot(title=title, proj_type=:ortho)

	hovertext=map(q_values, h_values, quantiles) do q, h, dev
		@sprintf "q=%.3f h=%.3f dev=%.3f" q h dev
	end
	
	# Show different graph depending on `mode`
	if mode == "q"
		scatter!(q_values, quantiles,
			# xlim=(qmin, qmax), 
			# ylim=(0, maximum(quantiles)),
			xlabel="Quantile", ylabel="Deviation", color=color,
			hover=hovertext,
			label=nothing
		)
		if !isnothing(confidence)
			plot!(q_values, quantiles,
				ribbon=(confidence[:,1], confidence[:,2]),
				label=nothing)
		end
	elseif mode == "period"
		scatter!(h_values, quantiles,
			# xlim=(hmin, hmax), 
			# ylim=(0, maximum(quantiles)),
			xlabel="Period (s)", ylabel="Deviation", color=color,
			label=nothing
		)
		if !isnothing(confidence)
			plot!(h_values, quantiles,
				ribbon=(confidence[:,1], confidence[:,2]),
				label="95% Confidence Interval")
		end
	elseif mode == "free"
	    if draw_surface
	        surface!(q_values, h_values, quantiles)
	    end
		scatter!(q_values, h_values, quantiles,
			markersize=2.5,
			# xlim=(qmin, qmax), ylim=(hmin, hmax), 
			# zlim=(0, maximum(quantiles)),
			xlabel="Quantile", ylabel="Period (s)", zlabel="Deviation",
			color=color)
			plot!(camera=(az, el))
	else
		throw(ArgumentError("`mode` has to be `free`, `q`, or `period`. $mode is given"))
	end
	plt
end

# ╔═╡ 4cdad6fd-9976-4bbf-b6de-95cb6a973fc7
md"""
## Plot Data

We select a distribution and plot its PDF/CDF below:
"""

# ╔═╡ aca1dd0e-89b9-4a66-ae0b-2af07cea3636
# dist = Normal(0, 0)
# dist = Normal(0.03, 0.005)
dist = Pareto(1.5, 0.03)
# dist = Uniform(0, 0.06)

# ╔═╡ 9c608fab-8fc9-4bf9-be68-96468257e6a8
let hs = range(quantile.(dist, [0.01, 0.99])..., 100)
	fig1 = plot(hs, pdf.(dist, hs), xlabel="T", ylabel="pdf", label="")
	fig2 = plot(hs, cdf.(dist, hs), xlabel="T", ylabel="cdf", label="")
	plot(fig1, fig2, plot_title="PDF/CDF for $dist")
end

# ╔═╡ 7b5d607e-da02-4060-8a9c-60c541fe34aa
begin
	b = 1_000_000
	p = 0.99
	qs = 0.01:0.01:0.99
	hs = quantile.(dist, qs)
end

# ╔═╡ d4549398-3ca3-44ba-b016-ca55bf4056cf
md"""
$(@bind go Button("Reset"))
"""

# ╔═╡ dde0af91-d0ff-4e4d-a4d4-f8fb770987dd
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
	| cap 			   | $(@bind cap  Slider(1:0.5:15, default=2, show_value=true)) |
	| all points       | $(@bind all_points CheckBox(default=true)) |
	"""
end

# ╔═╡ d12e5d80-7099-4564-8144-778c812b106c
plt = let
	SYS = "F1T"
	# Hack to remove the prefix "Distributions." from distribution name
	DISTNAME = join(split("$dist", ".")[2:end], ".")
	full_matrix = readdlm("../data-proxy/nmc-dist-$SYS-$DISTNAME.csv", ',')
	@info full_matrix
	colors = let
		res = fill(:lightblue, size(full_matrix, 1))
		mindev = cap
		for (i, row) in enumerate(eachrow(full_matrix))
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
	plot_results(full_matrix[:,1], full_matrix[:,2], full_matrix[:,4],
		confidence=[full_matrix[:,4]-full_matrix[:,3] full_matrix[:,5]-full_matrix[:,4]],
		filter_fn=filter_fn, 
		# title="Deviation vs. Period for the F1/10 model", 
		draw_surface=surf,
		mode=mode, az=az, el=el)
		# mode=mode, az=az, el=el, color=colors)
end

# ╔═╡ 5a426b8d-e062-43b2-81ef-77f8245292ee
savefig(plt, "fig1.pdf")

# ╔═╡ ae494460-6f76-4562-a6b4-42f0bc6df262
let
	dist = Pareto(1.5, 0.01)
	
	points_small, quantiles_small = deserialize("../data-proxy/dev-quantiles.jls")
	points_large, quantiles_large = deserialize("../data-proxy/dev-quantiles-zoom.jls")
	points = [points_small points_large]
	quantiles = vcat(quantiles_small, quantiles_large)
	
	pts, qts = points, quantiles
	filter_fn = (q, h, dev) -> 
		q < cdf(dist, h) &&
		# qmin <= q <= qmax && 
		# hmin <= h <= hmax &&
		dev <= cap
	colors = let
		res = fill(:lightblue, length(qts))
		mindev = cap
		for (i, row) in enumerate(zip(pts[1,:], pts[2,:], qts))
			if row[3] < mindev && row[1] < cdf(dist, row[2])
				mindev = row[3]
				res[i] = :red
				@info @sprintf "q=%.3f, h=%.4f, dev=%.3f" row[1] row[2] row[3]
			end
		end
		res
	end
	plot_results(pts[1,:], pts[2,:], qts,
		filter_fn=filter_fn, color=colors,
		title="$dist", draw_surface=surf, mode=mode, az=az, el=el)
end

# ╔═╡ 2e3b2b76-c100-4f6d-b5a6-8263da6ca724
md"""
## Explore Specific Data Points

This section requires full simulation results (in `../data/`), not just proxy
values from `../data-proxy/`.
"""

# ╔═╡ 730899e1-edd5-4ddc-be98-b193a22bb6f8
begin
	i_99 = round(Int64, b * p)
	i_low, i_high = find_intervals(b, p, 0.05, centered=true)[1]
end

# ╔═╡ 87de4ba4-7f8e-4d8f-9bda-900df72339f9
let
	q = 0.70
	prob_con_miss(n, q) = (1-q)^n*q
	sum(n -> prob_con_miss(n, q), 0:3) - 0.99
end

# ╔═╡ 0c83a619-0405-4039-90fd-a09676eb9326
# ╠═╡ skip_as_script = true
#=╠═╡
let
	path = "../data/nmc-dist/"
	# q, h = 0.7, 0.0225
	q, h = 0.75, 0.0275
	# q, h = 0.79, 0.028
	filename = generate_filename(b, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i][2] for i in (i_low, i_99, i_high)]
end
  ╠═╡ =#

# ╔═╡ b34323cf-913f-4a1b-b385-f3ae52f320ac
let
	path="../data/nmc-dist/Pareto{Float64}(α=1.5, θ=0.01)"
	q = 0.70
	h = quantile(dist, q)
	filename = generate_filename(b, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i] for i in (i_low, i_99, i_high)]
	function count_con_zeros(bv)
		# Note: learn about why `let` makes this assignment
		# using outside scope by default unless prefixed with `local`
		local cnt = 0
		for b in bv
			if b
				return cnt
			end
			cnt += 1
		end
		return cnt
	end
	cnt = 0
	for i in i_99:b
		if count_con_zeros(data[i][1]) < 4
			cnt += 1
			# @info "Found" i data[i]
		end
	end
	println(cnt/b)
end

# ╔═╡ Cell order:
# ╟─0becd1ea-ade9-43d4-9039-fb335de62c6d
# ╟─305d3f52-0208-4fe7-b362-1f9b9424b705
# ╠═b216a1aa-dc98-11ee-0312-21d71fee5020
# ╟─9bd08151-ce7f-4d57-9dd7-7af79ff42d3f
# ╠═ce5e566b-b1e3-4e79-9156-4237a170546e
# ╟─4cdad6fd-9976-4bbf-b6de-95cb6a973fc7
# ╠═aca1dd0e-89b9-4a66-ae0b-2af07cea3636
# ╠═9c608fab-8fc9-4bf9-be68-96468257e6a8
# ╠═7b5d607e-da02-4060-8a9c-60c541fe34aa
# ╟─dde0af91-d0ff-4e4d-a4d4-f8fb770987dd
# ╟─d4549398-3ca3-44ba-b016-ca55bf4056cf
# ╠═d12e5d80-7099-4564-8144-778c812b106c
# ╠═5a426b8d-e062-43b2-81ef-77f8245292ee
# ╠═ae494460-6f76-4562-a6b4-42f0bc6df262
# ╟─2e3b2b76-c100-4f6d-b5a6-8263da6ca724
# ╠═730899e1-edd5-4ddc-be98-b193a22bb6f8
# ╠═87de4ba4-7f8e-4d8f-9bda-900df72339f9
# ╠═0c83a619-0405-4039-90fd-a09676eb9326
# ╠═b34323cf-913f-4a1b-b385-f3ae52f320ac
