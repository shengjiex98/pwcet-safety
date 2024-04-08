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

# ╔═╡ b216a1aa-dc98-11ee-0312-21d71fee5020
begin
	import Pkg
	Pkg.activate("..")

	using PlutoUI
	using Printf
	using Serialization
	using DelimitedFiles
	using Distributions: Distribution, Normal, Pareto, Uniform, cdf, pdf
	# using PlotlyJS
	using Plots
	plotlyjs()

	push!(LOAD_PATH, "../src")
	using Experiments
	
end

# ╔═╡ 55e5c7d9-969c-415e-8822-de80174e8710
begin
	path = "../data/nmc-samenom"

	b = 100_000
	p = 0.99

	qs = append!(collect(0.1:0.05:0.9), [0.99, 0.999])
	hs = collect(0.005:0.0025:0.06)
	# hs = append!(collect(0.005:0.0025:0.06), [0.1, 0.2, 0.5, 1.0]
	@info "Number of quantiles and periods:" size(qs,1) size(hs,1)
	
	get_quantile(b, q, h; cap=Inf) = let
	    filename = generate_filename(b, q, h, th=16)
		if !isfile("$path/$filename.jls")
			@info "Parameters without valid data" b q h
			return Inf
		end
	    min(deserialize("$path/$filename.jls")[round(Int64, b * p)][2], cap)
	end

	readback = true
	if !readback
		let
			points = Base.product(qs, hs) |> collect |> vec |> stack
			quantiles = get_quantile.(b, points[1,:], points[2,:])
			serialize("points.jls", (points, quantiles))
		end
	end
	
	points_small, quantiles_small = deserialize("../data-proxy/points.jls")
end

# ╔═╡ d8e29adb-6cf3-4e0c-a077-72c6606c220e
begin
	b_large = 1_000_000

	readback_large = true
	if !readback_large
		let
			qs = collect(0.7:0.01:0.8)
			hs = append!(collect(0.02:0.001:0.03), [0.0225, 0.0275])
			@info "Number of quantiles and periods:" size(qs,1) size(hs,1)
			
			points = Base.product(qs, hs) |> collect |> vec |> stack
			quantiles = get_quantile.(b_large, points[1,:], points[2,:])
			serialize("points_large.jls", (points, quantiles))
		end
	end
	points_large, quantiles_large = deserialize("../data-proxy/points_large.jls")
end

# ╔═╡ 3f0c133d-6482-4c1d-95cf-d7802799e2f8
begin
	points = [points_small points_large]
	quantiles = vcat(quantiles_small, quantiles_large)
	# points = points_small
	# quantiles = quantiles_small
	# points = points_large
	# quantiles = quantiles_large
end

# ╔═╡ 730899e1-edd5-4ddc-be98-b193a22bb6f8
begin
	i_low, i_high = find_intervals(b_large, p, 0.05, centered=true)[1]
	i_99 = round(Int64, b_large*p)
end

# ╔═╡ 297e0882-2eb5-4f54-9b9f-91840659b34a
let
	q = 0.7
	h = 0.0225
	filename = generate_filename(b_large, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i][2] for i in (i_low, i_99, i_high)]
end

# ╔═╡ 0c83a619-0405-4039-90fd-a09676eb9326
let
	q = 0.75
	h = 0.0275
	filename = generate_filename(b_large, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i][2] for i in (i_low, i_99, i_high)]
end

# ╔═╡ f8cbd064-2eff-4ce3-8a2b-a942d8ccb47d
let
	q = 0.79
	h = 0.028
	filename = generate_filename(b_large, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i][2] for i in (i_low, i_99, i_high)]
end

# ╔═╡ d4549398-3ca3-44ba-b016-ca55bf4056cf
md"""
$(@bind go Button("Reset"))
$(@bind sv Button("Save"))
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
		az::Union{Real, AbstractRange{<:Real}}=55,
		el::Real=15)

	# Filter data if needed; `selection` is a BitVector representing points to show
	if !isnothing(filter_fn)
		selection = filter_fn.(q_values, h_values, quantiles)
	else
		selection = trues(length(q_values))
	end
	if color isa Vector
		color = color[selection]
	end
	# Set common properties
	plt = plot(title=title, legend=:topleft, proj_type=:ortho)
	
	# Show different graph depending on `mode`
	if mode == "q"
		if !isnothing(confidence)
			print(confidence[:,1][selection])
			print(confidence[:,2][selection])
            plot!(q_values[selection], quantiles[selection],
                    ribbon=(confidence[:,1][selection], confidence[:,2][selection]), 					fillalpha=0.2, color=:gray)
		end
		scatter!(q_values[selection], quantiles[selection], 
			xlim=(qmin, qmax), 
			# ylim=(0, maximum(quantiles[selection])),
			xlabel="quantile", ylabel="deviation", color=color)
	elseif mode == "period"
		scatter!(h_values[selection], quantiles[selection],
			xlim=(hmin, hmax), 
			# ylim=(0, maximum(quantiles[selection])),
			xlabel="period", ylabel="deviation", color=color)
	elseif mode == "free"
	    if draw_surface
	        surface!(q_values[selection], h_values[selection], quantiles[selection])
	    end
		scatter!(q_values[selection], h_values[selection], quantiles[selection],
			xlim=(qmin, qmax), ylim=(hmin, hmax), 
			# zlim=(0, maximum(quantiles[selection])),
			xlabel="quantile", ylabel="period", zlabel="deviation",
			color=color)
		if az isa Real
			plot!(camera=(az, el))
		else
			@gif for i in -50:10:150
				plot!(camera=(i,el))
			end fps=4
		end

	else
		throw(ArgumentError("`mode` has to be `free`, `q`, or `period`. $mode is given"))
	end
end

# ╔═╡ d12e5d80-7099-4564-8144-778c812b106c
let
	full_matrix = readdlm("../data-proxy/nmc-dist.csv", ',')
	colors = let
		res = fill(:lightblue, size(full_matrix, 1))
		mindev = cap
		for (i, row) in enumerate(eachrow(full_matrix))
			if row[4] < mindev
				mindev = row[4]
				res[i] = :red
				@info @sprintf "q=%.3f, h=%.5f, dev=%.3f" row[1] row[2] row[4]
			end
		end
		res
	end
	confidence = nothing
	# confidence = zeros(size(full_matrix, 1), 2) 
	# for (idx, elem) in enumerate(full_matrix[:, 1])
 #   		i_low, i_high = find_intervals(b_large, elem, 0.05, centered=true)[1]
	# 	filename = generate_filename(b_large, elem, full_matrix[:, 2][idx], th=16)
	# 	data = deserialize("$path/$filename.jls")
 #    	confidence[idx, :] = [data[i_low], data[i_high]]
	# end
	filter_fn = (q, h, dev) -> 
		qmin <= q <= qmax && 
		hmin <= h <= hmax &&
		dev <= cap
	plot_results(full_matrix[:,1], full_matrix[:,2], full_matrix[:,4],
		filter_fn=filter_fn, confidence = confidence, title="Pareto(1.5, 0.01)", draw_surface=surf,
		mode=mode, az=az, el=el, color=colors)
end

# ╔═╡ d0d10130-ef6e-4a7f-b3f9-4715f01b737c
if sv == "Save"
	savefig("illustration.pdf")
end

# ╔═╡ aca1dd0e-89b9-4a66-ae0b-2af07cea3636
# dist = Normal(0, 0)
# dist = Normal(0.02, 0.005)
dist = Pareto(1.5, 0.01)
# dist = Uniform(0, 0.06)

# ╔═╡ ae494460-6f76-4562-a6b4-42f0bc6df262
let
	filter_fn = (q, h, dev) -> 
		q < cdf(dist, h) &&
		qmin <= q <= qmax && 
		hmin <= h <= hmax &&
		dev <= cap
	colors = let
		res = fill(:lightblue, length(quantiles_large))
		mindev = cap
		for (i, row) in enumerate(zip(points_large[1,:], points_large[2,:], quantiles_large))
			if row[3] < mindev && row[1] < cdf(dist, row[2])
				mindev = row[3]
				res[i] = :red
				@info @sprintf "q=%.3f, h=%.4f, dev=%.3f" row[1] row[2] row[3]
			end
		end
		res
	end
	plot_results(points_large[1,:], points_large[2,:], quantiles_large,
		filter_fn=filter_fn, color=colors,
		title="$dist", draw_surface=surf, mode=mode, az=az, el=el)
end

# ╔═╡ 9c608fab-8fc9-4bf9-be68-96468257e6a8
let x = range(0, 0.06, 100)
	
	plot(x, cdf.(dist, x), title="CDF for $dist", xlabel="T")
end

# ╔═╡ Cell order:
# ╠═b216a1aa-dc98-11ee-0312-21d71fee5020
# ╠═55e5c7d9-969c-415e-8822-de80174e8710
# ╟─d8e29adb-6cf3-4e0c-a077-72c6606c220e
# ╠═3f0c133d-6482-4c1d-95cf-d7802799e2f8
# ╠═730899e1-edd5-4ddc-be98-b193a22bb6f8
# ╠═297e0882-2eb5-4f54-9b9f-91840659b34a
# ╠═0c83a619-0405-4039-90fd-a09676eb9326
# ╠═f8cbd064-2eff-4ce3-8a2b-a942d8ccb47d
# ╠═ce5e566b-b1e3-4e79-9156-4237a170546e
# ╟─dde0af91-d0ff-4e4d-a4d4-f8fb770987dd
# ╟─d4549398-3ca3-44ba-b016-ca55bf4056cf
# ╠═d12e5d80-7099-4564-8144-778c812b106c
# ╠═ae494460-6f76-4562-a6b4-42f0bc6df262
# ╠═d0d10130-ef6e-4a7f-b3f9-4715f01b737c
# ╠═aca1dd0e-89b9-4a66-ae0b-2af07cea3636
# ╠═9c608fab-8fc9-4bf9-be68-96468257e6a8
