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
	using Serialization
	using Distributions
	# using PlotlyJS
	using Plots
	plotlyjs()

	push!(LOAD_PATH, "../lib")
	using Experiments
	
end

# ╔═╡ 55e5c7d9-969c-415e-8822-de80174e8710
begin
	path = "../data/nmc-samenom"

	b = 100_000
	p = 0.99
	
	qs = append!(collect(0.1:0.05:0.9), [0.99, 0.999])
	hs = collect(0.005:0.0025:0.06)
	# hs = append!(collect(0.005:0.0025:0.06), [0.1, 0.2, 0.5, 1.0])
	@info size(qs), size(hs)
	
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
		# points = Base.product(qs, hs) |> collect |> vec |> stack
		# quantiles = get_quantile.(b, points[1,:], points[2,:])
		# serialize("points.jls", (points, quantiles))
	else
		points, quantiles = deserialize("points.jls")
	end
end

# ╔═╡ 297e0882-2eb5-4f54-9b9f-91840659b34a
let
	q = 0.7
	h = 0.0225
	filename = generate_filename(b, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i][2] for i in (99_000, 98937, 99060)]
end

# ╔═╡ f8cbd064-2eff-4ce3-8a2b-a942d8ccb47d
let
	q = 0.8
	h = 0.03
	filename = generate_filename(b, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i][2] for i in (99_000, 98937, 99060)]
end

# ╔═╡ 0c83a619-0405-4039-90fd-a09676eb9326
let
	q = 0.75
	h = 0.0275
	filename = generate_filename(b, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i][2] for i in (99_000, 98937, 99060)]
end

# ╔═╡ 730899e1-edd5-4ddc-be98-b193a22bb6f8
find_intervals(b, p, 0.05, centered=true)

# ╔═╡ ce5e566b-b1e3-4e79-9156-4237a170546e
function plot_results(filtered::BitVector, points::Matrix{<:Real}, 
		quantiles::Vector{<:Real};
		title::String,
        cap::Real=Inf, draw_surface::Bool=true, mode::String="free",
		az::Union{Real, AbstractRange{Real}}=55, el::Real=15)
	@boundscheck mode == "free" || mode == "q" || mode == "period" || 
	throw(ArgumentError("`mode` has to be `free`, `q`, or `period`. $mode is given"))
	if mode == "q"
		az, el = (0, 0)
	elseif mode =="period"
		az, el = (90, 0)
	end
	plot(xlabel="quantile", ylabel="period", zlabel="deviation",
		xlims=(0, 1), ylims=(0, hs[end]), 
		zlims=(0, min(cap, maximum(quantiles))),
		title=title, legend=:topleft, proj_type=:ortho)
    if draw_surface
        surface!(points[1,filtered], points[2,filtered], 
			# min.(cap, quantiles[filtered]))
			quantiles[filtered])
    end
    scatter!(points[1,filtered], points[2,filtered], 
		# min.(cap, quantiles[filtered]), 
		quantiles[filtered],
		label="")
    # Index of minimum deviation in filtered points.
    min_id_in_filtered = argmin(quantiles[filtered])
    # Index of minimum deviation in all points.
    min_id = findall(filtered)[min_id_in_filtered]
    println("Combination with lowest deviation: " *
        "q=$(points[1,min_id]), " *
        "period=$(points[2,min_id]), " *
        "dev=$(round(quantiles[min_id], sigdigits=3))")
    scatter!([points[1,min_id]], 
        [points[2,min_id]], 
        [min(quantiles[min_id], cap)], 
        label="Minimum Deviation")
	if az isa Real
		plot!(camera=(az, el))
	else
		@gif for i in -50:10:150
			plot!(camera=(i,el))
		end fps=4
	end
end

# ╔═╡ a88efb9c-e111-42bc-86c4-93015193bd02
function plot_results(dist::Distribution, points::Matrix{<:Real}, 
		quantiles::Vector{<:Real};
        cap::Real=Inf, title="$dist cap=$cap",
		draw_surface::Bool=true, mode::String="free",
		az::Union{Real, AbstractRange{Real}}=55, el::Real=15)
	filtered = points[1,:] .<= cdf.(dist, points[2,:])
	plot_results(filtered, points, quantiles, title=title,
		cap=cap, draw_surface=draw_surface, mode=mode, az=az, el=el)
end

# ╔═╡ 96ae6cbd-0a6d-4aaf-8881-1ff1471a8c9c
function plot_results(filter_fn::Function, points::Matrix{<:Real}, 
		quantiles::Vector{<:Real};
        cap::Real=Inf, title="cap=$cap",
		draw_surface::Bool=true, mode::String="free",
		az::Union{Real, AbstractRange{Real}}=15, el::Real=15)
	filtered = filter_fn.(points[1,:], points[2,:], quantiles)
	plot_results(filtered, points, quantiles, title=title,
		cap=cap, draw_surface=draw_surface, mode=mode, az=az, el=el)
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
	"""
end

# ╔═╡ d0d10130-ef6e-4a7f-b3f9-4715f01b737c
let 
	if sv == "Save"
		savefig("illustration.pdf")
	end
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
	plot_results(filter_fn, points, quantiles, cap=cap, 
		title="", draw_surface=surf,
	mode=mode, az=az, el=el)
end

# ╔═╡ 9c608fab-8fc9-4bf9-be68-96468257e6a8
let x = range(0, 0.06, 100)
	
	plot(x, pdf.(dist, x), title="CDF for $dist", xlabel="T")
end

# ╔═╡ Cell order:
# ╠═b216a1aa-dc98-11ee-0312-21d71fee5020
# ╠═55e5c7d9-969c-415e-8822-de80174e8710
# ╠═297e0882-2eb5-4f54-9b9f-91840659b34a
# ╠═f8cbd064-2eff-4ce3-8a2b-a942d8ccb47d
# ╠═0c83a619-0405-4039-90fd-a09676eb9326
# ╠═730899e1-edd5-4ddc-be98-b193a22bb6f8
# ╠═ce5e566b-b1e3-4e79-9156-4237a170546e
# ╠═a88efb9c-e111-42bc-86c4-93015193bd02
# ╠═96ae6cbd-0a6d-4aaf-8881-1ff1471a8c9c
# ╟─dde0af91-d0ff-4e4d-a4d4-f8fb770987dd
# ╟─d4549398-3ca3-44ba-b016-ca55bf4056cf
# ╠═ae494460-6f76-4562-a6b4-42f0bc6df262
# ╠═d0d10130-ef6e-4a7f-b3f9-4715f01b737c
# ╠═aca1dd0e-89b9-4a66-ae0b-2af07cea3636
# ╠═9c608fab-8fc9-4bf9-be68-96468257e6a8
