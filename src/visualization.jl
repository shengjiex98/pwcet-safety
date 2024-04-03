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
	using Distributions: Distribution, Normal, Pareto, Uniform, cdf, pdf
	# using PlotlyJS
	using Plots
	# plotlyjs()

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
	
	points_small, quantiles_small = deserialize("points.jls")
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
	points_large, quantiles_large = deserialize("points_large.jls")
end

# ╔═╡ 3f0c133d-6482-4c1d-95cf-d7802799e2f8
begin
	points = [points_small points_large]
	quantiles = append!(quantiles_small, quantiles_large)
	# points = points_small
	# quantiles = quantiles_small
	# points = points_large
	# quantiles = quantiles_large
end

# ╔═╡ b6177e6a-7587-4fac-8af5-ebf946b8434d


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
	q = 0.8
	h = 0.03
	filename = generate_filename(b_large, q, h, th=16)
	data = deserialize("$path/$filename.jls")
	[data[i][2] for i in (i_low, i_99, i_high)]
end

# ╔═╡ 6eaac4e2-f3b4-4005-bd53-35232577d44d
function set_indices_to_1(bool_list::Vector{Bool}, indices::Vector{Int})
    new_bool_list = falses(length(bool_list))
    new_bool_list[indices] .= true
    return new_bool_list
end

# ╔═╡ 0699d6f6-caf7-4f74-a3f4-fe7f3e82e0fb
function divide(filtered_data::Vector{Tuple{Bool, Float64}}, quantiles::AbstractVector{Float64})
    # Extract bool_list from filtered_data
    bool_list = [x[1] for x in filtered_data]
    
    # Group indices of filtered_data by their second element (Float64)
    grouped_indices = Dict{Float64, Vector{Int}}()
    for (idx, (_, value)) in enumerate(filtered_data)
		if bool_list[idx]
        	if haskey(grouped_indices, value)
            	push!(grouped_indices[value], idx)
        	else
            	grouped_indices[value] = [idx]
        	end
		end
    end
	return grouped_indices
end

# ╔═╡ 9ce24ee4-8b94-4d2f-98a5-30302b61a9ef
function argmin_quantiles(filtered_data::Vector{Tuple{Bool, Float64}}, quantiles::AbstractVector{Float64})
	grouped_indices = divide(filtered_data, quantiles)
    # Calculate quantiles for each group and find the argmin
    argmins = Dict{Float64, Int}()
    for (value, _) in grouped_indices
        list = grouped_indices[value]
        min_index = argmin(quantiles[list])
        argmins[value] = list[min_index]
    end
    return argmins
end

# ╔═╡ ce5e566b-b1e3-4e79-9156-4237a170546e
function plot_results(filtered_data::Vector{Tuple{Bool, Float64}}, 		  points::Matrix{<:Real}, 
		quantiles::Vector{<:Real};
		title::String,
        cap::Real=Inf, draw_surface::Bool=true, draw_all_points::Bool=true, 	mode::String="free",
		az::Union{Real, AbstractRange{Real}}=55, el::Real=15)
	@boundscheck mode == "free" || mode == "q" || mode == "period" || 
	throw(ArgumentError("`mode` has to be `free`, `q`, or `period`. $mode is given"))
	if mode == "q"
		az, el = (0, 0)
	elseif mode =="period"
		az, el = (90, 0)
	end
	filtered = [x[1] for x in filtered_data]
	plot(xlabel="quantile", ylabel="period", zlabel="deviation",
		# xlims=(0, 1), ylims=(0, hs[end]), 
		# zlims=(0, min(cap, maximum(quantiles))),
		title=title, legend=:topleft, proj_type=:ortho)
    if draw_surface
        surface!(points[1,filtered], points[2,filtered], 
			# min.(cap, quantiles[filtered]))
			quantiles[filtered])
    end
	if draw_all_points
    	scatter!(points[1,filtered], points[2,filtered], 
			# min.(cap, quantiles[filtered]), 
			quantiles[filtered],
			label="")
		grouped_indices = divide(filtered_data, quantiles)
		for (key,_) in grouped_indices
			list = grouped_indices[key]
			group_bool = set_indices_to_1(filtered, list)
			scatter!(points[1,group_bool], points[2,group_bool], 
			# min.(cap, quantiles[filtered]), 
			quantiles[group_bool],
			label="$key")
		end
	end
	
    # Index of minimum deviation in filtered points.
    # min_id_in_filtered = argmin(quantiles[filtered])
    # Index of minimum deviation in all points.
	# println(min_id_in_filtered)
 	# min_id = findall(filtered)[min_id_in_filtered]
	# println(min_id)
	
	min_id_in_filtered = argmin_quantiles(filtered_data, quantiles)
	for (key, min_id) in min_id_in_filtered 
    	println("Combination with lowest deviation: " *
        	"q=$(points[1,min_id]), " *
        	"period=$(points[2,min_id]), " *
        	"dev=$(round(quantiles[min_id], sigdigits=3))")
    	scatter!([points[1,min_id]], 
       		[points[2,min_id]], 
        	[min(quantiles[min_id], cap)], color=:red,
        	label="")
	end
	list = collect(values(min_id_in_filtered))
	filtered_min = set_indices_to_1(filtered, list)
	plot!(points[1,filtered_min], points[2,filtered_min], 
		# min.(cap, quantiles[filtered]), 
		quantiles[filtered_min], color=:red,
		label="Pareto Front")
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
		draw_surface::Bool=true, draw_all_points::Bool=true, mode::String="free",
		az::Union{Real, AbstractRange{Real}}=55, el::Real=15)
	filtered = points[1,:] .<= cdf.(dist, points[2,:])
	plot_results(filtered, points, quantiles, title=title,
		cap=cap, draw_surface=draw_surface,draw_all_points=draw_all_points,  mode=mode, az=az, el=el)
end

# ╔═╡ 96ae6cbd-0a6d-4aaf-8881-1ff1471a8c9c
function plot_results(filter_fn::Function, points::Matrix{<:Real}, 
		quantiles::Vector{<:Real};
        cap::Real=Inf, title="cap=$cap",
		draw_surface::Bool=true, draw_all_points::Bool=true, mode::String="free",
		az::Union{Real, AbstractRange{Real}}=15, el::Real=15)
	filtered = filter_fn.(points[1,:], points[2,:], quantiles)
	plot_results(filtered, points, quantiles, title=title,
		cap=cap, draw_surface=draw_surface, draw_all_points=draw_all_points, mode=mode, az=az, el=el)
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
		(q < cdf(dist, h) &&
		qmin <= q <= qmax && 
		hmin <= h <= hmax &&
		dev <= cap, q)
	plot_results(filter_fn, points, quantiles, cap=cap, 
		title="", draw_surface=surf, draw_all_points=all_points,
	mode=mode, az=az, el=el)
end

# ╔═╡ 9c608fab-8fc9-4bf9-be68-96468257e6a8
let x = range(0, 0.06, 100)
	
	plot(x, cdf.(dist, x), title="CDF for $dist", xlabel="T")
end

# ╔═╡ Cell order:
# ╠═b216a1aa-dc98-11ee-0312-21d71fee5020
# ╠═55e5c7d9-969c-415e-8822-de80174e8710
# ╠═d8e29adb-6cf3-4e0c-a077-72c6606c220e
# ╠═3f0c133d-6482-4c1d-95cf-d7802799e2f8
# ╠═b6177e6a-7587-4fac-8af5-ebf946b8434d
# ╠═730899e1-edd5-4ddc-be98-b193a22bb6f8
# ╠═297e0882-2eb5-4f54-9b9f-91840659b34a
# ╠═0c83a619-0405-4039-90fd-a09676eb9326
# ╠═f8cbd064-2eff-4ce3-8a2b-a942d8ccb47d
# ╠═ce5e566b-b1e3-4e79-9156-4237a170546e
# ╠═6eaac4e2-f3b4-4005-bd53-35232577d44d
# ╠═0699d6f6-caf7-4f74-a3f4-fe7f3e82e0fb
# ╠═9ce24ee4-8b94-4d2f-98a5-30302b61a9ef
# ╠═a88efb9c-e111-42bc-86c4-93015193bd02
# ╠═96ae6cbd-0a6d-4aaf-8881-1ff1471a8c9c
# ╟─dde0af91-d0ff-4e4d-a4d4-f8fb770987dd
# ╟─d4549398-3ca3-44ba-b016-ca55bf4056cf
# ╠═ae494460-6f76-4562-a6b4-42f0bc6df262
# ╠═d0d10130-ef6e-4a7f-b3f9-4715f01b737c
# ╠═aca1dd0e-89b9-4a66-ae0b-2af07cea3636
# ╠═9c608fab-8fc9-4bf9-be68-96468257e6a8
