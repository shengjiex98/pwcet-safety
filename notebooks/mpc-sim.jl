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
	using Colors
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
	available_colors = [:red, :blue, :yellow, :green, :purple, :orange, :pink, :brown, :gray, :black, :burlywood, :chartreuse, :coral, :cyan, :mistyrose1, :grey, :lightyellow3, :darkolivegreen, :paleturquoise1, :thistle, :darkgray, :darkorange, :ivory, :khaki, :lawngreen, :lightsalmon, :linen, :magenta, :mediumaquamarine, :mintcream, :red, :blue, :yellow, :green, :purple, :orange, :pink, :brown, :gray, :black, :burlywood, :chartreuse, :coral, :cyan, :mistyrose1, :grey, :lightyellow3, :darkolivegreen, :paleturquoise1, :thistle, :darkgray, :darkorange, :ivory, :khaki, :lawngreen, :lightsalmon, :linen, :magenta, :mediumaquamarine, :mintcream, :red, :blue, :yellow, :green, :purple, :orange, :pink, :brown, :gray, :black, :burlywood, :chartreuse, :coral, :cyan, :mistyrose1, :grey, :lightyellow3, :darkolivegreen, :paleturquoise1, :thistle, :darkgray, :darkorange, :ivory, :khaki, :lawngreen, :lightsalmon, :linen, :magenta, :mediumaquamarine, :mintcream, :red, :blue, :yellow, :green, :purple, :orange, :pink, :brown, :gray, :black, :burlywood, :chartreuse, :coral, :cyan, :mistyrose1, :grey, :lightyellow3, :darkolivegreen, :paleturquoise1, :thistle, :darkgray, :darkorange, :ivory, :khaki, :lawngreen, :lightsalmon, :linen, :magenta, :mediumaquamarine, :mintcream]
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
	| low_cap 			   | $(@bind low_cap  Slider(500:1:3000, default=1050, show_value=true)) |
	| cap 			   | $(@bind cap  Slider(1000:1:18000, default=5000, show_value=true)) |
	| all points       | $(@bind all_points CheckBox(default=false)) |
	"""
end

# ╔═╡ a2a8f41a-9357-448b-a2a6-c1f68c2e2621
let
	files = 1:100
	len = length(files)
	plots_list = []
	color_names = available_colors[1:len]
	x = 1
	for file_num in files
		color = color_names[x]
		JOB_ID = 44669211
		y = readdlm("../data-proxy/mpc-grid/y/$JOB_ID-$file_num.csv", ',')
		x_y = y[:, 2]
		q_y = y[:, 1]
		y_y = y[:, 4]
		num_rows = size(y, 1)
		file_index = fill(file_num, num_rows)
		color_index_y = map(Int, file_index)
		hover_text_y = map(q_y, x_y, y_y, color_index_y) do q, p, h, dev
			@sprintf "q=%.3f period=%.5f dev=%.3f flag=%.3f" q p h dev
		end
		plot_y = scatter(x_y, y_y, ylim=(low_cap, cap), markercolor=color, legend=false, xlabel="time", ylabel="deviation", title="y", hover=hover_text_y)
		plot!(x_y, y_y, label="Pareto-Front y")
		plots = plot(plot_y, plot_title="MPC_$file_num")
		push!(plots_list, plots)
		x += 1
	end
	plots_list
end

# ╔═╡ 64fe734f-db91-4324-9200-c321cab08915
let
	files = 1:100
	len = length(files)
	matrix_y = []
	for file_num in files
		JOB_ID = 44669211
		y = readdlm("../data-proxy/mpc-grid/y/$JOB_ID-$file_num.csv", ',')
		num_rows = size(y, 1)
		file_index = fill(file_num, num_rows)
		y = hcat(y, file_index)
		if size(matrix_y, 1) == 0 || size(matrix_y, 2) == 0
			matrix_y = y
		end
		min_dev = Inf
		x = Inf
		for i in 1:size(y, 1)
			if (y[i, :][1] != 0) & (y[i, :][1] != 1.0e-5)
				if y[i, :][4] < min_dev
					min_dev = y[i, :][4]
					x = i
				end
			end
		end
		if(x != Inf)
			matrix_y = vcat(matrix_y, reshape(y[x, :], 1, :))
		end
	end
	matrix_y = matrix_y[101:end, :]
	sorted_indices = sortperm(matrix_y[:, 2])
	matrix_y = matrix_y[sorted_indices, :]
	print(matrix_y)
	x_y = matrix_y[:, 2]
	q_y = matrix_y[:, 1]
	y_y = matrix_y[:, 4]
	color_index_y = map(Int, matrix_y[:, 6])
	colors_y = [available_colors[i] for i in color_index_y]
	hover_text_y = map(q_y, x_y, y_y, color_index_y) do q, p, h, dev
		@sprintf "q=%.3f period=%.5f dev=%.3f flag=%.3f" q p h dev
	end
	plot_y = scatter(x_y, y_y, ylim=(0, cap),markercolor=colors_y, legend=false, xlabel="time", ylabel="deviation", title="y", hover=hover_text_y)
	plot!(x_y, y_y, ylim=(0, cap),label="Pareto-Front y")
end

# ╔═╡ 6c11abb1-46b0-4f9d-8677-75fcbf9d3d8f
let
	files = 1:100
	combined_matrix_ref = []
	combined_matrix_y = []
	JOB_ID = 44181059
	JOB_ID_2 = 44181668
	for file_num in files
		full_matrix_ref = readdlm("../data-proxy/mpc-flags-uniform-period/ref/$JOB_ID/$JOB_ID-$file_num.csv", ',')
		num_rows = size(full_matrix_ref, 1)
		file_index = fill(file_num, num_rows)
		matrix_ref = hcat(full_matrix_ref, file_index)
		if size(combined_matrix_ref, 1) == 0 || size(combined_matrix_ref, 2) == 0
			combined_matrix_ref = matrix_ref
		end
		combined_matrix_ref = vcat(combined_matrix_ref, matrix_ref)
		full_matrix_y = readdlm("../data-proxy/mpc-flags-uniform-period/y/$JOB_ID_2/$JOB_ID_2-$file_num.csv", ',')
		num_rows = size(full_matrix_y, 1)
		file_index = fill(file_num, num_rows)
		matrix_y = hcat(full_matrix_y, file_index)
		if size(combined_matrix_y, 1) == 0 || size(combined_matrix_y, 2) == 0
			combined_matrix_y = matrix_y
		end
		combined_matrix_y = vcat(combined_matrix_y, matrix_y)
	end
	sorted_matrix_ref = sortslices(combined_matrix_ref,dims=1)
	n_columns = size(sorted_matrix_ref, 2)
	ref = reshape(sorted_matrix_ref[1, :], 1, :)
	min_dev = sorted_matrix_ref[1, :][4]
    for i in 2:size(sorted_matrix_ref, 1)
		if sorted_matrix_ref[i, :][4] < min_dev
			min_dev = sorted_matrix_ref[i, :][4]
			ref = vcat(ref, reshape(sorted_matrix_ref[i, :], 1, :))
		end
	end
	sorted_matrix_y = sortslices(combined_matrix_y,dims=1)
	n_columns = size(sorted_matrix_y, 2)
	y = reshape(sorted_matrix_y[1, :], 1, :)
	min_dev = sorted_matrix_y[1, :][4]
    for i in 2:size(sorted_matrix_y, 1)
		if sorted_matrix_y[i, :][4] < min_dev
			min_dev = sorted_matrix_y[i, :][4]
			y = vcat(y, reshape(sorted_matrix_y[i, :], 1, :))
		end
	end
	x_ref = ref[:, 2]
	q_ref = ref[:, 1]
	y_ref = ref[:, 4]
	color_index_ref = map(Int, ref[:, 6])
	colors_ref = [available_colors[i] for i in color_index_ref]
	hover_text_ref = map(q_ref, x_ref, y_ref, color_index_ref) do q, p, h, dev
		@sprintf "q=%.3f period=%.5f dev=%.3f flag=%.3f" q p h dev
	end
	plot_ref = scatter(x_ref, y_ref, ylim=(0, cap), markercolor=colors_ref, legend=false, xlabel="time", ylabel="deviation", title="ref", hover=hover_text_ref)
	plot!(x_ref, y_ref, ylim=(0, cap),label="Pareto-Front ref")
	if all_points
		all_color_index_ref = map(Int, sorted_matrix_ref[:, 6])
		all_hover_text_ref = map(sorted_matrix_ref[:, 1], sorted_matrix_ref[:, 2], sorted_matrix_ref[:, 4], all_color_index_ref) do q, p, h, dev
			@sprintf "q=%.3f period=%.5f dev=%.3f flag=%.3f" q p h dev
		end
		scatter!(sorted_matrix_ref[:, 2], sorted_matrix_ref[:, 4], ylim=(0, cap), markercolor=:white, alpha=0.1,  legend=false, hover=all_hover_text_ref)
	end
	x_y = y[:, 2]
	q_y = y[:, 1]
	y_y = y[:, 4]
	color_index_y = map(Int, y[:, 6])
	colors_y = [available_colors[i] for i in color_index_y]
	hover_text_y = map(q_y, x_y, y_y, color_index_y) do q, p, h, dev
		@sprintf "q=%.3f period=%.5f dev=%.3f flag=%.3f" q p h dev
	end
	plot_y = scatter(x_y, y_y,  ylim=(0, cap), markercolor=colors_y, legend=false, xlabel="time", ylabel="deviation", title="y", hover=hover_text_y)
	plot!(x_y, y_y, ylim=(0, cap), label="Pareto-Front y")
	if all_points
		all_color_index_y = map(Int, sorted_matrix_y[:, 6])
		all_hover_text_y = map(sorted_matrix_y[:, 1], sorted_matrix_y[:, 2], sorted_matrix_y[:, 4], all_color_index_y) do q, p, h, dev
			@sprintf "q=%.3f period=%.5f dev=%.3f flag=%.3f" q p h dev
		end
		scatter!(sorted_matrix_y[:, 2], sorted_matrix_y[:, 4],  ylim=(0, cap), markercolor=:white, alpha=0.1,  legend=false, hover=all_hover_text_y)
	end
	plots = plot(plot_ref, plot_y, plot_title="Pareto Front for all MPC flags")
end

# ╔═╡ 2e99fe56-345c-4d2a-94ca-941fe663b924
# ╠═╡ disabled = true
#=╠═╡
let
	files = 1:30
	combined_matrix_ref = []
	combined_matrix_y = []
	JOB_ID = 44181059
	JOB_ID_2 = 44181668
	for file_num in files
		full_matrix_ref = readdlm("../data-proxy/mpc-flags-uniform-period/ref/$JOB_ID/$JOB_ID-$file_num.csv", ',')
		if size(combined_matrix_ref, 1) == 0 || size(combined_matrix_ref, 2) == 0
			num_rows = size(full_matrix_ref, 1)
			file_index = fill(file_num, num_rows)
			combined_matrix_ref = hcat(full_matrix_ref, file_index)
		end
		for i in 1:size(full_matrix_ref, 1)
			if full_matrix_ref[i, :][4] < combined_matrix_ref[i, :][4]
				combined_matrix_ref[i, :] = [full_matrix_ref[i, :]; file_num]
			end
		end
		full_matrix_y = readdlm("../data-proxy/mpc-flags-uniform-period/y/$JOB_ID_2/$JOB_ID_2-$file_num.csv", ',')
		if size(combined_matrix_y, 1) == 0 || size(combined_matrix_y, 2) == 0
			num_rows = size(full_matrix_y, 1)
			file_index = fill(file_num, num_rows)
			combined_matrix_y = hcat(full_matrix_y, file_index)
		end
		for i in 1:size(full_matrix_y, 1)
			if full_matrix_y[i, :][4] < combined_matrix_y[i, :][4]
				combined_matrix_y[i, :] = [full_matrix_y[i, :]; file_num]
			end
		end
	end
	x_ref = combined_matrix_ref[:, 1]
	y_ref = combined_matrix_ref[:, 4]
	color_index_ref = map(Int, combined_matrix_ref[:, end])
	colors_ref = [available_colors[i] for i in color_index_ref]
	hover_text_ref = map(x_ref, y_ref, color_index_ref) do q, h, dev
		@sprintf "q=%.3f dev=%.3f flag=%.3f" q h dev
	end
	plot_ref = scatter(x_ref, y_ref, markercolor=colors_ref, legend=false, xlabel="quantile", ylabel="deviation", title="ref", hover=hover_text_ref)
	
	x_y = combined_matrix_y[:, 1]
	y_y = combined_matrix_y[:, 4]
	color_index_y = map(Int, combined_matrix_y[:, end])
	colors_y = [available_colors[i] for i in color_index_y]
	hover_text_y = map(x_y, y_y, color_index_y) do q, h, dev
		@sprintf "q=%.3f dev=%.3f flag=%.3f" q h dev
	end
	plot_y = scatter(x_y, y_y, markercolor=colors_y, legend=false, xlabel="quantile", ylabel="deviation", title="y", hover=hover_text_y)
	plots = plot(plot_ref, plot_y, plot_title="Pareto Front for all MPC flags")
end
  ╠═╡ =#

# ╔═╡ b75abdc8-f875-48b5-9497-cf67b6c305fe
# ╠═╡ disabled = true
#=╠═╡
let
	files = 1:10
	len = length(files)
	plots_list = []
	color_names = available_colors[1:len]
	x = 1
	for file_num in files
		color = color_names[x]
		JOB_ID = 41671552
		full_matrix_ref = readdlm("../data-proxy/mpc-flags/ref/mpc-$JOB_ID-$file_num.csv", ',')
		colors = let
			res = fill(:lightblue, size(full_matrix_ref, 1))
			mindev = cap
			for (i, row) in enumerate(eachrow(full_matrix_ref))
				if row[4] < mindev
					mindev = row[4]
					res[i] = color
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
					res[i] = color
				end
			end
			res
		end
		fig2 = plot_results(full_matrix_y[:,1], full_matrix_y[:,2], full_matrix_y[:,4],
			confidence=[full_matrix_y[:,4]-full_matrix_y[:,3] full_matrix_y[:,5]-full_matrix_y[:,4]],
			filter_fn=filter_fn, title="y", draw_surface=surf,
			mode=mode, az=az, el=el, color=colors)
		plots = plot(fig1, fig2, plot_title="MPC_$file_num")
		push!(plots_list, plots)
		x += 1
	end
	plots_list
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═e37c5e34-0191-11ef-0c30-379698507b64
# ╠═39bbf368-4d9b-447a-be29-a0b74f391cf3
# ╠═566cc4d1-e984-44e4-ab06-4d05d903335a
# ╟─07e901e0-5a41-4b02-8bd7-c50e67022b70
# ╟─b18059b1-4eb3-486c-ad09-11cf386fff20
# ╠═07df95d3-8763-47bd-a6a4-3b5a42e3ccb4
# ╠═a2a8f41a-9357-448b-a2a6-c1f68c2e2621
# ╠═64fe734f-db91-4324-9200-c321cab08915
# ╠═6c11abb1-46b0-4f9d-8677-75fcbf9d3d8f
# ╠═2e99fe56-345c-4d2a-94ca-941fe663b924
# ╠═b75abdc8-f875-48b5-9497-cf67b6c305fe
