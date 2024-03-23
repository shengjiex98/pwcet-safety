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

# ╔═╡ 0df1dd2a-e7ee-11ee-2b9f-05f3e0239116
begin
	import Pkg
	Pkg.activate("..")

	using Random: seed!
	using ControlSystems: lqr, time_scale
	using LinearAlgebra: I
	using Plots: plot
	using PlutoUI

	using ControlTimingSafety: evol, hold_kill
	
	push!(LOAD_PATH, "../lib")
	using Benchmarks: benchmarks, augment, c2d, delay_lqr
	using ContinuousSims: nominal_trajectory
	using Experiments: generate_samples
end

# ╔═╡ 553a9eb7-1394-4d00-8d73-226e37f3d59d
md"""
## Compare with previous implementation
"""

# ╔═╡ 337031c3-c5ec-4520-a29f-142d5030b27a
@bind period Slider(append!([0.00001, 0.0001, 0.001, 0.005], collect(0.01:0.01:0.1)), show_value=true)

# ╔═╡ 2f0e9332-f742-4631-9ce1-0698783caa87
F1 = time_scale(benchmarks[:F1T], period)

# ╔═╡ debd15e0-2807-46d4-a348-78d38495f247
K = lqr(F1, I, I)

# ╔═╡ c8bfacc3-386f-49d0-b754-b5b5fbff509e
a, z0 = let
	sys = benchmarks[:F1T]
	
	x0 = 1.
	u0 = 0.
	hold_kill(c2d(sys, period), delay_lqr(sys, period)),
		[fill(x0, size(sys.A, 1)); u0]
end

# ╔═╡ 302e495c-7550-480f-9a02-29a3fba6a79c
begin
	H_steps = floor(Int64, 100 * 0.02 / period)
end

# ╔═╡ 3b543ddc-d34d-4e8a-833f-edb1e0b8721e
z_nom = let u(x, t) = -K * x
	period = 1
	x0 = fill(1., F1.nx)
	nominal_trajectory(F1, u, period, period * H_steps, x0)
end

# ╔═╡ de095ebb-2569-43db-899c-bb2cb026deda
plot(0:H_steps, eachrow(z_nom), label=["x1" "x2" "u"])

# ╔═╡ c96499ee-adc9-4f83-a27b-316ca2a1e653
z_nom_d = evol(a, z0, 2 .- trues(H_steps))'

# ╔═╡ e10ece46-7be6-4b74-85b9-c1de072d0ff1
# plot(1:101, eachrow(F1.C * z_nom_d[1:2,:]))
# plot(1:101, [eachrow(F1.C * z_nom[1:2,:]), eachrow(F1.C * z_nom_d[1:2,:])], label=["c" "d"])
plot(1:H_steps+1, [eachrow(z_nom), eachrow(z_nom_d)], label=["x1 continuous" "x2 continuous" "u continuous" "x1 discrete" "x2 discrete" "u discrete"])

# ╔═╡ 0df2e281-a91a-4085-9443-1d28f35434c5
# ╠═╡ disabled = true
#=╠═╡
let 
	seed!(1234)
	samples1 = generate_samples(a, z0, 0.9, 15, H=H_steps)
end
  ╠═╡ =#

# ╔═╡ bd0d912d-2353-44bc-9646-dc23f7caf3f8
# ╠═╡ disabled = true
#=╠═╡
let 
	seed!(1234)
	samples2 = generate_samples(a, z0, 0.9, 15, H=H_steps, nominal_trajectory=z_nom_d)
end
  ╠═╡ =#

# ╔═╡ dbc9154e-a86a-4e72-9fe8-52b8ab7e47c1
# ╠═╡ disabled = true
#=╠═╡
let
	seed!(1234)
	samples3 = generate_samples(a, z0, 0.9, 15, H=H_steps, nominal_trajectory=z_nom)
end
  ╠═╡ =#

# ╔═╡ 80efb84e-61b4-4bdd-b5cf-15b26b10722d
# ╠═╡ disabled = true
#=╠═╡
let
	seed!(1234)
	samples3 = generate_samples(a, z0, 1.0, 15, H=H_steps, nominal_trajectory=z_nom)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═0df1dd2a-e7ee-11ee-2b9f-05f3e0239116
# ╠═2f0e9332-f742-4631-9ce1-0698783caa87
# ╠═debd15e0-2807-46d4-a348-78d38495f247
# ╠═3b543ddc-d34d-4e8a-833f-edb1e0b8721e
# ╠═de095ebb-2569-43db-899c-bb2cb026deda
# ╟─553a9eb7-1394-4d00-8d73-226e37f3d59d
# ╠═c8bfacc3-386f-49d0-b754-b5b5fbff509e
# ╠═c96499ee-adc9-4f83-a27b-316ca2a1e653
# ╠═337031c3-c5ec-4520-a29f-142d5030b27a
# ╠═302e495c-7550-480f-9a02-29a3fba6a79c
# ╠═e10ece46-7be6-4b74-85b9-c1de072d0ff1
# ╠═0df2e281-a91a-4085-9443-1d28f35434c5
# ╠═bd0d912d-2353-44bc-9646-dc23f7caf3f8
# ╠═dbc9154e-a86a-4e72-9fe8-52b8ab7e47c1
# ╠═80efb84e-61b4-4bdd-b5cf-15b26b10722d
