### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ faf59350-8d67-11ee-0bdd-2510e986118b
begin
    import Pkg; Pkg.activate()
    using CSV, Interpolations, DataFrames, CairoMakie, DifferentialEquations, Turing, StatsBase, PlutoUI, Distributions, Optim, Dierckx, MakieThemes, Printf, Colors, Random, ColorSchemes
end

# ╔═╡ c8da63f7-fa15-4b2b-9018-40790723c6a7
md"📉 first, settings for plot theme."

# ╔═╡ 4391f124-cbef-46e5-8462-e4e5126f5b38
begin
	# see https://makieorg.github.io/MakieThemes.jl/dev/themes/ggthemr
	gg_theme = :flat 
	set_theme!(ggthemr(gg_theme))
	
	figsize = (0.9 * 500, 0.9 * 380)
	
	update_theme!(
		backgroundcolor="white",
		fontsize=20, 
		linewidth=4,
		markersize=14,
		# titlefont=AlgebraOfGraphics.firasans("Light"),
		size=figsize,
		Axis=(; backgroundcolor=
			(MakieThemes.GGThemr.ColorTheme[:flat][:background], 0.5)
		)
	)

	theme_colors = MakieThemes.GGThemr.ColorTheme[gg_theme][:swatch]
	colors = Dict(
		zip(
			["data", "model", "distn", "other"], 
			parse.(Color, theme_colors[1:4])
		)
	)
end

# ╔═╡ 245836a9-6b44-4639-9209-e7ad9035e293
TableOfContents()

# ╔═╡ 7752316d-9dd0-4403-aa08-22c977ff3727
md"""
# tank geometry and measurements

we characterize the cross-sectional area (from a helicopter view) of the liquid holding tank, as a function of height $h$. the cross-sectional area we model as a [rounded rectangle](https://mathworld.wolfram.com/RoundedRectangle.html). see notes in `data/tank_length_measurements.txt`.

💡 we can't measure the radius of the circle in the corner, so we measure perimeter instead and infer it.
"""

# ╔═╡ 76624080-150a-4783-b675-794365dcecee
md"📏 length-measurements"

# ╔═╡ cccf3dfb-8b3c-45e8-bb1c-e9579afc7e1a
@kwdef struct LengthMeasurements # units: cm
	# top and bottom base, rounded rectangles
	l̂_t::Float64 # length, l + 2 r
	ŵ_t::Float64 # width, 2 + 2r
 	p_t::Float64 # perimeter
	
	l̂_b::Float64
	ŵ_b::Float64
	p_b::Float64

	# height
	h_max::Float64
	
	# orifice height and radius
	hₒ::Float64
	rₒ::Float64
end

# ╔═╡ 0f8e955e-fe9e-4b9d-ac24-1d2fb403a46c
function l̂ŵp_to_lwr(l̂, ŵ, p)
	# x = [l, w, r]
	A = [
		1 0 2;
		0 1 2;
		2 2 2*π
	]
	b = [l̂, ŵ, p]
	x = A \ b
	l, w, r = x[1], x[2], x[3]
	return l, w, r
end

# ╔═╡ 0e485727-495c-444d-9fb4-f20bdaac2676
# our raw measurements
length_measurements = LengthMeasurements(
	# top cross-section
	14.6, 9.0, 44.3,
	# bottom cross-section
	13.4, 7.8, 40.1,
	# measured height of tank [cm]
	28.6,
	# height of orifice, cm
	0.9,
	# radius of orifice, cm
	5 / 64 * 2.54 / 2 # in diam -> cm diam -> radius
)

# ╔═╡ 85bac8de-3224-4f02-a625-fdf55ec611a1
l̂ŵp_to_lwr(
	length_measurements.l̂_b, length_measurements.ŵ_b, length_measurements.p_b
)

# ╔═╡ 6d94f549-645d-4c27-9e4f-046542b5fb16
begin
	struct TankGeometry
		# rounded rectangle (top)
		l_t::Float64
		w_t::Float64
		r_t::Float64 
		
		# rounded rectangle (bottom)
		l_b::Float64
		w_b::Float64
		r_b::Float64 

		# height of tank
		h_max::Float64

		# orifice height
		hₒ::Float64

		# orifice radius
		rₒ::Float64
	end
	
	function TankGeometry(lm::LengthMeasurements)
		# infer l, w, r of rounded rectangle (top and bottom)
		l_t, w_t, r_t = l̂ŵp_to_lwr(lm.l̂_t, lm.ŵ_t, lm.p_t)
		l_b, w_b, r_b = l̂ŵp_to_lwr(lm.l̂_b, lm.ŵ_b, lm.p_b)
		
		return TankGeometry(
			l_t, w_t, r_t,
			l_b, w_b, r_b,
			lm.h_max, lm.hₒ, lm.rₒ
		)
	end
end

# ╔═╡ 109a382d-8d41-4bc3-a23b-439a987b17c7
# area of rounded rectangle given length, width, radius
lwr_to_a(l, w, r) = l * w + # main rectangle
		2 * r * (l + w) + # four strips
		π * r ^ 2 # four circles

# ╔═╡ a6687107-7448-451e-a3cf-04a3d2c3d7a5
# cross-sectional area of water in the tank, from helicopter view, 
#   as a function of liquid level, h.
function A_of_h(h::Float64, tg::TankGeometry)
	# check for over/underflow
	h < 0.0      ? error("tank underflow!") : nothing
	h > tg.h_max ? error("tank overflow!")  : nothing

	# fraction tank is full
	θ = h / tg.h_max

	# l, w, r of rounded rectangle here
	#   note, these vary linearly with height.
	l = tg.l_b * (1 - θ) + tg.l_t * θ
	w = tg.w_b * (1 - θ) + tg.w_t * θ
	r = tg.r_b * (1 - θ) + tg.r_t * θ

	return lwr_to_a(l, w, r) # area of rounded rectange at this height
end

# ╔═╡ a391cd0a-f752-4efd-92de-43e7cec656d4
tank_geometry = TankGeometry(length_measurements)

# ╔═╡ dd4556a3-d6af-4aec-a59d-d6fcd4f4144d
lwr_to_a(tank_geometry.l_t, tank_geometry.w_t, tank_geometry.r_t)

# ╔═╡ c728c5d5-7aa8-437d-a6d8-8dc75974c86e
lwr_to_a(tank_geometry.l_b, tank_geometry.w_b, tank_geometry.r_b)

# ╔═╡ 9a7e5903-69be-4e0a-8514-3e05feedfed5
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1],
		xlabel="water height, h [cm]", 
		ylabel="cross-sectional area, A(h) [cm²]"
	)
	lines!(range(0, tank_geometry.h_max), 
		[A_of_h(hᵢ, tank_geometry) for hᵢ in range(0, tank_geometry.h_max)]
	)
	vlines!(tank_geometry.h_max, linestyle=:dash, color="gray")
	xlims!(0, nothing)
	ylims!(0, nothing)
	fig
end

# ╔═╡ 418525b7-c358-41da-b865-5df3feb15855
md"
# calibration of liquid level sensor

read in data characterizing the calibration curve of the liquid level sensor.
"

# ╔═╡ a95e371e-9319-4c7e-b5d9-4c4a50d12cd7
begin
	calibration_data = CSV.read(
		joinpath("data", "level_sensor_calibration.csv"), DataFrame
	)
	sort!(calibration_data, "level sensor reading")
end

# ╔═╡ 6ebe0cb0-ba35-411c-9a7a-a8b6eecf326f
md"compute the true liquid level."

# ╔═╡ ef43e50a-5af8-4733-88a4-cd159d173034
md"fit spline to calibration data to construct calibration curve."

# ╔═╡ 9dabad13-cfa4-4e06-950d-f7c7d96c1147
begin
	# map liquid level sensor reading to liquid level h.
	level_sensor_to_h = Spline1D(
		calibration_data[:, "level sensor reading"],
		calibration_data[:, "h [cm]"];
		k=2, s=4.0, bc="error"
	)
end

# ╔═╡ e040094c-7511-4831-b94a-1c1185868202
md"viz calibration data as well as interpolator."

# ╔═╡ 23ee0e85-a84b-4b63-b432-5526559efcee
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="level sensor reading", 
		ylabel="liquid level, h [cm]"
	)
	scatter!(
		calibration_data[:, "level sensor reading"], 
		calibration_data[:, "h [cm]"]
	)
	lsr = range(
		minimum(calibration_data[:, "level sensor reading"]),
		maximum(calibration_data[:, "level sensor reading"]),
		length=100
	)
	lines!(lsr, level_sensor_to_h.(lsr), color=Cycled(2))
	fig
end

# ╔═╡ 078c01f7-e47e-4af0-be1c-ac4527b735fd
md"
# time series data processing
"

# ╔═╡ 8b00d2b3-9182-42ab-8393-91707b813f60
function read_h_time_series(filename::String)
	#=
	read in file
	=#
	data = CSV.read(joinpath("data", filename), DataFrame, comment="#")
	rename!(data, 
		" liquid_level_reading" => "liquid_level_reading",
		"Time [s]" => "t [s]"
	)

	#=
	find starting point, t = 0, as when max water level occurs.
	=#
	# index of highest water level
	id_start = argmax(data[:, "liquid_level_reading"])
	# remove data before that
	data = data[id_start:end, :]
	# shift time
	data[:, "t [s]"] = data[:, "t [s]"] .- data[1, "t [s]"]

	#=
	apply calibration curve to get h
	=#
	data[:, "h [cm]"] = level_sensor_to_h.(data[:, "liquid_level_reading"])
	
	return select(data, ["t [s]", "h [cm]"])
end

# ╔═╡ 7899f488-9c48-466f-857d-f5a31b5820ab
md"""
## read in train/test time series data
"""

# ╔═╡ 96f26378-846c-4964-935c-0372e2e86e91
md"time series data for two experiments without any object in the tank"

# ╔═╡ 661feb84-339c-4bbe-a8a5-65de74ed58c8
all_train_data = read_h_time_series("liq_level_data_empty_train.csv")

# ╔═╡ 1e8a535e-25ea-490b-b545-e532c4fbc0f3
all_test_data = read_h_time_series("liq_level_data_empty_test.csv")

# ╔═╡ a0849611-23b3-4a91-a054-f390bc6c9f0a
md"""
## visualize experimental data
"""

# ╔═╡ 33f889d7-e875-40d8-9d6d-cc87b0fbaf22
function viz_data(data::DataFrame)
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="time, t [s]", 
		ylabel="water level, h [cm]"
	)
	scatter!(
		data[:, "t [s]"], 
		data[:, "h [cm]"],
		label="experiment",
		color=colors["data"]
	)
	xlims!(0, nothing)
	ylims!(0, nothing)
	fig
end

# ╔═╡ 0e4b3c36-2f09-405e-912b-22c893cd1715
viz_data(all_train_data)

# ╔═╡ ddd8f749-3126-4563-8177-4941b6b0447b
viz_data(all_test_data)

# ╔═╡ 710ddf42-2397-4d96-9a61-ff5c600ccd43
md"## downsampling"

# ╔═╡ 46084a31-e591-42f2-b8e6-02183ddfc6ac
@bind n_data_sample Select([20, 50, 100], default=20)

# ╔═╡ f21dc58e-d4e8-4314-b5dd-abbcb29efe86
function downsample(data::DataFrame, n::Int)
	# filter out sensor data that is below 508
	# this is the liquid level where the jet shooting outside the tank
	# dies, and it just runs down the tank (Toricelli's law invalid)
	data = filter("h [cm]" => (x -> x >= level_sensor_to_h(508)), data)
	
	ids = trunc.(Int, collect(range(1, nrow(data), length=n)))
	return data[ids, :]
end

# ╔═╡ 14b713b9-70c6-4506-a7b3-c67d021f8fce
train_data = downsample(all_train_data, n_data_sample)

# ╔═╡ 86b34265-6f2d-48d7-95ea-d10c8ae29ea9
viz_data(train_data)

# ╔═╡ 6bbb7a2c-464f-4e92-ab15-66745bac03ef
test_data = downsample(all_test_data, n_data_sample)

# ╔═╡ bcf21d4b-cc76-4b9e-8cdb-7a17eb2af605
viz_data(test_data)

# ╔═╡ e379461f-8896-4e2a-a71b-1871a8a37eb5
md"""
# dynamic model of the liquid level

```math
A(h)\frac{dh}{dt} = - c\pi r_{\rm hole}^2 \sqrt{2g [h(t)-h_{\rm hole}]}
```
"""

# ╔═╡ 7b7baa41-0185-4ed6-8fae-3a44e9912016
g = 980.665 # cm/s²

# ╔═╡ f2f236f4-59f3-4c05-811d-078cd04ddd79
md"""
## dynamic model 
"""

# ╔═╡ c6a263eb-cb45-4ee7-9c02-549c89298652
function f!(dh, h, params, t) # use in-place to prevent ODE error
	# if liquid level drops below height of orifice, none can leave.
	if h[1] <= params.hₒ
		dh[1] = 0.0
		return 0.0
	end
	dh[1] = - π * params.rₒ ^ 2 * params.c * 
		sqrt(2 * g * (h[1] .- params.hₒ)) / params.A_of_h(h[1])
end

# ╔═╡ 05ed4187-a01a-4a16-a0e7-b3867d252578
function simulate(h₀, params::NamedTuple, tf::Float64, callback=nothing)
	prob = ODEProblem(f!, [h₀], (0.0, tf), params)
	h_of_t = solve(
		prob, Tsit5(), saveat=0.5, 
		reltol=1e-6, abstol=1e-6, 
		callback=callback
	)
	return h_of_t
end

# ╔═╡ 6f7d4335-d9a6-4896-9d69-bfc1c2c1c3d0
md"""
## classical approach

### identify $c$
"""

# ╔═╡ 66815e8e-09d9-4b43-9f45-9379b3d34f78
function loss(data::DataFrame, c::Float64, tg::TankGeometry)	
	params = (
		# radius of the hole
		rₒ = tg.rₒ,
		# coefficient of determination
		c = c, 
		# height of the hole
		hₒ = tg.hₒ,
		# area as a function of h
		A_of_h = h -> A_of_h(h, tg)
	)

	h₀ = data[1, "h [cm]"]
	h_of_t = simulate(
		h₀, # h₀
		params,
		data[end, "t [s]"] # end of time
	)
	
	cost = 0.0
	for i in 1:nrow(data)
		tᵢ = data[i, "t [s]"]
		cost += (data[i, "h [cm]"] - h_of_t(tᵢ)[1]) ^ 2
	end
	
	return cost
end

# ╔═╡ 5e79d8e1-429c-414f-b3a6-8cf4b93d1336
md"maximum likelihood estimate of $c$"

# ╔═╡ 0b75c073-f167-4553-b746-539a14cfcf25
loss(train_data, 0.4, tank_geometry)

# ╔═╡ 4ed31219-9ce0-4f1b-8152-0002e64649ad
function compute_mle(data::DataFrame, tg::TankGeometry)
	res = optimize(c -> loss(data, c, tg), 0.2, 0.9)
	return Optim.minimizer(res)
end

# ╔═╡ 0763f1f3-dfee-4d1f-a934-bf387f9c80ff
c_opt = compute_mle(train_data, tank_geometry)

# ╔═╡ c57c808a-297c-4887-bf20-5ad0207d055e
params = (
	# radius of the hole
	rₒ = tank_geometry.rₒ,
	# coefficient of discharge
	c = c_opt,
	# height of the hole
	hₒ = tank_geometry.hₒ,
	# area as a function of h
	A_of_h = h -> A_of_h(h, tank_geometry)
)

# ╔═╡ d3307918-1fdb-4f87-bb92-67330d22e58b
sim_train_data = DataFrame(
	simulate(
		train_data[1, "h [cm]"], # h₀
		params,
		1.25 * train_data[end, "t [s]"] # end of time
	)
)

# ╔═╡ d25cda2f-6ec6-4c93-8860-f5ce9c3ee629
md"""
## visualize dynamic model (classically fit)
"""

# ╔═╡ 444f6d74-273e-486d-905a-1443ec0e98df
function viz_sim_fit(
	data::DataFrame, sim_data::DataFrame, c::Float64; 
	savename::Union{Nothing, String}=nothing
)
	 fig = Figure()
	 ax = Axis(
		fig[1, 1], 
		xlabel="time, t [s]", 
		ylabel="liquid level, h[cm]"
	)
	
	scatter!(
		data[:, "t [s]"], 
		data[:, "h [cm]"], 
		label="experiment"
	)
	
	lines!(sim_data[:, "timestamp"], sim_data[:, "value1"], 
		label="model", color=Cycled(2)
	)

	axislegend(@sprintf("c*=%.2f", c), titlefont=:regular)
	if ! isnothing(savename)
		save("paper/$savename.pdf", fig)
	end
	fig
end

# ╔═╡ 8cfdc784-4060-48b8-8d1a-3b8d11f7a9a7
viz_sim_fit(train_data, sim_train_data, c_opt, savename="classic_MLE_fit")

# ╔═╡ 6c010734-e8a0-4000-88eb-f2a85d25ed99
function viz_toy_h(sim_data::DataFrame; savename::String="toy_h")
	 fig = Figure(resolution=(300, 300), backgroundcolor=:transparent)
	 ax = Axis(
		fig[1, 1], 
		xlabel="time, t [s]", 
		ylabel="liquid level, h(t) [cm]",
		yticks=(
			[0, sim_data[1, "value1"]], 
			["0", rich("h", subscript("0"))]
			),
		xticks=[0]
	)
	hlines!(sim_data[end, "value1"], color="gray", linestyle=:dash, linewidth=1)
	hidedecorations!(ax, label=false, ticklabels=false)

	lines!(sim_data[:, "timestamp"], sim_data[:, "value1"], 
		label="model", color=Cycled(2)
	)
	xlims!(0, sim_data[end, "timestamp"])
	ylims!(0, nothing)
	
	save(savename * ".pdf", fig)
	fig
end

# ╔═╡ b7e27b45-78d4-41b5-9770-9632057413c6
viz_toy_h(sim_train_data)

# ╔═╡ a1a10e2f-1b78-4b93-9295-7c0055e32692
md"""
# Bayesian inference for model parameters
"""

# ╔═╡ 58eff13c-44b5-4f19-8a42-cf9907ac9515
@bind n_MC_sample Select([25, 50, 100, 250, 1000], default=25)

# ╔═╡ 68c9d88a-99b7-49be-9ac4-1e06c694c1a6
@bind n_chains Select([3, 5], default=3)

# ╔═╡ 8a21fa0f-d3c3-4aa2-8b8b-74001d921c4a
md"""
## infer model parameters for object-free experiment (training)
"""

# ╔═╡ c5754b4b-b576-4257-95d2-8888bbd063ec
const σₗ = 0.1 # cm [gives precision of our length measurements via tape]

# ╔═╡ 0bc7df52-4ac9-42ac-9094-ecaf3c27da31
const σᵣ = 0.001 # cm [precision of our drill]

# ╔═╡ 8f5b8859-6b8c-4f2a-af3a-b13c2d33fe2a
@model function forward_model(
	data::DataFrame, 
	lm::LengthMeasurements;
	prior_only::Bool=false
)
	σs_trunc = 2.0 
	#=
	prior distributions
	=#
	# length-measurements
	l̂_t ~ Normal(lm.l̂_t, σₗ)
	ŵ_t ~ Normal(lm.ŵ_t, σₗ)
	p_t ~ Normal(lm.p_t, σₗ)
	l_t, w_t, r_t = l̂ŵp_to_lwr(l̂_t, ŵ_t, p_t)

	l̂_b ~ Normal(lm.l̂_b, σₗ)
	ŵ_b ~ Normal(lm.ŵ_b, σₗ)
	p_b ~ Normal(lm.p_b, σₗ)
	l_b, w_b, r_b = l̂ŵp_to_lwr(l̂_b, ŵ_b, p_b)

	h_max ~ Normal(lm.h_max, σₗ)

	hₒ ~ Normal(lm.hₒ, σₗ)
	rₒ ~ Normal(lm.rₒ, σᵣ)

	function sampled_A_of_h(h)
		θ = h / h_max # fraction tank is full
		# calcualte l, w, r of rounded rectangle here.
		l = l_b * (1 - θ) + l_t * θ
		w = w_b * (1 - θ) + w_t * θ
		r = r_b * (1 - θ) + r_t * θ
		return lwr_to_a(l, w, r)
	end

	# discharge coefficient
	c ~ Truncated(Normal(0.65, 0.25), 0.0, 1.0) # unitless
	
	# measurement noise
	σ ~ Uniform(0.0, 0.5) # cm

	# initial liquid level
	h₀_obs = data[1, "h [cm]"] # cm
	h₀ ~ Truncated(
		Normal(h₀_obs, σ),
		0.0, h_max
	)

	# do not use the rest of the data if doing prior only.
	if prior_only
		return
	end
	
	# parameters for ODE solver
	params = (
			  rₒ=rₒ,
			  c=c,
			  hₒ=hₒ,
			  A_of_h=sampled_A_of_h
			)
	
	# set up and solve ODE
	h_of_t = simulate(h₀, params, 1.1 * data[end, "t [s]"])
	
	#=
	code up likelihood
	=#
	for i in 2:nrow(data) # start at 2 b/c IC handled with informative prior
		tᵢ = data[i, "t [s]"]
		ĥᵢ = h_of_t(tᵢ)[1]
		data[i, "h [cm]"] ~ Normal(ĥᵢ, σ)
	end

	return nothing
end

# ╔═╡ 8082559e-a5b0-41a8-b8ed-aec3b09e5b2b
begin
	nb_data_train_omit = 3 # surface tension prevents flow
	
	train_model = forward_model(
		train_data[1:end-nb_data_train_omit, :], length_measurements
	)
	
	train_posterior = DataFrame(
		sample(
			train_model, 
			NUTS(0.65), MCMCSerial(), n_MC_sample, n_chains; progress=true
		)
	)
end

# ╔═╡ b04ad0dc-10b5-433e-abc3-e87b4aa4f7eb
md"compute tank areas"

# ╔═╡ 2a973d9c-8b33-4c67-8177-73fd826c8dac
function infer_tank_radius_and_area!(data::DataFrame)
	for tb in ["_t", "_b"] # top or bottom
		lwrs = [
			l̂ŵp_to_lwr(row["l̂"*tb], row["ŵ"*tb], row["p"*tb])
			for row in eachrow(data)
		]
		
		data[:, "l"*tb] = [lwr[1] for lwr in lwrs]
		data[:, "w"*tb] = [lwr[2] for lwr in lwrs]
		data[:, "r"*tb] = [lwr[3] for lwr in lwrs]

		data[:, "a"*tb] = [
			lwr_to_a(row["l"*tb], row["w"*tb], row["r"*tb]) 
				for row in eachrow(data)
		]
	end
end

# ╔═╡ c31a2d3f-902b-4be9-a64c-b04cb83ffaa4
infer_tank_radius_and_area!(train_posterior)

# ╔═╡ d03cf081-8ebe-4f5b-a81f-abc843e1bb65
train_posterior

# ╔═╡ 2ee1ca40-141f-40ad-b4c1-a2e025f69f95
md"make sure never $h_0>H$."

# ╔═╡ c2d877b5-d309-4868-925d-dab8d7d23403
@assert all(train_posterior[:, "h_max"] .>= train_posterior[:, "h₀"])

# ╔═╡ c239deed-8291-45aa-95cf-94df26e0136d
md"""
## visualize posterior distribution of model parameters
"""

# ╔═╡ 7979b889-4782-45be-9a4f-91375f22f26f
params_to_title = Dict(
					"h_max" => rich("h", subscript("max")), 
	
					"a_b" => rich("a", subscript("b")), 
					"a_t" => rich("a", subscript("t")),
	
					"r_t" => rich("r", subscript("t")), 
					"l_t" => rich("l", subscript("t")),
					"w_t" => rich("w", subscript("t")),
					"r_b" => rich("r", subscript("b")), 
					"l_b" => rich("l", subscript("b")),
					"w_b" => rich("w", subscript("b")),
	
					"rₒ" => rich("r", subscript("o")),
					"c" => "c",
					"hₒ" => rich("h", subscript("o")),
					"h₀" => rich("h", subscript("0")),
					"σₗ" => rich("σ", subscript("ℓ")),
					"H" => rich("h", subscript("max")), 
					"σ" => "σ"
)

# ╔═╡ 73d702b9-cdc0-4ce9-802d-89443c8412ab
params_to_units = Dict(
	"a_b" => "cm²", "a_t" => "cm²", "rₒ" => "cm", "c" => "",
	"hₒ" => "cm", "h₀" => "cm", "h_max" => "cm", "σ" => "cm",
	"l_t" => "cm", "w_t" => "cm", "r_t" => "cm",
	"l_b" => "cm", "w_b" => "cm", "r_b" => "cm",
)

# ╔═╡ 5bb0b72a-8c77-4fcb-bbde-d144986d9c1e
function viz_posterior(
	posterior::DataFrame, params::Matrix{String},
	tg::TankGeometry, h₀_obs::Float64
)
	fig = Figure(size=(850, 600))
	axs = [Axis(fig[i, j]) for i = 1:3, j = 1:4]

	# l_t, w_t, r_t = l̂ŵp_to_lwr(lm.l̂_t, lm.ŵ_t, lm.p_t)
	# l_b, w_b, r_b = l̂ŵp_to_lwr(lm.l̂_b, lm.ŵ_b, lm.p_b)

	for i = 1:3
		for j = 1:4
			if j != 4
				colgap!(fig.layout, j, Relative(0.03))
			end
			p = params[i, j]
			if p == ""
				hidedecorations!(axs[i, j])
				hidespines!(axs[i, j])
				axs[i, j].backgroundcolor = "white"
				continue
			end

			if p in ["c", "σ"]
				println(p)
				println("\tstd: ", std(posterior[:, p]))
				println("\tmean: ", mean(posterior[:, p]))
			end

			# vizualize the distribution
			hist!(axs[i, j], posterior[:, p], color=colors["distn"])
		
			# plot equal-tailed 80% interval
			lo, hi = quantile(posterior[:, p], [0.1, 0.9])
			lines!(axs[i, j], [lo, hi], [1, 1], linewidth=8, color="black")
			
			axs[i, j].xlabel = rich(
				params_to_title[p], " [" * params_to_units[p] * "]"
			)
		
			if j == 1
				axs[i, j].ylabel = "# samples"
			else
				hideydecorations!(axs[i, j], grid=false)
			end

			p_obs = nothing
			if ! (p in ["h₀", "c", "σ"])
				p_obs = getfield(tg, Symbol(p))
				# plot measured value
			elseif p == "h₀"
				p_obs = h₀_obs
			end
			if ! isnothing(p_obs)
				vlines!(axs[i, j], p_obs, linestyle=:dash, color=Cycled(4), 
						label="measurement")
			end
			ylims!(axs[i, j], 0, nothing)

			axs[i, j].xticks = WilkinsonTicks(3)
			if p in ["σ", "c"]
				axs[i, j].xticks = LinearTicks(3)
			elseif p in ["rₒ", "h_max"]
				axs[i, j].xticks = LinearTicks(2)
			end

			μ, σ = mean(posterior[:, p]), std(posterior[:, p])
			if p == "rₒ"
				Label(fig[i, j], justification=:left,
					@sprintf("mean:\n\t%.2f %s\nSTD:\n\t%.3f %s", μ, params_to_units[p], σ, params_to_units[p]), font=:regular, fontsize=11.0,
					tellwidth=false, tellheight=false, 
					halign=p in ["a_b", "a_t", "h₀"] ? 0.025 : 0.99, 
					valign=0.9
				)
			else
				Label(fig[i, j], justification=:left,
					@sprintf("mean:\n\t%.2f %s\nSTD:\n\t%.2f %s", μ, params_to_units[p], σ, params_to_units[p]), font=:regular, fontsize=11.0,
					tellwidth=false, tellheight=false, 
					halign=p in ["a_b", "a_t", "h₀"] ? 0.025 : 0.99, 
					valign=0.9
				)
			end
		end
	end
	linkyaxes!(axs...)
	save("paper/posterior_train_theta.pdf", fig)
	return fig
end

# ╔═╡ d49936ff-c8c0-4a8d-a804-0fb56908b383
tank_geometry

# ╔═╡ ded5b462-06dd-43a4-93b0-c52ad87174eb
viz_posterior(
	train_posterior,
	["l_t" "w_t" "r_t" "hₒ"; "l_b" "w_b" "r_b" "h_max"; "rₒ" "h₀" "σ" "c"],
	tank_geometry,
	train_data[1, "h [cm]"]
)

# ╔═╡ 86b56683-c80e-4c0f-8b03-a4869860d04f
md"## posterior predictive check"

# ╔═╡ 323f3fd7-e9a9-4598-ad2e-c1790cf4a264
function mean_abs_residual(data::DataFrame, sim_data::DataFrame)
	# build h(t; θ)
	h_of_t = linear_interpolation(sim_data[:, "timestamp"], sim_data[:, "value1"])
	# compute mean residual
	r = 0.0
	for row in eachrow(data)
		tᵢ, hᵢ = row["t [s]"], row["h [cm]"]
		ĥᵢ = h_of_t(tᵢ)
		r += abs(hᵢ - ĥᵢ)
	end
	return r / nrow(data)
end

# ╔═╡ 2ab35999-3615-4f5c-8d89-36d77802fe9b
function viz_fit(
	posterior::DataFrame, data::DataFrame; 
	savename::Union{String, Nothing}=nothing, 
	n_sample::Int=50, only_ic::Bool=false,
	n_data_end_omit::Int=0, incl_legend::Bool=true, legend_pos=:lb
)
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="time, t [s]", 
		ylabel="water height, h [cm]"
	)
	
	ts = range(0, 1.05 * maximum(data[:, "t [s]"]), length=500)
	tspan = (0.0,  maximum(ts) * 1.05)
	
	# sample posterior models
	mar = 0.0 # mean absolute residual
	for i in sample(1:nrow(posterior), n_sample)
		# area of tank
		function a_of_h_tank(h)
			θ = h / posterior[i, "h_max"] # fraction tank is full
			# calcualte l, w, r of rounded rectangle here.
			l = posterior[i, "l_b"] * (1 - θ) + posterior[i, "l_t"] * θ
			w = posterior[i, "w_b"] * (1 - θ) + posterior[i, "w_t"] * θ
			r = posterior[i, "r_b"] * (1 - θ) + posterior[i, "r_t"] * θ
			return lwr_to_a(l, w, r)
		end

		# area of object
		if "sqrt_a_obj[1]" in names(posterior)
			N = sum(contains.(names(posterior), "sqrt_a_obj"))
			
			sqrt_a_obj = [posterior[i, "sqrt_a_obj[$n]"] for n in 1:N]
			hs = range(
				0.0, 0.999 * posterior[i, "h_max"], length=N
			)
			
			A_of_h_object = extrapolate(
				interpolate(
					hs, sqrt_a_obj .^ 2,
					FritschButlandMonotonicInterpolation()
				),
				Interpolations.Flat()
			)
		else
			A_of_h_object(h) = 0.0
		end
		
		params = (
			  rₒ=posterior[i, "rₒ"],
			  c=posterior[i, "c"],
			  hₒ=posterior[i, "hₒ"],
			  A_of_h=h -> a_of_h_tank(h) - A_of_h_object(h)
			)

		# set up, solve ODE
		sim_data = DataFrame(
			simulate(posterior[i, "h₀"], params, 1150.0)
		)
			
		lines!(
			sim_data[:, "timestamp"], sim_data[:, "value1"], 
			label="model", color=(colors["model"], 0.1)
		)

		# mean abs residual
		mar += mean_abs_residual(data, sim_data)

		# h hole
		hlines!(ax, posterior[i, "hₒ"], color=("gray", 0.1), linestyle=:dash, label="hₒ")
	end
	mar /= n_sample
	println("mean abs residual: [cm] ", mar)

	# data for fitting
	scatter!(
		data[only_ic ? 1 : 1:end-n_data_end_omit, "t [s]"], 
		data[only_ic ? 1 : 1:end-n_data_end_omit, "h [cm]"],
		label="data",
		color=colors["data"]
	)
	# data without fitting
	if n_data_end_omit > 0
		scatter!(
			data[end-n_data_end_omit:end, "t [s]"], 
			data[end-n_data_end_omit:end, "h [cm]"],
			strokewidth=1,
			strokecolor=colors["data"],
			color=("white", 0.0)
		)
	end
	if incl_legend
		axislegend(unique=true, position=legend_pos)
	end
	ylims!(0, 28.0)
	xlims!(0, 1150.0)
	if savename!=nothing
		save( "$savename.pdf", fig)
	end
	return fig
end

# ╔═╡ 2a01b228-f281-46c4-9764-fac6cc1b4217
viz_fit(train_posterior, train_data, savename="paper/posterior_train", n_data_end_omit=nb_data_train_omit, incl_legend=false)

# ╔═╡ 18b5c6d2-2230-4881-8066-51eff42125ae
train_data

# ╔═╡ a5ae695b-bfc0-4425-9b64-bbeeba7da015
md"""
## validate posterior with test data set
"""

# ╔═╡ eaf470e9-2898-41d5-a6d5-4cd846e9c0de
function viz_test(posterior::DataFrame, test_data::DataFrame;
				  savename::Union{String, Nothing}=nothing, n_sample::Int=100
)
	fig = Figure(size=(figsize[1], figsize[2] * 1.25))
	
	ax_stopping = Axis(
		fig[1, 1], 
		ylabel="# samples", 
		xticks=([],[]),
		height=100
	)
	ax = Axis(
		fig[2, 1], 
		xlabel="time, t [s]", 
		ylabel="water height, h [cm]",
	)

	# sample from the train posterior
	emptying_time = Float64[]
	mar = 0.0 # mean absolute residual
	for i in sample(1:nrow(posterior), n_sample)
		# check for first instance when the liquid level
		#  is the same as the height of the hole in the base
		condition(h, t, integrator) = h[1] - posterior[i, "hₒ"]
		
		# retrive the emptying time [t] when h(t) = h_hole
		function affect!(integrator)
			push!(emptying_time, integrator.t)
		end
		cb = ContinuousCallback(condition, affect!)

		function a_of_h_tank(h)
			θ = h / posterior[i, "h_max"] # fraction tank is full
			# calcualte l, w, r of rounded rectangle here.
			l = posterior[i, "l_b"] * (1 - θ) + posterior[i, "l_t"] * θ
			w = posterior[i, "w_b"] * (1 - θ) + posterior[i, "w_t"] * θ
			r = posterior[i, "r_b"] * (1 - θ) + posterior[i, "r_t"] * θ
			return lwr_to_a(l, w, r)
		end

		params = (
			  rₒ=posterior[i, "rₒ"],
			  c=posterior[i, "c"],
			  hₒ=posterior[i, "hₒ"],
			  A_of_h=a_of_h_tank
		)

		# sample an initial condtion
		h₀_obs = test_data[1, "h [cm]"]
		h₀_distn = Truncated(
			Normal(h₀_obs, posterior[i, "σ"]),
			0.0, posterior[i, "h_max"]
		)
		h₀ = rand(h₀_distn)

		# simulate trajectory
		sim_data = DataFrame(
			simulate(h₀, params, 1.25 * test_data[end, "t [s]"], cb)
		)
		
		# mean abs residual
		mar += mean_abs_residual(test_data, sim_data)
		
		# plot trajectories
		lines!(ax, sim_data[:, "timestamp"], sim_data[:, "value1"], 
			label="model", color=(colors["model"], 0.1))
		# h hole
		hlines!(ax, posterior[i, "hₒ"], color=("gray", 0.1), 
			linestyle=:dash)
	end
	mar /= n_sample
	println("mean abs residual: [cm] ", mar)
	
	# plot dist'n of emptying times
	hist!(ax_stopping, emptying_time, color=colors["other"], bins=12)
	# text!(ax_stopping, 1000, 11, 
	# 	text="emptying\ntime\ndistribution", align=(:center, :center),
	# 	fontsize=10, color=colors["other"]
	# )

	scatter!(
		ax,
		test_data[:, "t [s]"], 
		test_data[:, "h [cm]"],
		label="data",
		color=colors["data"]
	)
	# @assert 1.25 * test_data[end, "t [s]"] > maximum(emptying_time)
	xlims!(0, 1.01 * maximum(emptying_time))
	ylims!(ax, 0, 28.0)
	ylims!(ax_stopping, 0, nothing)

	linkxaxes!(ax, ax_stopping)
	axislegend(ax, unique=true)
	
	if ! isnothing(savename)
		save("$savename.pdf", fig)
	end
	
	return fig 
end	

# ╔═╡ a3ba0c9d-5f81-4023-9ce0-ff29536aa968
viz_test(train_posterior, test_data, savename="paper/test")

# ╔═╡ 67c46219-2183-47fd-bd3c-82facff98d53
md"## visualize prior too"

# ╔═╡ 38a03e22-d596-4227-a9de-3ef54dc7256e
begin
	train_model_prior = forward_model(train_data, length_measurements, prior_only=true)
	
	train_prior = DataFrame(
		sample(
			train_model_prior, 
			NUTS(0.65), MCMCSerial(), n_MC_sample, n_chains; progress=true
		)
	)
	
	infer_tank_radius_and_area!(train_prior)
end

# ╔═╡ 2148cbb2-b41f-4df3-ab5c-55a89eff7bf1
viz_fit(train_prior, train_data, savename="paper/prior_train", only_ic=true)

# ╔═╡ 968de6ad-eb48-4d70-b431-209e609904aa
md"## covariance matrix of posterior for reconstruction problem"

# ╔═╡ 63977532-9afa-454c-9f51-af6f4b238120
function compute_mean_cov(posterior_data::DataFrame, var_list::Array{String})
	Σ = cov(Matrix(posterior_data[:, var_list]))
	μ = mean(Matrix(posterior_data[:, var_list]), dims=1)[:]
	return μ, Σ
end

# ╔═╡ aa9c9d45-bb7c-4eee-af87-6fbc01df271d
var_list = [
	"h_max", 
	"l_t", "w_t", "r_t",
	"l_b", "w_b", "r_b",
	"rₒ", "hₒ", 
	"c", 
	"σ"
]

# ╔═╡ 868189ef-6a3f-425f-94ac-dbb0e1847b2e
var_list

# ╔═╡ bb0a7df4-7e84-472a-ab00-e3dd801daf8e
μ_train, Σ_train = compute_mean_cov(train_posterior, var_list)

# ╔═╡ 3f640581-edcc-4c7a-86ba-b168f31fe4a3
function viz_cov_matrix(
	Σ::Matrix{Float64}, var_list::Array{String}; incl_values::Bool=false
)
	cmap = ColorSchemes.diverging_cwm_80_100_c22_n256
	cbar_limits = (-1.0, 1.0)
	
	fig = Figure(size=(600, 600))
	ax = Axis(
		fig[1, 1], 
		xticks=(1:length(var_list), [params_to_title[p] for p in var_list]),
		xticklabelrotation=π/2,
		yticks=(1:length(var_list), [params_to_title[p] for p in var_list]),
		aspect=DataAspect(),
		title="correlation matrix, C",
		titlefont=:regular
	)
	hm = heatmap!(Σ, colormap=cmap, colorrange=cbar_limits)
	if incl_values
		for i = 1:size(Σ)[1]
			for j = 1:size(Σ)[1]
				text!(ax, i, j-.05, 
					text=@sprintf("%0.1e", Σ[i,j]), 
					align=(:center, :center), fontsize=6
				)
			end
		end
	end
	Colorbar(
		fig[1, 2], label="correlation", limits=cbar_limits, 
		colormap=cmap, 
		#highclip=cmap[end], lowclip=cmap[1],
		height=350
	)
	save("paper/posterior_corr_matrix.pdf", fig)
	fig
end

# ╔═╡ 3bb65b71-d191-498b-81bf-40ffff4df1f4
begin
	local C = cor(Matrix(train_posterior[:, var_list]))
	viz_cov_matrix(C, var_list)
end

# ╔═╡ 7c8608ba-759a-4bbf-be14-32813dcbf79b
Σ_train

# ╔═╡ 2fc78ec1-bf53-49da-b31a-6b5bf165eb81
function viz_mean_matrix(μ::Vector{Float64}, var_list::Array{String})
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xticks=([], []),
		yticks=(1:length(var_list), [params_to_title[p] for p in var_list]),
		aspect=DataAspect(),
		title=rich("mean vector, ", rich("m", font="bold")),
		titlefont=:regular
	)
	hm = heatmap!(reshape(μ, (1, length(μ))))
	Colorbar(fig[1, 2], hm, label="mean")
	# save("posterior_mean.pdf", fig)
	fig
end

# ╔═╡ 34f621c0-207e-41e5-8193-ad81d2d21a01
viz_mean_matrix(μ_train, var_list)

# ╔═╡ 29518390-5d57-4f2d-b617-d7699468caf9
[Σ_train[i, i] for i = 1:length(var_list)]

# ╔═╡ e041f6b5-92b1-46f0-b383-606625d59a4a
var_list

# ╔═╡ b4fbbef8-3389-4fe9-9ebf-7b356bedd705
@assert var(train_posterior[:, var_list[1]]) ≈ Σ_train[1, 1]

# ╔═╡ 22adb08e-4d8c-40c5-97ef-cb1f3f0f6d90
@assert mean(train_posterior[:, var_list[1]]) ≈ μ_train[1]

# ╔═╡ a515407b-f749-48f2-b8ea-62940a186cce
@assert mean(train_posterior[:, var_list[3]]) ≈ μ_train[3]

# ╔═╡ 692a7f74-ca65-494b-b683-2d30e34e4c1e
rand(MvNormal(μ_train, Σ_train)) # e.g. how to sample

# ╔═╡ 9533c662-80af-4dd4-bf25-02e894867360
md"""
# Bayesian inference of object shape

## read experimental data
"""

# ╔═╡ 8b6d766a-8f7b-4b9a-9a15-0f7375087120
all_data_w_object = read_h_time_series("liq_level_data_w_bottle.csv")

# ╔═╡ 807222ba-5ff8-4f33-a9a0-7c69b1dccf52
data_w_object = downsample(all_data_w_object, n_data_sample)[1:end-5, :] # cut end to avoid overfitting to h_hole

# ╔═╡ 16158266-36ed-44c3-a418-0c454955ce78
begin
	function viz_free_vs_occupied(data_w_object::DataFrame, data::DataFrame)
		fig = Figure()
		ax = Axis(
			fig[1, 1], 
			xlabel="time, t [s]", 
			ylabel="water level, h [cm]"
		)
		scatter!(
			data[:, "t [s]"], 
			data[:, "h [cm]"],
			label="no object",
			color=colors["data"]
		)
		scatter!(
			data_w_object[:, "t [s]"], 
			data_w_object[:, "h [cm]"],
			label="with object",
			color=colors["other"]
		)
		axislegend()
		ylims!(0, nothing)
		xlims!(0, nothing)
		# save("h_of_t_with_without_object.pdf", fig)
		fig
	end
	
	viz_free_vs_occupied(data_w_object, train_data)
end

# ╔═╡ 580de17a-625d-420e-974c-86766197025e
md"## measured area of object

a ground truth"

# ╔═╡ cb59f55b-c748-4a94-b344-e50a8fa7c690
begin
	object_true_area = CSV.read(joinpath("data", "bottle_area.csv"), DataFrame)
	rename!(object_true_area, "area " => "area [cm²]")
	select!(object_true_area, ["h [cm]", "area [cm²]"])
end

# ╔═╡ 89cced40-f24e-499e-8bfd-19c3964f689b
A_of_object = Spline1D(object_true_area[:, "h [cm]"], 
					   object_true_area[:, "area [cm²]"]; k=3, s=10.0, bc="zero")

# ╔═╡ b9515b3a-b254-49ae-8c2c-b8ce7ced4d3a
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="water level, h [cm]", 
		ylabel="cross-sectional area, Aₒ(h) [cm²]"
	)
	
	h_range = range(0.0, object_true_area[end, "h [cm]"], length=100)
	
	lines!(h_range, A_of_object.(h_range), color=Cycled(2))
	scatter!(object_true_area[:, "h [cm]"], object_true_area[:, "area [cm²]"])
	ylims!(0, nothing)
	xlims!(0, nothing)
	fig
end

# ╔═╡ a8861082-2214-45f1-bc49-733efe74c949
md"## simulate model with knowledge of true area

in practice, couldn't do this. just a check to see if this has any hope...
"

# ╔═╡ b59fa654-6946-4687-b14b-c2ef1f766f5c
object_params = (
		# area of the hole
		rₒ = length_measurements.rₒ,
		# fudge factor
		c = c_opt,
		# height of the hole
		hₒ = length_measurements.hₒ,
		# area as a function of h"
		A_of_h = h -> A_of_h(h, tank_geometry) - A_of_object(h)
	)

# ╔═╡ b12963ae-bf7d-4ef7-b1a8-e2d1e24f9b4b
sim_object_data = DataFrame(
	simulate(
		data_w_object[1, "h [cm]"], object_params, data_w_object[end, "t [s]"]
	)
)

# ╔═╡ cfbe753d-85a8-445f-9eda-14a376d7e0c6
viz_sim_fit(data_w_object, sim_object_data, c_opt)

# ╔═╡ 54e9eda2-d564-453a-8ea8-4c8395be9ed6
viz_toy_h(sim_object_data, savename="paper/toy_h_w_object")

# ╔═╡ c53edeef-324a-418f-907d-aaf557cb8d24
md"## classical (baseline) method


```math
[A(h)-A_o(h)]\frac{dh}{dt} = - c\pi r_{\rm hole}^2 \sqrt{2g [h(t)-h_{\rm hole}]}
```

thus
```math
A_o(h) = A(h) + \dfrac{c\pi r_{\rm hole}^2 \sqrt{2g [h(t)-h_{\rm hole}]}}{dh/dt}
```
(don't be alarmed at the sign, dh/dt < 0)

💡 fit splines to data, differentiate splines to get dh/dt.
"

# ╔═╡ 4831a3be-35d3-420c-8463-bb14a597cc6a
begin
	# fitting to data past drainage introduces spurious curvature
	data_w_object_spline_fit = filter(
		row -> row["h [cm]"] >= length_measurements.hₒ + 1.0, data_w_object
	)
	
	h_of_t_w_object = Spline1D(
		data_w_object_spline_fit[:, "t [s]"], data_w_object_spline_fit[:, "h [cm]"]; k=4, s=150.0, bc="nearest"
	)
end

# ╔═╡ 08cee0a9-d358-4954-8e46-74de23d48d86
t_end_classical = data_w_object[end-2, "t [s]"]

# ╔═╡ f3f886d6-3010-4dd9-b42a-d5309463beb6
h_of_t_w_object(30.0)

# ╔═╡ 5003d6c6-fa30-423e-80c1-a0f82e4085b9
derivative(h_of_t_w_object, 30.0) # WARNING: doesn't extrapolate.

# ╔═╡ 36b1822e-fe08-494b-a57d-5888163a7b54
function viz_spline_fit(
	data_w_object::DataFrame, h_of_t_w_object, h_hole::Float64
)
	fig = Figure(size=(500, 500))
	axs = [Axis(fig[i, 1]) for i = 1:2]
	linkxaxes!(axs...)
	axs[2].xlabel = "t [s]"
	axs[1].ylabel = "h(t) [cm]"
	axs[2].ylabel = "dh/dt [cm/s]"
	
	ts = range(0.0, t_end_classical, length=150)

	# top
	scatter!(
		axs[1], data_w_object[:, "t [s]"], data_w_object[:, "h [cm]"], label="data"
	)
	hs = h_of_t_w_object.(ts)
	lines!(axs[1], ts, hs, color=Cycled(3), label="spline fit")
	ylims!(axs[1], 0, nothing)
	axislegend(axs[1])

	# bottom
	h′ = [derivative(h_of_t_w_object, tᵢ) for tᵢ in ts]
	# flatten derivative since it doesn't handle extrapolation
	# h′[h′ .> 0.0] .= 0.0
	# shut off derivative
	lines!(axs[2], ts, h′, color=Cycled(3))
	
	for ax in axs
		xlims!(ax, 0, nothing)
	end
	save("paper/classical_spline_fit.pdf", fig)
	fig
end

# ╔═╡ 29ccc89b-f76a-44fa-8a89-e3ca10742ba1
viz_spline_fit(data_w_object, h_of_t_w_object, length_measurements.hₒ)

# ╔═╡ 3385b22e-85ef-4bb0-8b1b-d03411c89b4f
tank_geometry

# ╔═╡ 29d3cc8f-780b-449a-87ba-8d543ad2473b
function classical_soln_A_object(
	h_of_t_w_object,
	c_opt::Float64,
	tank_geometry::TankGeometry
)
	ts = range(0.0, t_end_classical, length=150)

	# dh/dt
	h′ = [derivative(h_of_t_w_object, tᵢ) for tᵢ in ts]

	# h
	hs = h_of_t_w_object.(ts)

	# inferred area in tank
	A_tank = [A_of_h(hᵢ, tank_geometry) for hᵢ in hs]
	α = A_tank .+ π * tank_geometry.rₒ ^ 2 * c_opt * 
			sqrt.(2 * g * (hs .- tank_geometry.hₒ)) ./ h′
	
	inferred_object_data = DataFrame("h [cm]" => hs, "α [cm²]" => α)
	
	# filter data where...
	filter!(row -> ! isinf(row["α [cm²]"]), inferred_object_data) # dh/dt=0
	filter!(row -> row["α [cm²]"] > 0.0, inferred_object_data) # α > 0

	sort!(inferred_object_data)
	
	return inferred_object_data
end

# ╔═╡ 0a48e016-2fba-47cc-a212-47b4a3324b20
classical_α = classical_soln_A_object(
	h_of_t_w_object, c_opt, tank_geometry
)

# ╔═╡ 57328e4f-e945-46eb-aba1-8b75dfeb4575
# reconstruction error
function reconstruction_error(
	object_true_area::DataFrame,
	classical_α::DataFrame,
	data_w_object_spline_fit::DataFrame
)		
	residual = 0.0
	n_pts = 0
	
	# estimate of radius
	r̂ = Spline1D(
		classical_α[:, "h [cm]"], sqrt.(classical_α[:, "α [cm²]"] / π)
	)
	for j = 1:nrow(object_true_area)
		# get true radius.
		hᵢ, aᵢ = object_true_area[j, "h [cm]"], object_true_area[j, "area [cm²]"]
		rᵢ = sqrt(aᵢ / π)
		
		# get predicted radius
		r̂ᵢ = r̂(hᵢ)

		# add to residual
		if hᵢ > minimum(data_w_object_spline_fit[:, "h [cm]"]) && hᵢ < maximum(data_w_object_spline_fit[:, "h [cm]"])
			n_pts += 1
			residual += abs(r̂ᵢ - rᵢ)
		end
	end
	@show n_pts
	return residual / n_pts
end

# ╔═╡ 8e785e48-f7ea-4d27-8977-f1879c8bd74b
reconstruction_error(
	object_true_area, classical_α, data_w_object_spline_fit
)

# ╔═╡ 5feb46c0-3888-4586-8b12-f990d4d38912
begin
	local fig = Figure(size=(500, 500))
	local ax = Axis(
		fig[1, 1], 
		xlabel="√(area), √(α) [cm]",
		ylabel="height, h [cm]",
		aspect=DataAspect()
	)

	# true bottle area
	scatterlines!(
		sqrt.(object_true_area[:, "area [cm²]"]), object_true_area[:, "h [cm]"],
		label="measured", color="black"
	)
	scatterlines!(
		-sqrt.(object_true_area[:, "area [cm²]"]), object_true_area[:, "h [cm]"],
		color="black"
	)

	# predicted area
	lines!(
		sqrt.(classical_α[:, "α [cm²]"]), classical_α[:, "h [cm]"], 
		label="predicted", color=theme_colors[8]
	)
	lines!(
		-sqrt.(classical_α[:, "α [cm²]"]), classical_α[:, "h [cm]"], 
		color=theme_colors[8]
	)
			
	# plot area of tank for reference
	local a_t = lwr_to_a(tank_geometry.l_t, tank_geometry.w_t, tank_geometry.r_t)
	local a_b = lwr_to_a(tank_geometry.l_b, tank_geometry.w_b, tank_geometry.r_b)

	r_tank = sqrt.([a_b, a_t])
	lines!(r_tank, [0, tank_geometry.h_max], color="gray", label="tank")
	lines!(-r_tank, [0, tank_geometry.h_max], color="gray")
	lines!([-r_tank[2], r_tank[2]], [tank_geometry.h_max, tank_geometry.h_max], color="gray")
	lines!([-r_tank[1], r_tank[1]], [0, 0], color="gray")
	
	# xlims!(0, nothing)
	ylims!(0, nothing)
	axislegend(position=:cb)
	save("paper/classical_soln.pdf", fig)
	fig
end

# ╔═╡ b23dc763-d91f-4d66-94d2-dcf96cb07f54
md"## Bayesian inference

### forward model

play with interpolator used to represent the $\alpha(h)$ function.
"

# ╔═╡ 8897acea-5efb-47a6-83a2-0c70fccfdb46
begin
	x_i = range(0.0, 5.0, length=10)
	y_i = x_i .^ 2 .+ 3.0 .+ 3 * sin.(x_i) .+ exp.(x_i)
	
	itp = extrapolate(
		interpolate(x_i, y_i, FritschButlandMonotonicInterpolation()), Line()
	)
	x_d = range(0.0, 5.2, length=100)
	local fig = Figure()
	local ax = Axis(fig[1, 1])
	lines!(x_d, itp.(x_d))
	scatter!(x_i, y_i)
	fig
	# sitp = Interpolations.scale(itp, x_i)
end

# ╔═╡ 798d8d16-1c19-400d-8a94-e08c7f991e33
@model function forward_model_object(
	data_w_object::DataFrame, train_posterior::DataFrame, 
	N::Int, lm::LengthMeasurements, γ::Float64, var_list::Vector{String};
	prior_only::Bool=false
)
	#=
	prior distributions
	yesterday's posterior is today's prior
	=#
	μ_pr, Σ_pr = compute_mean_cov(train_posterior, var_list)
	# sample param vector from prior
	Ψ ~ MvNormal(μ_pr, Σ_pr)
	l_t    = Ψ[1] # hard-coded based on new_var_list: warning!
	w_t    = Ψ[2]
	r_t    = Ψ[3]
	l_b    = Ψ[4]
	w_b    = Ψ[5]
	r_b    = Ψ[6]
	h_max  = Ψ[7]
	rₒ     = Ψ[8]
	hₒ     = Ψ[9]
	c      = Ψ[10]
	σ      = Ψ[11]
	
	# initial liquid level
	h₀_obs = data_w_object[1, "h [cm]"]
	h₀ ~ Truncated(
		Normal(h₀_obs, σ), 0, h_max
	)

	# tank geometry
	function a_of_tank(h)
		θ = h / h_max # fraction tank is full
		# calcualte l, w, r of rounded rectangle here.
		l = l_b * (1 - θ) + l_t * θ
		w = w_b * (1 - θ) + w_t * θ
		r = r_b * (1 - θ) + r_t * θ
		return lwr_to_a(l, w, r)
	end

	# solid geometry
	hs = range(0.0, 0.999 * h_max, length=N)
	sqrt_a_obj = Vector{Float64}(undef, N)
	a_b = lwr_to_a(l_b, w_b, r_b)
	sqrt_a_obj[1] ~ Uniform(0.0, sqrt(a_b))
	for i = 2:N
		sqrt_a_obj[i] ~ Truncated(
			Normal(sqrt_a_obj[i - 1], γ),
			0.0, sqrt(a_of_tank(hs[i]))
		)
	end

	# for prior, do not show the algo the data :)
	if prior_only
		return nothing
	end

	#=
	set up dynamic model for h(t)
	=#
	a_of_object = extrapolate(
		interpolate(hs, sqrt_a_obj .^ 2, FritschButlandMonotonicInterpolation()), Interpolations.Flat()
	)
	
	# parameter for ODE solver
	params = (
			  rₒ=rₒ,
			  c=c,
			  hₒ=hₒ,
			  A_of_h=h -> a_of_tank(h) - a_of_object(h)
			)

	# checks before simulation
	@assert all([params.A_of_h(hᵢ) ≥ 0 for hᵢ in hs])
	# @assert h₀ ≤ H
	# @assert r_hole > 0.0
	# @assert σ > 0.0
	# @assert h_hole > 0.0
	# @assert c > 0.0
	
	# set up, solve ODE
	h_of_t = simulate(h₀, params, 1.01 * data_w_object[end, "t [s]"])

	# observations.
	for i in 2:nrow(data_w_object)
		tᵢ = data_w_object[i, "t [s]"]
		ĥᵢ = h_of_t(tᵢ)[1]
		data_w_object[i, "h [cm]"] ~ Normal(ĥᵢ, σ)
	end
	
	return nothing
end

# ╔═╡ da44647a-36e4-4116-9698-df1cb059c2b7
md"### posterior"

# ╔═╡ fb3ece76-f85c-41e1-a332-12c71d9d3cc0
N = 10  # number of points to infer area on

# ╔═╡ 743e74ec-9c66-4665-845d-75ede418616b
@bind infer_shape CheckBox()

# ╔═╡ 1aca6b92-7754-4cb3-b9e8-5d486e3bfcf8
begin
	if infer_shape
		new_var_list = [
			"l_t", "w_t", "r_t",
			"l_b", "w_b", "r_b",
			"h_max", "rₒ", "hₒ", "c", "σ"
		] # DO NOT CHANGE unless you change forward_model_object too.
		
		γ = 1.0 # smoothness param
		nb_data_object_omit = 2 # surface tension prevents flow
		
		object_tank_model = forward_model_object(
			data_w_object[1:end-nb_data_object_omit, :], train_posterior, N, length_measurements, γ, new_var_list
		)
		
		object_posterior = DataFrame(
			sample(object_tank_model, NUTS(0.65), MCMCSerial(), 
				n_MC_sample, n_chains, progress=true
			)
		)
		rename!(
			object_posterior, 
			["Ψ[$i]" => new_var_list[i] for i = 1:length(new_var_list)]...
		)

		# compute top, bottom area for posterior
		for tb in ["_t", "_b"]
			object_posterior[:, "a"*tb] = lwr_to_a.(
				object_posterior[:, "l"*tb], 
				object_posterior[:, "w"*tb], 
				object_posterior[:, "r"*tb]
			)
		end
	end
end

# ╔═╡ 57b18cdd-27a8-44df-b959-f5e7c7eb7413
for tb in ["_t", "_b"]
	object_posterior[:, "a"*tb] = lwr_to_a.(
		object_posterior[:, "l"*tb], 
		object_posterior[:, "w"*tb], 
		object_posterior[:, "r"*tb]
	)
end

# ╔═╡ 3c9a219f-74ef-45fb-83e7-c497e0bee362
@assert all(object_posterior[:, "h₀"] .< object_posterior[:, "h_max"])

# ╔═╡ e1264f57-f675-4f37-b4db-313cfc52ab8e
viz_fit(object_posterior, data_w_object, savename="paper/posterior_object", n_data_end_omit=nb_data_object_omit, incl_legend=false)

# ╔═╡ 7127fc35-a0af-4463-9448-a948f229fd47
function viz_inferred_radius(
	object_posterior::DataFrame, 
	object_true_area::DataFrame, 
	length_measurements::LengthMeasurements;
	savename::Union{Nothing, String}=nothing,
	show_legend::Bool=true,
	viz_measurements::Bool=true,
	n_sample::Int=50
)
	N = sum(contains.(names(object_posterior), "sqrt_a_obj"))
	
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="√(area), √(α) [cm]",
		ylabel="height, h [cm]",
		aspect=DataAspect()
	)

	residuals = zeros(nrow(object_posterior), nrow(object_true_area))
	fig_ia = nothing, nothing
	n_plotted = 1
	for i in shuffle(1:nrow(object_posterior))
		# unpack samples
		h_max = object_posterior[i, "h_max"]
		a_t, a_b = object_posterior[i, "a_t"], object_posterior[i, "a_b"]
		hₒ = object_posterior[i, "hₒ"]
		h₀ = object_posterior[i, "h₀"]
		
		# characterize A(h) and r(h)
		hs = range(
			0.0, h_max, length=N
		)
		sqrt_a_obj = [object_posterior[i, "sqrt_a_obj[$n]"] for n in 1:N]
		a_of_object = extrapolate(
			interpolate(hs, sqrt_a_obj .^ 2, FritschButlandMonotonicInterpolation()), Interpolations.Flat()
		)
		
		# compute residuals
		for j = 1:nrow(object_true_area)
			hᵢ, aᵢ = object_true_area[j, "h [cm]"], object_true_area[j, "area [cm²]"]
			r̂ᵢ = sqrt(a_of_object(hᵢ) / π)
			if hᵢ > mean(object_posterior[:, "hₒ"]) && hᵢ < mean(object_posterior[:, "h₀"])
				residuals[i, j] += abs(sqrt(aᵢ / π) - r̂ᵢ)
			else
				residuals[i, j] = NaN
			end
		end

		if n_plotted <= n_sample
			n_plotted += 1
			# plot inferred area of the object
			hs_dense = range(0.0, h_max, length=100)
			sqrt_a_obj_dense = sqrt.(a_of_object.(hs_dense))
			fig_ia = lines!(
				sqrt_a_obj_dense, hs_dense, label="model", color=(theme_colors[8], 0.1)
			)
			lines!(
				-sqrt_a_obj_dense, hs_dense, label="model", color=(theme_colors[8], 0.1)
			)
			
			# plot area of tank for reference
			r_tank = sqrt.([a_b, a_t])
			lines!(r_tank, [0, length_measurements.h_max], color="gray")
			lines!(-r_tank, [0, length_measurements.h_max], color="gray")
			lines!([-r_tank[2], r_tank[2]], [length_measurements.h_max, length_measurements.h_max], color="gray")
			lines!([-r_tank[1], r_tank[1]], [0, 0], color="gray")
		end
	end

	# plot hₒ and h₀
	hₒ = mean(object_posterior[:, "hₒ"])
	h₀ = mean(object_posterior[:, "h₀"])
	a_t = mean(object_posterior[:, "a_t"])
	a_b = mean(object_posterior[:, "a_b"])
	h_max = mean(object_posterior[:, "h_max"])
	sqrt_a_of_h(h) = sqrt(h / h_max * a_t + (1 - h / h_max) * a_b)
	
	lines!([sqrt_a_of_h(hₒ), 1.1 * sqrt_a_of_h(hₒ)], [hₒ, hₒ],
		color="blue", label="hₒ"
	)
	text!([sqrt_a_of_h(hₒ)*1.1], [hₒ], 
		text=rich("h", subscript("o")), align=(:left, :center)
	)

	
	lines!([-sqrt_a_of_h(h₀), sqrt_a_of_h(h₀)], [h₀, h₀],
		color="blue", label="h₀", linestyle=:dash
	)
	text!([sqrt_a_of_h(h₀)*1.1], [h₀], 
		text=rich("h", subscript("0")), align=(:left, :center)
	)

	if viz_measurements
		# measured area
		fig_ma = scatterlines!(
			sqrt.(object_true_area[:, "area [cm²]"]), 
			object_true_area[:, "h [cm]"], markersize=10,
			label="measurment", color=colors["data"]
		)
		scatterlines!(
			-sqrt.(object_true_area[:, "area [cm²]"]), 
			object_true_area[:, "h [cm]"], markersize=10,
			label="measurment", color=colors["data"]
		)
	end

	my_xlim = 1.35 * sqrt(mean(object_posterior[:, "a_t"]))
	xlims!(-my_xlim, my_xlim)
	ylims!(-1, length_measurements.h_max * 1.05)

	if show_legend
		# axislegend(#"γ=$γ; N=$N", 
		# 	unique=true, position=(0.8, 0.8), titlefont="normal", labelsize=16)
		Legend(fig[1, 2], [fig_ma, fig_ia], ["measured", "posterior"])
		colgap!(fig.layout, 1)
	end

	println("mean residual radius: ", mean(residuals[.! isnan.(residuals)]))
	println("avg actual radius: ", 
		mean(sqrt.(object_true_area[:, "area [cm²]"] ./ π))
	)
	
	if ! isnothing(savename)
		save("$savename.pdf", fig)
	end
	
	return fig
end

# ╔═╡ 40157899-dffb-4e3a-b5ca-be3c23a465ae
viz_inferred_radius(
	object_posterior, object_true_area, length_measurements, savename="paper/posterior_area"
)

# ╔═╡ 0c624592-df44-413e-9101-eefc3913d658
object_posterior

# ╔═╡ 3d0c3999-77a2-40d6-923b-78d3329e2154
0.31/3.22

# ╔═╡ bd95428d-1077-4417-bfca-0c5da7378af2
md"### prior"

# ╔═╡ 65d81268-9ff2-4a18-b0ce-4b105740dc8b
begin
	object_tank_model_prior = forward_model_object(
		data_w_object, train_posterior, N, length_measurements, γ, new_var_list, prior_only=true
	)
	
	object_prior = DataFrame(
		sample(object_tank_model_prior, NUTS(0.65), MCMCSerial(), 
			n_MC_sample, n_chains; progress=true
		)
	)
	
	rename!(
		object_prior, 
		["Ψ[$i]" => new_var_list[i] for i = 1:length(new_var_list)]...
	)
	# compute top, bottom area for posterior
	for tb in ["_t", "_b"]
		object_prior[:, "a"*tb] = lwr_to_a.(
			object_prior[:, "l"*tb], 
			object_prior[:, "w"*tb], 
			object_prior[:, "r"*tb]
		)
	end
end

# ╔═╡ 8c1d1401-bc6b-4be3-8481-1c9a8f86f63d
viz_inferred_radius(object_prior, object_true_area, length_measurements, savename="paper/prior_area", show_legend=false, viz_measurements=false)

# ╔═╡ 5b9d558a-2991-489e-be58-f5a5db0479f8
viz_fit(object_prior, data_w_object, savename="paper/prior_object", only_ic=true, legend_pos=:rt)

# ╔═╡ 6d4b0c74-4228-41e4-a8d0-98e0d71333b9
hist(object_prior[:, "sqrt_a_obj[1]"]) # check prior

# ╔═╡ 67b3c66c-b3eb-438f-96e4-c09e117cde87
lines(object_prior[:, "sqrt_a_obj[1]"])

# ╔═╡ Cell order:
# ╠═faf59350-8d67-11ee-0bdd-2510e986118b
# ╟─c8da63f7-fa15-4b2b-9018-40790723c6a7
# ╠═4391f124-cbef-46e5-8462-e4e5126f5b38
# ╠═245836a9-6b44-4639-9209-e7ad9035e293
# ╟─7752316d-9dd0-4403-aa08-22c977ff3727
# ╟─76624080-150a-4783-b675-794365dcecee
# ╠═cccf3dfb-8b3c-45e8-bb1c-e9579afc7e1a
# ╠═0f8e955e-fe9e-4b9d-ac24-1d2fb403a46c
# ╠═85bac8de-3224-4f02-a625-fdf55ec611a1
# ╠═0e485727-495c-444d-9fb4-f20bdaac2676
# ╠═6d94f549-645d-4c27-9e4f-046542b5fb16
# ╠═109a382d-8d41-4bc3-a23b-439a987b17c7
# ╠═a6687107-7448-451e-a3cf-04a3d2c3d7a5
# ╠═a391cd0a-f752-4efd-92de-43e7cec656d4
# ╠═dd4556a3-d6af-4aec-a59d-d6fcd4f4144d
# ╠═c728c5d5-7aa8-437d-a6d8-8dc75974c86e
# ╠═9a7e5903-69be-4e0a-8514-3e05feedfed5
# ╟─418525b7-c358-41da-b865-5df3feb15855
# ╠═a95e371e-9319-4c7e-b5d9-4c4a50d12cd7
# ╟─6ebe0cb0-ba35-411c-9a7a-a8b6eecf326f
# ╟─ef43e50a-5af8-4733-88a4-cd159d173034
# ╠═9dabad13-cfa4-4e06-950d-f7c7d96c1147
# ╟─e040094c-7511-4831-b94a-1c1185868202
# ╠═23ee0e85-a84b-4b63-b432-5526559efcee
# ╟─078c01f7-e47e-4af0-be1c-ac4527b735fd
# ╠═8b00d2b3-9182-42ab-8393-91707b813f60
# ╟─7899f488-9c48-466f-857d-f5a31b5820ab
# ╟─96f26378-846c-4964-935c-0372e2e86e91
# ╠═661feb84-339c-4bbe-a8a5-65de74ed58c8
# ╠═1e8a535e-25ea-490b-b545-e532c4fbc0f3
# ╟─a0849611-23b3-4a91-a054-f390bc6c9f0a
# ╠═33f889d7-e875-40d8-9d6d-cc87b0fbaf22
# ╠═0e4b3c36-2f09-405e-912b-22c893cd1715
# ╠═ddd8f749-3126-4563-8177-4941b6b0447b
# ╟─710ddf42-2397-4d96-9a61-ff5c600ccd43
# ╠═46084a31-e591-42f2-b8e6-02183ddfc6ac
# ╠═f21dc58e-d4e8-4314-b5dd-abbcb29efe86
# ╠═14b713b9-70c6-4506-a7b3-c67d021f8fce
# ╠═86b34265-6f2d-48d7-95ea-d10c8ae29ea9
# ╠═6bbb7a2c-464f-4e92-ab15-66745bac03ef
# ╠═bcf21d4b-cc76-4b9e-8cdb-7a17eb2af605
# ╟─e379461f-8896-4e2a-a71b-1871a8a37eb5
# ╠═7b7baa41-0185-4ed6-8fae-3a44e9912016
# ╟─f2f236f4-59f3-4c05-811d-078cd04ddd79
# ╠═c6a263eb-cb45-4ee7-9c02-549c89298652
# ╠═05ed4187-a01a-4a16-a0e7-b3867d252578
# ╟─6f7d4335-d9a6-4896-9d69-bfc1c2c1c3d0
# ╠═66815e8e-09d9-4b43-9f45-9379b3d34f78
# ╟─5e79d8e1-429c-414f-b3a6-8cf4b93d1336
# ╠═0b75c073-f167-4553-b746-539a14cfcf25
# ╠═4ed31219-9ce0-4f1b-8152-0002e64649ad
# ╠═0763f1f3-dfee-4d1f-a934-bf387f9c80ff
# ╠═c57c808a-297c-4887-bf20-5ad0207d055e
# ╠═d3307918-1fdb-4f87-bb92-67330d22e58b
# ╟─d25cda2f-6ec6-4c93-8860-f5ce9c3ee629
# ╠═444f6d74-273e-486d-905a-1443ec0e98df
# ╠═8cfdc784-4060-48b8-8d1a-3b8d11f7a9a7
# ╠═6c010734-e8a0-4000-88eb-f2a85d25ed99
# ╠═b7e27b45-78d4-41b5-9770-9632057413c6
# ╟─a1a10e2f-1b78-4b93-9295-7c0055e32692
# ╠═58eff13c-44b5-4f19-8a42-cf9907ac9515
# ╠═68c9d88a-99b7-49be-9ac4-1e06c694c1a6
# ╟─8a21fa0f-d3c3-4aa2-8b8b-74001d921c4a
# ╠═c5754b4b-b576-4257-95d2-8888bbd063ec
# ╠═0bc7df52-4ac9-42ac-9094-ecaf3c27da31
# ╠═8f5b8859-6b8c-4f2a-af3a-b13c2d33fe2a
# ╠═8082559e-a5b0-41a8-b8ed-aec3b09e5b2b
# ╟─b04ad0dc-10b5-433e-abc3-e87b4aa4f7eb
# ╠═2a973d9c-8b33-4c67-8177-73fd826c8dac
# ╠═c31a2d3f-902b-4be9-a64c-b04cb83ffaa4
# ╠═d03cf081-8ebe-4f5b-a81f-abc843e1bb65
# ╟─2ee1ca40-141f-40ad-b4c1-a2e025f69f95
# ╠═c2d877b5-d309-4868-925d-dab8d7d23403
# ╟─c239deed-8291-45aa-95cf-94df26e0136d
# ╠═7979b889-4782-45be-9a4f-91375f22f26f
# ╠═73d702b9-cdc0-4ce9-802d-89443c8412ab
# ╠═5bb0b72a-8c77-4fcb-bbde-d144986d9c1e
# ╠═d49936ff-c8c0-4a8d-a804-0fb56908b383
# ╠═ded5b462-06dd-43a4-93b0-c52ad87174eb
# ╠═868189ef-6a3f-425f-94ac-dbb0e1847b2e
# ╟─86b56683-c80e-4c0f-8b03-a4869860d04f
# ╠═323f3fd7-e9a9-4598-ad2e-c1790cf4a264
# ╠═2ab35999-3615-4f5c-8d89-36d77802fe9b
# ╠═2a01b228-f281-46c4-9764-fac6cc1b4217
# ╠═18b5c6d2-2230-4881-8066-51eff42125ae
# ╟─a5ae695b-bfc0-4425-9b64-bbeeba7da015
# ╠═eaf470e9-2898-41d5-a6d5-4cd846e9c0de
# ╠═a3ba0c9d-5f81-4023-9ce0-ff29536aa968
# ╟─67c46219-2183-47fd-bd3c-82facff98d53
# ╠═38a03e22-d596-4227-a9de-3ef54dc7256e
# ╠═2148cbb2-b41f-4df3-ab5c-55a89eff7bf1
# ╟─968de6ad-eb48-4d70-b431-209e609904aa
# ╠═63977532-9afa-454c-9f51-af6f4b238120
# ╠═aa9c9d45-bb7c-4eee-af87-6fbc01df271d
# ╠═bb0a7df4-7e84-472a-ab00-e3dd801daf8e
# ╠═3f640581-edcc-4c7a-86ba-b168f31fe4a3
# ╠═3bb65b71-d191-498b-81bf-40ffff4df1f4
# ╠═7c8608ba-759a-4bbf-be14-32813dcbf79b
# ╠═2fc78ec1-bf53-49da-b31a-6b5bf165eb81
# ╠═34f621c0-207e-41e5-8193-ad81d2d21a01
# ╠═29518390-5d57-4f2d-b617-d7699468caf9
# ╠═e041f6b5-92b1-46f0-b383-606625d59a4a
# ╠═b4fbbef8-3389-4fe9-9ebf-7b356bedd705
# ╠═22adb08e-4d8c-40c5-97ef-cb1f3f0f6d90
# ╠═a515407b-f749-48f2-b8ea-62940a186cce
# ╠═692a7f74-ca65-494b-b683-2d30e34e4c1e
# ╟─9533c662-80af-4dd4-bf25-02e894867360
# ╠═8b6d766a-8f7b-4b9a-9a15-0f7375087120
# ╠═807222ba-5ff8-4f33-a9a0-7c69b1dccf52
# ╠═16158266-36ed-44c3-a418-0c454955ce78
# ╟─580de17a-625d-420e-974c-86766197025e
# ╠═cb59f55b-c748-4a94-b344-e50a8fa7c690
# ╠═89cced40-f24e-499e-8bfd-19c3964f689b
# ╠═b9515b3a-b254-49ae-8c2c-b8ce7ced4d3a
# ╟─a8861082-2214-45f1-bc49-733efe74c949
# ╠═b59fa654-6946-4687-b14b-c2ef1f766f5c
# ╠═b12963ae-bf7d-4ef7-b1a8-e2d1e24f9b4b
# ╠═cfbe753d-85a8-445f-9eda-14a376d7e0c6
# ╠═54e9eda2-d564-453a-8ea8-4c8395be9ed6
# ╟─c53edeef-324a-418f-907d-aaf557cb8d24
# ╠═4831a3be-35d3-420c-8463-bb14a597cc6a
# ╠═08cee0a9-d358-4954-8e46-74de23d48d86
# ╠═f3f886d6-3010-4dd9-b42a-d5309463beb6
# ╠═5003d6c6-fa30-423e-80c1-a0f82e4085b9
# ╠═36b1822e-fe08-494b-a57d-5888163a7b54
# ╠═29ccc89b-f76a-44fa-8a89-e3ca10742ba1
# ╠═3385b22e-85ef-4bb0-8b1b-d03411c89b4f
# ╠═29d3cc8f-780b-449a-87ba-8d543ad2473b
# ╠═0a48e016-2fba-47cc-a212-47b4a3324b20
# ╠═57328e4f-e945-46eb-aba1-8b75dfeb4575
# ╠═8e785e48-f7ea-4d27-8977-f1879c8bd74b
# ╠═5feb46c0-3888-4586-8b12-f990d4d38912
# ╟─b23dc763-d91f-4d66-94d2-dcf96cb07f54
# ╠═8897acea-5efb-47a6-83a2-0c70fccfdb46
# ╠═798d8d16-1c19-400d-8a94-e08c7f991e33
# ╟─da44647a-36e4-4116-9698-df1cb059c2b7
# ╠═fb3ece76-f85c-41e1-a332-12c71d9d3cc0
# ╠═743e74ec-9c66-4665-845d-75ede418616b
# ╠═1aca6b92-7754-4cb3-b9e8-5d486e3bfcf8
# ╠═57b18cdd-27a8-44df-b959-f5e7c7eb7413
# ╠═3c9a219f-74ef-45fb-83e7-c497e0bee362
# ╠═e1264f57-f675-4f37-b4db-313cfc52ab8e
# ╠═7127fc35-a0af-4463-9448-a948f229fd47
# ╠═40157899-dffb-4e3a-b5ca-be3c23a465ae
# ╠═0c624592-df44-413e-9101-eefc3913d658
# ╠═3d0c3999-77a2-40d6-923b-78d3329e2154
# ╟─bd95428d-1077-4417-bfca-0c5da7378af2
# ╠═65d81268-9ff2-4a18-b0ce-4b105740dc8b
# ╠═8c1d1401-bc6b-4be3-8481-1c9a8f86f63d
# ╠═5b9d558a-2991-489e-be58-f5a5db0479f8
# ╠═6d4b0c74-4228-41e4-a8d0-98e0d71333b9
# ╠═67b3c66c-b3eb-438f-96e4-c09e117cde87
