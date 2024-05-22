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

# ╔═╡ faf59350-8d67-11ee-0bdd-2510e986118b
begin
    import Pkg; Pkg.activate()
    using CSV, Interpolations, DataFrames, CairoMakie, DifferentialEquations, Turing, StatsBase, PlutoUI, Distributions, Optim, Dierckx
end

# ╔═╡ 4391f124-cbef-46e5-8462-e4e5126f5b38
begin
	# cool plot theme
	import AlgebraOfGraphics
	
	AlgebraOfGraphics.set_aog_theme!(
		fonts=[
			AlgebraOfGraphics.firasans("Light"), 
			AlgebraOfGraphics.firasans("Light")
		]
	)
	
	resolution = (0.9 * 500, 0.9 * 380)
	
	update_theme!(
		fontsize=20, 
		linewidth=4,
		markersize=14,
		titlefont=AlgebraOfGraphics.firasans("Light"),
		size=resolution
	)

	colors = Dict(zip(
		["data", "model", "distn", "other"], 
		AlgebraOfGraphics.wongcolors()[1:4])
	)
end

# ╔═╡ 245836a9-6b44-4639-9209-e7ad9035e293
TableOfContents()

# ╔═╡ 7752316d-9dd0-4403-aa08-22c977ff3727
md"""
# Tank Geometry and Hole Measurements

we approximate the area from a helicopter view as a rounded rectangle.

see [here](https://mathworld.wolfram.com/RoundedRectangle.html).

see notes in `tank_geometry/cory_measurements.txt`.

"""

# ╔═╡ c976bc08-97b2-45c0-b1ff-1819e7290a68
struct TankMeasurements
	# top, bottom areas
	A_t::Float64
	A_b::Float64
	
	# height (cm)
	H::Float64

	# radius of hole
	r_hole::Float64

	# height of hole
    h_hole::Float64
end

# ╔═╡ 20fa4266-be80-4d8e-b1d0-155a40a1241f
begin
	# height of tank along slant [cm]
	H★ = 28.6

	# measurements of top cross-section
	L_t = 14.6 #cm
	W_t = 9.0 # cm
	p_t = 44.3 # perimeter, cm
	
	# measurements of bottom cross-section
	L_b = 13.4 # cm
	W_b = 7.8 # cm
	p_b = 40.1 # perimeter, cm
	
	# height of tank (from perpendicular) [cm]
	local δ = (L_t - L_b) / 2 # overhang
	H = sqrt(H★ ^ 2 - δ ^ 2)

	# solve for the r consistent with a rounded rectangle [cm]
	#  https://mathworld.wolfram.com/RoundedRectangle.html
	local r(p, L, W) = (p / 2 - (L + W)) / (π - 4)
	r_b = r(p_b, L_b, W_b)
	r_t = r(p_b, L_b, W_b)
	@show r_b, r_t

	# finally, calculate areas.
	local A(L, W, r) = (L - 2 * r) * (W - 2 * r) + 2 * r * (L + W - 4 * r) + π * r ^ 2
	A_b = A(L_b, W_b, r_b)
	A_t = A(L_t, W_t, r_t)
	@show A_b, A_t

	# hole in tank
	r_hole = 5 / 64 * 2.54 / 2 # cm
    a_hole = π * r_hole ^ 2 # cm²

    h_hole = 0.9 # cm

	tank_measurements = TankMeasurements(
		A_t, A_b, H, r_hole, h_hole
	)
end

# ╔═╡ 48d7273e-a48b-49fd-991b-6e29f64a0760
"""
cross-sectional area of water from helicopter view, as a function of liquid level, h.
"""
function A_of_h(h::Float64, tm::TankMeasurements)
	if h < 0 || h > tm.H
		error("tank over/under-flow!")
	end
	return h / tm.H * tm.A_t + (1 - h / tm.H) * tm.A_b
end

# ╔═╡ 9a7e5903-69be-4e0a-8514-3e05feedfed5
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="water level, h [cm]", 
		ylabel="area of tank\nfrom top-view, A(h) [cm²]"
	)
	lines!(range(0, H), 
		[A_of_h(hᵢ, tank_measurements) for hᵢ in range(0, tank_measurements.H)]
	)
	ylims!(0, nothing)
	fig
end

# ╔═╡ 418525b7-c358-41da-b865-5df3feb15855
md"
# Calibration of Liquid Level Sensor

read in data characterizing the calibration curve of the liquid level sensor.
note liquid level here is actually slanted h★.
"

# ╔═╡ a95e371e-9319-4c7e-b5d9-4c4a50d12cd7
begin
	calibration_data = CSV.read("calibration_curve.csv", DataFrame)
	sort!(calibration_data, "level sensor reading")
	rename!(calibration_data, "h [cm]" => "h★ [cm]")
end

# ╔═╡ 6d48d04b-0786-4d15-a07b-8941805a5b09
md"we correct for the slant, though slight.

this function maps the slanted h★ to the height h perpendicular to the ground.
"

# ╔═╡ 9af216b7-4bf2-42fb-bd95-5b2040d019a7
h★_to_h(h★) = H * h★ / H★

# ╔═╡ d3eeccff-39b3-429d-aec6-e7f1a500b729
h★_to_h(1.0) # note, this correction is TINY!

# ╔═╡ 6ebe0cb0-ba35-411c-9a7a-a8b6eecf326f
md"compute the true liquid level."

# ╔═╡ 385442da-f101-4ddf-8293-46d71d6a48fc
begin
	calibration_data[:, "h [cm]"] = h★_to_h.(calibration_data[:, "h★ [cm]"])
	calibration_data
end

# ╔═╡ 9dabad13-cfa4-4e06-950d-f7c7d96c1147
begin
	# map liquid level sensor reading to liquid level h.
	level_sensor_to_h = Spline1D(
		calibration_data[:, "level sensor reading"],
		calibration_data[:, "h [cm]"];
		k=1, s=2, bc="error"
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
	# xlims!(700, 800)
	# ylims!(20, 30)
	lines!(lsr, level_sensor_to_h.(lsr))
	fig
end

# ╔═╡ 078c01f7-e47e-4af0-be1c-ac4527b735fd
md"""
# Data Preprocessing
"""

# ╔═╡ 8b00d2b3-9182-42ab-8393-91707b813f60
function read_h_time_series(file::String)
	#=
	read in file
	=#
	data = CSV.read(file, DataFrame)
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
## Load in data
"""

# ╔═╡ 96f26378-846c-4964-935c-0372e2e86e91
md"time series data for two experiments without any object in the tank"

# ╔═╡ 759852ee-50e7-4deb-ac7e-c4693103c2a7
begin
	train_experiment = "no_obs_4_18_2.csv"
	test_experiment  = "no_obs_4_18_3.csv"
end

# ╔═╡ 661feb84-339c-4bbe-a8a5-65de74ed58c8
train_data = read_h_time_series(train_experiment)

# ╔═╡ 1e8a535e-25ea-490b-b545-e532c4fbc0f3
test_data = read_h_time_series(test_experiment)

# ╔═╡ 8d4d88da-a48a-4018-b6c6-3a2f041fcfb4
colors

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
	lines!(
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
viz_data(train_data)

# ╔═╡ ddd8f749-3126-4563-8177-4941b6b0447b
viz_data(test_data)

# ╔═╡ e379461f-8896-4e2a-a71b-1871a8a37eb5
md"""
# Dynamic model of the tank

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
function f(h, params, t)
	if h < params.h_hole
		return 0.0
	end
	return - π * params.r_hole ^ 2 * params.c * sqrt(
		2 * g * (h .- params.h_hole)) / params.A_of_h(h)
end

# ╔═╡ f21dc58e-d4e8-4314-b5dd-abbcb29efe86
function downsample(data::DataFrame, n::Int)
	# filter out sensor data that is below 508
	# this is the liquid level where the jet shooting outside the tank
	# dies, and it just runs down the tank (Toricelli's law invalid)
	data = filter("h [cm]" => (x -> x >= level_sensor_to_h(508)), data)
	
	ids = trunc.(Int, collect(range(1, nrow(data), length=n)))
	return data[ids, :]
end

# ╔═╡ 6f7d4335-d9a6-4896-9d69-bfc1c2c1c3d0
md"""
## minimize dynamic model loss
"""

# ╔═╡ 46084a31-e591-42f2-b8e6-02183ddfc6ac
@bind n_data_sample Select([20, 50, 100], default=20)

# ╔═╡ 66815e8e-09d9-4b43-9f45-9379b3d34f78
function loss(data::DataFrame, c::Float64, tm::TankMeasurements)	
	params = (
		# area of the hole
		r_hole = tm.r_hole,
		# fudge factor
		c = c, 
		# height of the hole
		h_hole = tm.h_hole,
		# area as a function of h
		A_of_h = h -> A_of_h(h, tm)
	)

	data = downsample(data, n_data_sample)
	h₀ = data[1, "h [cm]"]
	tspan = data[end, "t [s]"]
	
	prob = ODEProblem(f, h₀, tspan, params)
	h_of_t = solve(prob)
	
	cost = 0.0
	for i in 1:nrow(data)
		tᵢ = data[i, "t [s]"]
		cost += (data[i, "h [cm]"] - h_of_t(tᵢ)) ^ 2
	end
	
	return cost
end

# ╔═╡ 0b75c073-f167-4553-b746-539a14cfcf25
loss(train_data, 0.4, tank_measurements)

# ╔═╡ b0657fe7-29d8-4887-bb4c-cf72cfd434d4
md"""
## compute mle
"""

# ╔═╡ 4ed31219-9ce0-4f1b-8152-0002e64649ad
function compute_mle(data::DataFrame, tm::TankMeasurements)
	res = optimize(c -> loss(data, c, tm), 0.2, 0.9)
	return Optim.minimizer(res)
end

# ╔═╡ 0763f1f3-dfee-4d1f-a934-bf387f9c80ff
c_opt = compute_mle(train_data, tank_measurements)

# ╔═╡ c57c808a-297c-4887-bf20-5ad0207d055e
params = (
	# area of the hole
	r_hole = tank_measurements.r_hole,
	# fudge factor
	c = c_opt,
	# height of the hole
	h_hole = tank_measurements.h_hole,
	# area as a function of h
	A_of_h = h -> A_of_h(h, tank_measurements)
)

# ╔═╡ d3307918-1fdb-4f87-bb92-67330d22e58b
begin
	train_h₀ = train_data[1, "h [cm]"]
	prob = ODEProblem(f, train_h₀, 1.1 * train_data[end, "t [s]"], params)
	sol = solve(prob, saveat=1.0)
	sim_data = DataFrame(sol)
end

# ╔═╡ d25cda2f-6ec6-4c93-8860-f5ce9c3ee629
md"""
## visualize dynamic model
"""

# ╔═╡ 444f6d74-273e-486d-905a-1443ec0e98df
function viz_sim_fit(data::DataFrame, sim_data::DataFrame; 
		         savename::Union{Nothing, String}=nothing)
	 fig = Figure()
	 ax = Axis(
		fig[1, 1], 
		xlabel="time, t [s]", 
		ylabel="liquid level, h[cm]"
	)
	
	ds_data = downsample(data, n_data_sample)
	
	scatter!(
		ds_data[:, "t [s]"], 
		ds_data[:, "h [cm]"], 
		label="experiment"
	)
	
	lines!(sim_data[:, "timestamp"], sim_data[:, "value"], 
		label="model", color=Cycled(2)
	)

	axislegend()
	if ! isnothing(savename)
		save("$savename.pdf", fig)
	end
	fig
end

# ╔═╡ 8cfdc784-4060-48b8-8d1a-3b8d11f7a9a7
viz_sim_fit(train_data, sim_data)

# ╔═╡ a1a10e2f-1b78-4b93-9295-7c0055e32692
md"""
# Bayesian inference for model parameters
"""

# ╔═╡ 58eff13c-44b5-4f19-8a42-cf9907ac9515
@bind n_MC_sample Select([10, 50, 100, 250], default=10)

# ╔═╡ 8a21fa0f-d3c3-4aa2-8b8b-74001d921c4a
md"""
## infer model parameters for object-free experiment (training)
"""

# ╔═╡ 8f5b8859-6b8c-4f2a-af3a-b13c2d33fe2a
@model function infer_params(data::DataFrame, tm::TankMeasurements)
	#=
	prior distributions
	=#
	# defines variance for measuring length.
	σ_ℓ ~ Uniform(0.0, 0.5) # cm
	
	# bottom, top tank area measurements
	# std of product of two Guassians
	#   https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
	A_b ~ Normal(tm.A_b, σ_ℓ ^ 2 / 2) # cm²
	A_t ~ Normal(tm.A_t, σ_ℓ ^ 2 / 2) # cm²

	# height of tank
	H ~ Truncated(
		Normal(tm.H, σ_ℓ),
		tm.H - 0.25, tm.H + 0.25
	) # cm

	# drill precision
	dp ~ Uniform(0.0, 0.02) # cm
	# radius of the hole
	r_hole ~ Truncated(
		Normal(tm.r_hole, dp),
		tm.r_hole - 2 * dp, 
		tm.r_hole + 2 * dp
	) # cm

	# discharge coefficient. Wikipedia says 0.65 for water.
	dcv ~ Uniform(0.0, 0.15) # discharge coeff variance
	c ~ Truncated(
		Normal(0.65, dcv), 
		0.3, 1.0
	) # unitless

	# height of the hole
	h_hole ~ Truncated(
		Normal(tm.h_hole, σ_ℓ),
		tm.h_hole - σ_ℓ, 
		tm.h_hole + σ_ℓ
	) # cm
	# TODO: think if this makes sense. hard truncations?
	
	# defines variance of liquid level sensor
	#   (Treated as an unknown and inferred)
	σ ~ Uniform(0.0, 1.0) # cm

	# initial liquid level
	h₀_obs = data[1, "h [cm]"] # cm
	h₀ ~ Normal(h₀_obs, σ)

	#=
	set up dynamic model for h(t)
	=#
	function my_A_of_h(h)
		if h < 0.0
			error("h < 0")
		elseif h > H
			error("h > H")
		end
		return h / H * A_t + (1 - h / H) * A_b
	end
	
	# parameter for ODE solver
	params = (
			  r_hole=r_hole,
			  c=c,
			  h_hole=h_hole,
			  A_of_h=my_A_of_h 
			)
	
	# set up and solve ODE
	tspan = (0.0, 1.1 * data[end, "t [s]"])
	prob = ODEProblem(f, h₀, tspan, params, saveat=1.0)
	h_of_t = solve(prob)
	
	#=
	code up likelihood
	=#
	for i in 2:nrow(data) # start at 2 b/c IC handled with informative prior
		tᵢ = data[i, "t [s]"]
		# continuity=:right
		data[i, "h [cm]"] ~ Normal(h_of_t(tᵢ, continuity=:right)[1], σ)
	end

	return nothing
end

# ╔═╡ 8082559e-a5b0-41a8-b8ed-aec3b09e5b2b
begin
	ds_train_data = downsample(train_data, n_data_sample)
	train_model = infer_params(ds_train_data, tank_measurements)
	
	train_chain = sample(
		train_model, 
		NUTS(0.65), MCMCSerial(), n_MC_sample, 3; progress=true
	)
end

# ╔═╡ 7ebe4680-c583-4f92-8bae-dd84c3fb5139
train_posterior = DataFrame(train_chain)

# ╔═╡ c239deed-8291-45aa-95cf-94df26e0136d
md"""
## visualize posterior distribution of model parameters
"""

# ╔═╡ ccb1f005-567d-47f8-bec1-8db268d878ec
inferred_params = ["σ_ℓ", "A_b", "A_t", "H", "r_hole", "c", "h_hole", "σ", "h₀", 
					"dp"]

# ╔═╡ 5bb0b72a-8c77-4fcb-bbde-d144986d9c1e
function viz_posterior(posterior::DataFrame, params::Vector{String},
			           tm::TankMeasurements, h₀_obs::Float64
)
	params_to_title = Dict(
						"A_b" => "area of the tank bottom", 
						"A_t" => "area of the tank top",
						"r_hole" => "radius of the oriface",
						"c" => "discharge coefficient",
						"h_hole" => "height of orifice from the base of tank",
						"h₀" => "initial water level",
						"σ_ℓ" => "variance for measuring length",
						"H" => "height of tank", 
						"σ" => "variance of liquid level sensor", 
						"dcv" => "variance of discharge coefficient",
						"dp" => "drill precision"
	)
	
	fig = Figure(resolution=(1350, 600))
	j, i = 1, 1
	 
	for p in params
		# vizualize the distribution
		ax = Axis(fig[i, j])
		hist!(ax, posterior[:, p], color=Cycled(1))

		# plot equal-tailed 80% interval
		lo, hi = quantile(posterior[:, p], [0.1, 0.9])
		lines!(ax, [lo, hi], [0, 0], linewidth=5, color=Cycled(2))
		
		ax.xlabel = params_to_title[p]

		if j == 1
			ax.ylabel = "density"
		end
		
		if ! (p in ["h₀", "c", "σ_ℓ", "σ", "dp"])
			p_obs = getfield(tm, Symbol(p))
			# plot measured value
			vlines!(ax, p_obs, linestyle=:dash, color=Cycled(3), 
					label="measurement")
		elseif p == "h₀"
			vlines!(ax, h₀_obs, linestyle=:dash, color=Cycled(3), label="true value")
		end
		
		i += 1
		if i > 2
			j += 1
			i = 1
		end
	end
	
	return fig	
end

# ╔═╡ ded5b462-06dd-43a4-93b0-c52ad87174eb
viz_posterior(train_posterior, inferred_params, tank_measurements, train_h₀)

# ╔═╡ 86b56683-c80e-4c0f-8b03-a4869860d04f
md"## posterior predictive check"

# ╔═╡ 2ab35999-3615-4f5c-8d89-36d77802fe9b
function viz_fit(posterior::DataFrame, data::DataFrame; 
				savename::Union{String, Nothing}=nothing)
	fig = Figure(resolution=(700, 500))
	ax = Axis(
		fig[1, 1], 
		xlabel="time, t [s]", 
		ylabel="liquid level, h [cm]"
	)
	
	ts = range(0, maximum(data[:, "t [s]"]), length=500)
	tspan = (0.0,  maximum(ts) * 1.05)
	

	# sample posterior models
	for i in 1:nrow(posterior)
		params = (
			  r_hole=posterior[i, "r_hole"],
			  c=posterior[i, "c"],
			  h_hole=posterior[i, "h_hole"],
			  A_of_h=h -> h / posterior[i, "H"] * posterior[i, "A_t"] + 
			  	    (1 - h / posterior[i, "H"]) * posterior[i, "A_b"]
			)

		# set up, solve ODE
		prob = ODEProblem(f, posterior[i, "h₀"], tspan, params)
		sim_data = DataFrame(solve(prob, saveat=1.0))
			
		
		lines!(sim_data[:, "timestamp"], sim_data[:, "value"], 
			label="model", color=(colors["model"], 0.1))
	end	
	scatter!(
		data[:, "t [s]"], 
		data[:, "h [cm]"],
		label="experiment",
		color=colors["data"]
		)
	axislegend(unique=true)
	ylims!(0, nothing)
	xlims!(0, nothing)
	if savename!=nothing
		save( "$savename.pdf", fig)
	end
	return fig
end

# ╔═╡ 2a01b228-f281-46c4-9764-fac6cc1b4217
viz_fit(train_posterior, ds_train_data)

# ╔═╡ a5ae695b-bfc0-4425-9b64-bbeeba7da015
md"""
## validate posterior with test data set
"""

# ╔═╡ eaf470e9-2898-41d5-a6d5-4cd846e9c0de
function viz_test(posterior::DataFrame, test_data::DataFrame;
				savename::Union{String, Nothing}=nothing
)
	fig = Figure(resolution=(700, 700))
	
	ax_stopping = Axis(
					fig[1, 1], 
					ylabel="# samples", 
					xticks=([],[]),
					height=100
	)
	ax = Axis(
		fig[2, 1], 
		xlabel="time, t [s]", 
		ylabel="liquid level, h [cm]",
	)

	# sample from the train posterior
	emptying_time = zeros(nrow(posterior))
	for i in 1:nrow(posterior)
		## Check for first instance when the liquid level is the same as the height of the oriface from the base
		condition(h, t, integrator) = h[1] - posterior[i, "h_hole"]
		# condition(h, t, integrator) = h[1] - posterior[i, "h_hole"]
		
		## Retrive the emptying time [t] when h(t) = h_hole
		function affect!(integrator)
			emptying_time[i] = integrator.t
		end
		# 
		cb = ContinuousCallback(condition, affect!)

		params = (
			  r_hole=posterior[i, "r_hole"],
			  c=posterior[i, "c"],
			  h_hole=posterior[i, "h_hole"],
			  A_of_h=h -> h / posterior[i, "H"] * posterior[i, "A_t"] + 
			  	    (1 - h / posterior[i, "H"]) * posterior[i, "A_b"]
		)

		# sample an initial condtion
		h₀ = test_data[1, "h [cm]"] + posterior[i, "σ"] * rand()
		tspan =  (0.0, 3.0 * test_data[end, "t [s]"])
		
		prob = ODEProblem(f, h₀, tspan, params)
		sim_data = DataFrame(solve(prob, callback=cb, saveat=1.0))
			
		# plot trajectories
		lines!(ax, sim_data[:, "timestamp"], sim_data[:, "value"], 
			label="model", color=(colors["model"], 0.1))
	end

	#Plot emptying time
	hist!(ax_stopping, emptying_time, color=Cycled(3))
	# True emptying time
	vlines!(ax_stopping, test_data[end, "t [s]"], linestyle=:dash, color=:black)

	scatter!(
		ax,
		test_data[:, "t [s]"], 
		test_data[:, "h [cm]"],
		label="experiment",
		color=colors["data"]
		)
	xlims!(0, maximum(emptying_time))
	ylims!(ax, 0, nothing)

	linkxaxes!(ax, ax_stopping)
	axislegend(ax, unique=true)
	
	if isnothing(savename)
		save("$savename.pdf", fig)
	end
	
	return fig 
end	

# ╔═╡ dd5ee69d-fb0a-484b-b1bd-32a589a9ed11
ds_test_data = downsample(test_data, n_data_sample)

# ╔═╡ a3ba0c9d-5f81-4023-9ce0-ff29536aa968
viz_test(train_posterior, ds_test_data, savename="test")

# ╔═╡ 9533c662-80af-4dd4-bf25-02e894867360
md"""
# Bayesian inference of object shape

## read experimental data
"""

# ╔═╡ b06a1c07-6250-4324-8802-010e5d847edb
begin
	obstruction_data = ["obs_4_18_2.csv", "obs_4_18_1.csv"]
	@bind obstruction Select(obstruction_data)
end

# ╔═╡ 8b6d766a-8f7b-4b9a-9a15-0f7375087120
block_data = read_h_time_series(obstruction)

# ╔═╡ 16158266-36ed-44c3-a418-0c454955ce78
begin
	function viz_obstruction_w_no_obstruction(obstruction_data, no_obstruction_data)
		fig = Figure()
		ax = Axis(
			fig[1, 1], 
			xlabel="time, t [s]", 
			ylabel="water level, h [cm]"
		)
		lines!(
			no_obstruction_data[:, "t [s]"], 
			no_obstruction_data[:, "h [cm]"],
			label="no object",
			color=colors["data"]
		)
		lines!(
			obstruction_data[:, "t [s]"], 
			obstruction_data[:, "h [cm]"],
			label="with object",
			color=colors["other"]
		)
		axislegend()
		ylims!(0, nothing)
		xlims!(0, nothing)
		save("h_of_t_with_without_object.pdf", fig)
		fig
	end
	
	viz_obstruction_w_no_obstruction(block_data, train_data)
end

# ╔═╡ 580de17a-625d-420e-974c-86766197025e
md"## measured area"

# ╔═╡ cb59f55b-c748-4a94-b344-e50a8fa7c690
begin
	object_true_area = CSV.read("obstacle_area.csv", DataFrame)
	rename!(object_true_area, "area " => "area [cm²]")
	select!(object_true_area, ["h [cm]", "area [cm²]"])
end

# ╔═╡ 89cced40-f24e-499e-8bfd-19c3964f689b
A_of_object = Spline1D(object_true_area[:, "h [cm]"], 
					  object_true_area[:, "area [cm²]"]; k=3, s=5, bc="zero")

# ╔═╡ b9515b3a-b254-49ae-8c2c-b8ce7ced4d3a
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="water level, h [cm]", 
		ylabel="area of object\nfrom top-view, A(h) [cm²]"
	)
	
	h_range = range(0.0, object_true_area[end, "h [cm]"])
	
	lines!(h_range, A_of_object.(h_range))
	scatter!(object_true_area[:, "h [cm]"], object_true_area[:, "area [cm²]"])
	ylims!(0, nothing)
	fig
end

# ╔═╡ c56a1461-d359-4aec-9564-b1abfcee8b6b
block_data

# ╔═╡ a8861082-2214-45f1-bc49-733efe74c949
md"## simulate model with knowledge of true area"

# ╔═╡ 9c2eb6c9-fcd6-49d4-bcc6-dd0c774261b5
A(h) = A_of_h(h, tank_measurements) - A_of_object(h)

# ╔═╡ b59fa654-6946-4687-b14b-c2ef1f766f5c
object_params = (
		# area of the hole
		r_hole = tank_measurements.r_hole,
		# fudge factor
		c = c_opt,
		# height of the hole
		h_hole = tank_measurements.h_hole,
		# area as a function of h"
		A_of_h = A
	)

# ╔═╡ b12963ae-bf7d-4ef7-b1a8-e2d1e24f9b4b
begin
	local h₀ = block_data[1, "h [cm]"]
	local prob = ODEProblem(f, h₀, block_data[end, "t [s]"], object_params)
	local sol = solve(prob, saveat=1.0)
	obs_sim_data = DataFrame(sol)
end

# ╔═╡ cfbe753d-85a8-445f-9eda-14a376d7e0c6
viz_sim_fit(block_data, obs_sim_data)

# ╔═╡ 798d8d16-1c19-400d-8a94-e08c7f991e33
@model function infer_object_area(
	data::DataFrame, train_posterior::DataFrame, γ::Float64; N=10
)
	#=
	yesterday's posterior is today's prior
	=#
	# variance for measuring length.
	σ_ℓ = mean(train_posterior.σ_ℓ)
	σ = mean(train_posterior.σ)
	
	#=
	prior distributions
	=#
	# defines variance for measuring length.
	# height of tank
	H ~ Normal(mean(train_posterior.H), std(train_posterior.H))

	# radius of the hole
	r_hole ~ Normal(mean(train_posterior.r_hole), std(train_posterior.r_hole))

	# discharge coefficient. 
	c ~ Normal(mean(train_posterior.c), std(train_posterior.c))

	# height of the hole
	h_hole ~ Normal(mean(train_posterior.h_hole), std(train_posterior.h_hole))

	# initial liquid level
	h₀_obs = data[1, "h [cm]"]
	h₀ ~ Normal(h₀_obs, σ)

	# bottom, top tank area measurements
	# std of product of two Guassians
	#   https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
	A_b ~ Normal(mean(train_posterior.A_b), σ_ℓ ^ 2 / 2)
	A_t ~ Normal(mean(train_posterior.A_t), σ_ℓ ^ 2 / 2)

	# Random distribution of obstacle area grid at points N 
	As ~ filldist(Normal(), N) # cm²

	function A_of_tank(h)
		if h < 0
			h = 0
		elseif h > H
			h = H
		end
		return h / H * A_t + (1 - h / H) * A_b
	end

	# corresponding h's for the unknown A's
	hs = range(h_hole - 1.0, h₀_obs + 1.0, length=N)
	# prior: object could be ANY size
	As[1] ~ Uniform(0.1, 0.9 * A_of_tank(hs[1]))
	for i in 2:N
		As[i] ~ Truncated(
						  As[i-1] + γ * Normal(0.0, 1.0), 
						  0.0, 
						  0.99 * A_of_tank(hs[i])
						 )
	end

	#=
	set up dynamic model for h(t)
	=#
	A_of_object = linear_interpolation(hs, As, extrapolation_bc=Throw())

	# parameter for ODE solver
	params = (
			  r_hole=r_hole,
			  c=c,
			  h_hole=h_hole,
			  A_of_h=h ->  A_of_tank(h) - A_of_object(h)
			)
	
	# set up ODE
	tspan = (0.0, 4000)
	prob = ODEProblem(f, h₀, tspan, params, saveat=1.0)
	h_of_t = solve(prob)
	
	# Observations.
	for i in 2:nrow(data)
		tᵢ = data[i, "t [s]"]
		data[i, "h [cm]"] ~ Normal(h_of_t(tᵢ, continuity=:right)[1], σ)
	end
	return nothing
end



# ╔═╡ daa1e6c7-867b-4cf7-b7b8-dc24f859ec96
ds_block_data = downsample(block_data, 10)

# ╔═╡ 02939a87-e811-4ae4-8b6b-173370029889
begin
	γ = 6.5
	object_tank_model = infer_object_area(ds_block_data, train_posterior, γ)
	object_posterior = DataFrame(
		sample(object_tank_model, NUTS(0.65), MCMCSerial(), n_MC_sample, 3; 
									  progress=true)
	)
end

# ╔═╡ a127225a-5b79-4074-a16b-cecd11030800
function viz_area(posterior::DataFrame, data::DataFrame; N=10, 
				  γ::Float64=γ)
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="height [cm]", ylabel="Area [cm²]")
	
	for j in 1:nrow(posterior)
		As = Vector{Float64}()
		[append!(As, posterior[j, "As[$(i)]"]) for i in 1:N]

		hs = range(posterior[j, "h_hole"], data[1, "h [cm]"], length=N) 
		lines!(hs, As, label="model", color=(:green, 0.1))
	end

	scatter!(obstacle_true_area[:, "h [cm]"], obstacle_true_area[:, "area [cm²]"], 
			label="true_obstacle_area") 

	axislegend(unique=true, position=:rb)

	ax.title = "γ=$γ; smooth_grid=$N"
	save("area_inferred_3.png", fig)
	return fig
end

# ╔═╡ b43f9f58-94fd-4c92-8e91-9a6b86cfc041
viz_area(area_posterior, test_infer)

# ╔═╡ 327ba163-4c4f-4702-ab5e-dcce96904fc5
h_of_t = Spline1D(block_data[:, "t [s]"], block_data[:, "h [cm]"]; 
		w=ones(nrow(block_data)), k=1, s=10)

# ╔═╡ 58f6b039-fec4-4ffd-b3e4-6b8145e05a8d
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="time, t [s]", 
		ylabel="water level, h [cm]"
	)
	local t = block_data[:, "t [s]"]
	scatter!(block_data[:, "t [s]"], block_data[:, "h [cm]"])
	lines!(t, h_of_t.(t), color=:red, linewidth=1)
	ylims!(0, nothing)
	fig
end

# ╔═╡ 0c8cfb2b-c836-4697-a452-405c0d17caff
begin
	ts = range(0.0, block_data[end, "t [s]"], length=100)
	dh_dt = derivative(h_of_t, ts)
end

# ╔═╡ aaa8f514-bf02-48ba-b384-353c9b58c794
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="water level, h [cm]", 
		ylabel="dh/dt [cm/s]"
	)
	# local t = block_data[:, "t [s]"]
	# scatter!(block_data[:, "t [s]"], block_data[:, "h [cm]"])
	lines!(block_data[:, "h [cm]"], derivative(h_of_t, block_data[:, "t [s]"]))
	# ylims!(0, nothing)
	fig
end

# ╔═╡ 6ee72380-9bcb-4c21-b979-f49d7fc98d1d
params

# ╔═╡ e204ab6d-2be5-47e8-acb1-4b661285b03f
begin
	local ts = range(0.0, block_data[end, "t [s]"], length=100)
	local hs = h_of_t.(ts)
	
	A_obs = - π * params.r_hole ^ 2 * params.c * 
			sqrt.(2 * g * (hs .- params.h_hole)) ./ derivative(h_of_t, ts)
end

# ╔═╡ 7374ac4b-c2ac-4d6c-9bca-737f71dfb4f0
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="water level, h [cm]", 
		ylabel="Area of obstruction [cm²]"
	)
	local ts = range(0.0, block_data[end, "t [s]"], length=100)
	local hs = h_of_t.(ts)
	
	lines!(hs, A_of_h.(hs) - A_obs)
	print(maximum(A_of_h.(hs)))
	# ylims!(0, nothing)
	fig
end

# ╔═╡ 9d957844-8e48-4f16-8a78-8332bf5f9e2b
minimum(A_obs)

# ╔═╡ Cell order:
# ╠═faf59350-8d67-11ee-0bdd-2510e986118b
# ╠═4391f124-cbef-46e5-8462-e4e5126f5b38
# ╠═245836a9-6b44-4639-9209-e7ad9035e293
# ╟─7752316d-9dd0-4403-aa08-22c977ff3727
# ╠═c976bc08-97b2-45c0-b1ff-1819e7290a68
# ╠═20fa4266-be80-4d8e-b1d0-155a40a1241f
# ╠═48d7273e-a48b-49fd-991b-6e29f64a0760
# ╠═9a7e5903-69be-4e0a-8514-3e05feedfed5
# ╟─418525b7-c358-41da-b865-5df3feb15855
# ╠═a95e371e-9319-4c7e-b5d9-4c4a50d12cd7
# ╟─6d48d04b-0786-4d15-a07b-8941805a5b09
# ╠═9af216b7-4bf2-42fb-bd95-5b2040d019a7
# ╠═d3eeccff-39b3-429d-aec6-e7f1a500b729
# ╟─6ebe0cb0-ba35-411c-9a7a-a8b6eecf326f
# ╠═385442da-f101-4ddf-8293-46d71d6a48fc
# ╠═9dabad13-cfa4-4e06-950d-f7c7d96c1147
# ╟─e040094c-7511-4831-b94a-1c1185868202
# ╠═23ee0e85-a84b-4b63-b432-5526559efcee
# ╟─078c01f7-e47e-4af0-be1c-ac4527b735fd
# ╠═8b00d2b3-9182-42ab-8393-91707b813f60
# ╟─7899f488-9c48-466f-857d-f5a31b5820ab
# ╟─96f26378-846c-4964-935c-0372e2e86e91
# ╠═759852ee-50e7-4deb-ac7e-c4693103c2a7
# ╠═661feb84-339c-4bbe-a8a5-65de74ed58c8
# ╠═1e8a535e-25ea-490b-b545-e532c4fbc0f3
# ╠═8d4d88da-a48a-4018-b6c6-3a2f041fcfb4
# ╟─a0849611-23b3-4a91-a054-f390bc6c9f0a
# ╠═33f889d7-e875-40d8-9d6d-cc87b0fbaf22
# ╠═0e4b3c36-2f09-405e-912b-22c893cd1715
# ╠═ddd8f749-3126-4563-8177-4941b6b0447b
# ╟─e379461f-8896-4e2a-a71b-1871a8a37eb5
# ╠═7b7baa41-0185-4ed6-8fae-3a44e9912016
# ╟─f2f236f4-59f3-4c05-811d-078cd04ddd79
# ╠═c6a263eb-cb45-4ee7-9c02-549c89298652
# ╠═f21dc58e-d4e8-4314-b5dd-abbcb29efe86
# ╟─6f7d4335-d9a6-4896-9d69-bfc1c2c1c3d0
# ╠═46084a31-e591-42f2-b8e6-02183ddfc6ac
# ╠═66815e8e-09d9-4b43-9f45-9379b3d34f78
# ╠═0b75c073-f167-4553-b746-539a14cfcf25
# ╟─b0657fe7-29d8-4887-bb4c-cf72cfd434d4
# ╠═4ed31219-9ce0-4f1b-8152-0002e64649ad
# ╠═0763f1f3-dfee-4d1f-a934-bf387f9c80ff
# ╠═c57c808a-297c-4887-bf20-5ad0207d055e
# ╠═d3307918-1fdb-4f87-bb92-67330d22e58b
# ╟─d25cda2f-6ec6-4c93-8860-f5ce9c3ee629
# ╠═444f6d74-273e-486d-905a-1443ec0e98df
# ╠═8cfdc784-4060-48b8-8d1a-3b8d11f7a9a7
# ╟─a1a10e2f-1b78-4b93-9295-7c0055e32692
# ╠═58eff13c-44b5-4f19-8a42-cf9907ac9515
# ╟─8a21fa0f-d3c3-4aa2-8b8b-74001d921c4a
# ╠═8f5b8859-6b8c-4f2a-af3a-b13c2d33fe2a
# ╠═8082559e-a5b0-41a8-b8ed-aec3b09e5b2b
# ╠═7ebe4680-c583-4f92-8bae-dd84c3fb5139
# ╟─c239deed-8291-45aa-95cf-94df26e0136d
# ╠═ccb1f005-567d-47f8-bec1-8db268d878ec
# ╠═5bb0b72a-8c77-4fcb-bbde-d144986d9c1e
# ╠═ded5b462-06dd-43a4-93b0-c52ad87174eb
# ╟─86b56683-c80e-4c0f-8b03-a4869860d04f
# ╠═2ab35999-3615-4f5c-8d89-36d77802fe9b
# ╠═2a01b228-f281-46c4-9764-fac6cc1b4217
# ╟─a5ae695b-bfc0-4425-9b64-bbeeba7da015
# ╠═eaf470e9-2898-41d5-a6d5-4cd846e9c0de
# ╠═dd5ee69d-fb0a-484b-b1bd-32a589a9ed11
# ╠═a3ba0c9d-5f81-4023-9ce0-ff29536aa968
# ╟─9533c662-80af-4dd4-bf25-02e894867360
# ╠═b06a1c07-6250-4324-8802-010e5d847edb
# ╠═8b6d766a-8f7b-4b9a-9a15-0f7375087120
# ╠═16158266-36ed-44c3-a418-0c454955ce78
# ╟─580de17a-625d-420e-974c-86766197025e
# ╠═cb59f55b-c748-4a94-b344-e50a8fa7c690
# ╠═89cced40-f24e-499e-8bfd-19c3964f689b
# ╠═b9515b3a-b254-49ae-8c2c-b8ce7ced4d3a
# ╠═c56a1461-d359-4aec-9564-b1abfcee8b6b
# ╟─a8861082-2214-45f1-bc49-733efe74c949
# ╠═9c2eb6c9-fcd6-49d4-bcc6-dd0c774261b5
# ╠═b59fa654-6946-4687-b14b-c2ef1f766f5c
# ╠═b12963ae-bf7d-4ef7-b1a8-e2d1e24f9b4b
# ╠═cfbe753d-85a8-445f-9eda-14a376d7e0c6
# ╠═798d8d16-1c19-400d-8a94-e08c7f991e33
# ╠═daa1e6c7-867b-4cf7-b7b8-dc24f859ec96
# ╠═02939a87-e811-4ae4-8b6b-173370029889
# ╠═a127225a-5b79-4074-a16b-cecd11030800
# ╠═b43f9f58-94fd-4c92-8e91-9a6b86cfc041
# ╠═327ba163-4c4f-4702-ab5e-dcce96904fc5
# ╠═58f6b039-fec4-4ffd-b3e4-6b8145e05a8d
# ╠═0c8cfb2b-c836-4697-a452-405c0d17caff
# ╠═aaa8f514-bf02-48ba-b384-353c9b58c794
# ╠═6ee72380-9bcb-4c21-b979-f49d7fc98d1d
# ╠═e204ab6d-2be5-47e8-acb1-4b661285b03f
# ╠═7374ac4b-c2ac-4d6c-9bca-737f71dfb4f0
# ╠═9d957844-8e48-4f16-8a78-8332bf5f9e2b
