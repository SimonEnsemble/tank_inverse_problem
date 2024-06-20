### A Pluto.jl notebook ###
# v0.19.42

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
    using CSV, Interpolations, DataFrames, CairoMakie, DifferentialEquations, Turing, StatsBase, PlutoUI, Distributions, Optim, Dierckx, MakieThemes
end

# ╔═╡ 4391f124-cbef-46e5-8462-e4e5126f5b38
begin
	# see https://makieorg.github.io/MakieThemes.jl/dev/themes/ggthemr
	gg_theme = :fresh 
	set_theme!(ggthemr(gg_theme))
	
	figsize = (0.9 * 500, 0.9 * 380)
	
	update_theme!(
		fontsize=20, 
		linewidth=4,
		markersize=14,
		# titlefont=AlgebraOfGraphics.firasans("Light"),
		size=figsize
	)

	theme_colors = MakieThemes.GGThemr.ColorTheme[gg_theme][:swatch]
	colors = Dict(zip(
		["data", "model", "distn", "other"], 
		theme_colors[1:4])
	)
end

# ╔═╡ 245836a9-6b44-4639-9209-e7ad9035e293
TableOfContents()

# ╔═╡ 7752316d-9dd0-4403-aa08-22c977ff3727
md"""
# tank geometry and measurements

we characterize the area of the liquid holding tank, as a function of height $h$, from a helicopter view. the cross-sectional area we model as a [rounded rectangle](https://mathworld.wolfram.com/RoundedRectangle.html). see notes in `tank_geometry/cory_measurements.txt`.
"""

# ╔═╡ c976bc08-97b2-45c0-b1ff-1819e7290a68
struct TankMeasurements
	# top, bottom areas (cm²)
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
	# height of tank (along slant) [cm]
	H★ = 28.6

	# measurements of top cross-section
	L_t = 14.6 #cm
	W_t = 9.0  # cm
	p_t = 44.3 # perimeter, cm
	
	# measurements of bottom cross-section
	L_b = 13.4 # cm
	W_b = 7.8  # cm
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
	my_A(L, W, r) = (L - 2 * r) * (W - 2 * r) + 2 * r * (L + W - 4 * r) + π * r ^ 2
	A_b = my_A(L_b, W_b, r_b)
	A_t = my_A(L_t, W_t, r_t)
	@show A_b, A_t

	# hole in tank (radius, area, height)
	r_hole = 5 / 64 * 2.54 / 2 # cm
    a_hole = π * r_hole ^ 2    # cm²
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
	# check for over/underflow
	h < 0.0  ? error("tank underflow!") : nothing
	h > tm.H ? error("tank overflow!")  : nothing
	# linearly interpolate top and bottom areas
	return h / tm.H * tm.A_t + (1 - h / tm.H) * tm.A_b
end

# ╔═╡ 9a7e5903-69be-4e0a-8514-3e05feedfed5
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="water height, h [cm]", 
		ylabel="cross-sectional area, A(h) [cm²]"
	)
	lines!(range(0, H), 
		[A_of_h(hᵢ, tank_measurements) for hᵢ in range(0, tank_measurements.H)]
	)
	vlines!(tank_measurements.H, linestyle=:dash, color="gray")
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
	calibration_data = CSV.read("calibration_curve.csv", DataFrame)
	sort!(calibration_data, "level sensor reading")
	rename!(calibration_data, "h [cm]" => "h★ [cm]")
end

# ╔═╡ 6d48d04b-0786-4d15-a07b-8941805a5b09
md"we correct for the slant, though slight, in the tank, which makes the level strip slanted. this function maps the slanted h★ to the height h perpendicular to the ground.
"

# ╔═╡ 9af216b7-4bf2-42fb-bd95-5b2040d019a7
h★_to_h(h★) = H * h★ / H★

# ╔═╡ 6ebe0cb0-ba35-411c-9a7a-a8b6eecf326f
md"compute the true liquid level."

# ╔═╡ 385442da-f101-4ddf-8293-46d71d6a48fc
begin
	calibration_data[:, "h [cm]"] = h★_to_h.(calibration_data[:, "h★ [cm]"])
	calibration_data
end

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
## read in train/test time series data
"""

# ╔═╡ 96f26378-846c-4964-935c-0372e2e86e91
md"time series data for two experiments without any object in the tank"

# ╔═╡ 2ddf387c-5a61-4490-9746-96e1589c7a74
train_experiment = "no_obs_4_18_2.csv"

# ╔═╡ b2228d5c-16b4-4fee-b9b6-1112d7cf391c
test_experiment  = "no_obs_4_18_3.csv"

# ╔═╡ 661feb84-339c-4bbe-a8a5-65de74ed58c8
all_train_data = read_h_time_series(train_experiment)

# ╔═╡ 1e8a535e-25ea-490b-b545-e532c4fbc0f3
all_test_data = read_h_time_series(test_experiment)

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
function f(h, params, t)
	if h < params.h_hole
		return 0.0
	end
	return - π * params.r_hole ^ 2 * params.c * 
		sqrt(2 * g * (h .- params.h_hole)) / params.A_of_h(h)
end

# ╔═╡ 05ed4187-a01a-4a16-a0e7-b3867d252578
function simulate(h₀, params::NamedTuple, tf::Float64, callback=nothing)
	prob = ODEProblem(f, h₀, (0.0, tf), params)
	h_of_t = solve(
		prob, Tsit5(), saveat=0.5, 
		reltol=1e-6, abstol=1e-6, 
		callback=callback
	)
	return h_of_t
end

# ╔═╡ 6f7d4335-d9a6-4896-9d69-bfc1c2c1c3d0
md"""
## minimize loss (classical approach) to identify $c$
"""

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

	h_of_t = simulate(
		data[1, "h [cm]"], # h₀
		params,
		data[end, "t [s]"] # end of time
	)
	
	cost = 0.0
	for i in 1:nrow(data)
		tᵢ = data[i, "t [s]"]
		cost += (data[i, "h [cm]"] - h_of_t(tᵢ)) ^ 2
	end
	
	return cost
end

# ╔═╡ 5e79d8e1-429c-414f-b3a6-8cf4b93d1336
md"maximum likelihood estimate of $c$"

# ╔═╡ 0b75c073-f167-4553-b746-539a14cfcf25
loss(train_data, 0.4, tank_measurements)

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
sim_train_data = DataFrame(
	simulate(
		train_data[1, "h [cm]"], # h₀
		params,
		1.1 * train_data[end, "t [s]"] # end of time
	)
)

# ╔═╡ d25cda2f-6ec6-4c93-8860-f5ce9c3ee629
md"""
## visualize dynamic model (classically fit)
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
	
	scatter!(
		data[:, "t [s]"], 
		data[:, "h [cm]"], 
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
viz_sim_fit(train_data, sim_train_data)

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
@model function forward_model(
	data::DataFrame, 
	tm::TankMeasurements; 
	prior_only::Bool=false
)
	#=
	prior distributions
	=#
	# defines variance for measuring length with measuring tape
	σ_ℓ ~ Uniform(0.0, 0.5) # cm
	
	# bottom, top tank area measurements
	# std of product of two Guassians
	#   https://ccrma.stanford.edu/~jos/sasp/Product_Two_Gaussian_PDFs.html
	A_b ~ Normal(tm.A_b, σ_ℓ ^ 2 / 2) # cm²
	A_t ~ Normal(tm.A_t, σ_ℓ ^ 2 / 2) # cm²

	# height of tank
	H ~ Truncated(
		Normal(tm.H, σ_ℓ),
		tm.H - σ_ℓ, tm.H + σ_ℓ
	) # cm

	# radius of the hole. std 2%
	r_hole ~ Truncated(
		Normal(tm.r_hole, 0.02 * tank_measurements.r_hole),
		0.9 * tm.r_hole, 1.1 * tm.r_hole
	) # cm

	# discharge coefficient. Wikipedia says 0.65 for water.
	c ~ Truncated(Normal(0.65, 0.25), 0.1, 1.0) # unitless

	# height of the hole
	h_hole ~ Truncated(
		Normal(tm.h_hole, σ_ℓ),
		tm.h_hole - σ_ℓ, tm.h_hole + σ_ℓ
	) # cm
	
	# defines variance of liquid level sensor
	#   (treated as an unknown and inferred)
	σ ~ Uniform(0.0, 1.0) # cm

	# initial liquid level
	h₀_obs = data[1, "h [cm]"] # cm
	h₀ ~ Truncated(
		Normal(h₀_obs, σ),
		0.0, H
	)

	# do not use the rest of the data if doing prior only.
	if prior_only
		return
	end

	#=
	set up dynamic model for h(t)
	=#
	function my_A_of_h(h)
		if h < 0.0 error("h < 0") end
		if h > H error("h > H") end
		return h / H * A_t + (1 - h / H) * A_b
	end
	
	# parameters for ODE solver
	params = (
			  r_hole=r_hole,
			  c=c,
			  h_hole=h_hole,
			  A_of_h=my_A_of_h
			)
	
	# set up and solve ODE
	h_of_t = simulate(h₀, params, 1.1 * data[end, "t [s]"])
	
	#=
	code up likelihood
	=#
	for i in 2:nrow(data) # start at 2 b/c IC handled with informative prior
		tᵢ = data[i, "t [s]"]
		ĥᵢ = h_of_t(tᵢ, continuity=:right)[1]
		data[i, "h [cm]"] ~ Normal(ĥᵢ, σ)
	end

	return nothing
end

# ╔═╡ 8082559e-a5b0-41a8-b8ed-aec3b09e5b2b
begin
	train_model = forward_model(train_data, tank_measurements)
	
	train_posterior = DataFrame(
		sample(
			train_model, 
			NUTS(0.65), MCMCSerial(), n_MC_sample, 3; progress=true
		)
	)
end

# ╔═╡ 2ee1ca40-141f-40ad-b4c1-a2e025f69f95
md"make sure never $h_0>H$."

# ╔═╡ c2d877b5-d309-4868-925d-dab8d7d23403
@assert all(train_posterior[:, "H"] .>= train_posterior[:, "h₀"])

# ╔═╡ c239deed-8291-45aa-95cf-94df26e0136d
md"""
## visualize posterior distribution of model parameters
"""

# ╔═╡ ccb1f005-567d-47f8-bec1-8db268d878ec
inferred_params = ["σ_ℓ", "A_b", "A_t", "H", "r_hole", "c", "h_hole", "σ", "h₀"]

# ╔═╡ 5bb0b72a-8c77-4fcb-bbde-d144986d9c1e
function viz_posterior(posterior::DataFrame, params::Vector{String},
			           tm::TankMeasurements, h₀_obs::Float64
)
	params_to_title = Dict(
						"A_b" => "A, tank bottom", 
						"A_t" => "A, tank top",
						"r_hole" => "hole radius",
						"c" => "discharge coefficient",
						"h_hole" => "hole height",
						"h₀" => "h₀",
						"σ_ℓ" => "std length measurement",
						"H" => "tank height", 
						"σ" => "std level sensor"
	)
	
	fig = Figure(size=(600, 600))
	i, j = 1, 1 # row, col
	 
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
		
		if j == 3
			i += 1
			j = 1
		else
			j += 1
		end
	end
	
	return fig
end

# ╔═╡ ded5b462-06dd-43a4-93b0-c52ad87174eb
viz_posterior(train_posterior, inferred_params, tank_measurements, train_data[1, "h [cm]"])

# ╔═╡ 86b56683-c80e-4c0f-8b03-a4869860d04f
md"## posterior predictive check"

# ╔═╡ 2ab35999-3615-4f5c-8d89-36d77802fe9b
function viz_fit(posterior::DataFrame, data::DataFrame; 
				savename::Union{String, Nothing}=nothing, 
				n_sample::Int=75, only_ic::Bool=false
)
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="time, t [s]", 
		ylabel="liquid height, h [cm]"
	)
	
	ts = range(0, 1.05 * maximum(data[:, "t [s]"]), length=500)
	tspan = (0.0,  maximum(ts) * 1.05)
	

	# sample posterior models
	for i in sample(1:nrow(posterior), n_sample)
		# area of tank
		A_of_h_tank = h -> h / posterior[i, "H"] * posterior[i, "A_t"] + 
			  	    (1 - h / posterior[i, "H"]) * posterior[i, "A_b"]

		# area of object
		if "As[1]" in names(posterior)
			N = sum(contains.(names(posterior), "As"))
			
			Aₒs = [posterior[i, "As[$n]"] for n in 1:N]
			hs = range(
				posterior[i, "h_hole"], posterior[i, "h₀"], length=N
			)
			
			A_of_h_object = linear_interpolation(hs, Aₒs)
		else
			A_of_h_object(h) = 0.0
		end
		
		params = (
			  r_hole=posterior[i, "r_hole"],
			  c=posterior[i, "c"],
			  h_hole=posterior[i, "h_hole"],
			  A_of_h=h -> A_of_h_tank(h) - A_of_h_object(h)
			)

		# set up, solve ODE
		sim_data = DataFrame(
			simulate(posterior[i, "h₀"], params, 1.05 * maximum(ts))
		)
			
		lines!(
			sim_data[:, "timestamp"], sim_data[:, "value"], 
			label="model", color=(colors["model"], 0.1)
		)

		# h hole
		hlines!(ax, posterior[i, "h_hole"], color=("gray", 0.1), linestyle=:dash)
	end	
	
	scatter!(
		data[only_ic ? 1 : 1:end, "t [s]"], 
		data[only_ic ? 1 : 1:end, "h [cm]"],
		label="data",
		color=colors["data"]
	)
	axislegend(unique=true)
	ylims!(0, nothing)
	xlims!(0, maximum(ts))
	if savename!=nothing
		save( "$savename.pdf", fig)
	end
	return fig
end

# ╔═╡ 2a01b228-f281-46c4-9764-fac6cc1b4217
viz_fit(train_posterior, train_data, savename="posterior_train")

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
		ylabel="liquid level, h [cm]",
	)

	# sample from the train posterior
	emptying_time = Float64[]
	for i in sample(1:nrow(posterior), n_sample)
		# check for first instance when the liquid level
		#  is the same as the height of the hole in the base
		condition(h, t, integrator) = h[1] - posterior[i, "h_hole"]
		
		# retrive the emptying time [t] when h(t) = h_hole
		function affect!(integrator)
			push!(emptying_time, integrator.t)
		end
		cb = ContinuousCallback(condition, affect!)

		params = (
			  r_hole=posterior[i, "r_hole"],
			  c=posterior[i, "c"],
			  h_hole=posterior[i, "h_hole"],
			  A_of_h=h -> h > posterior[i, "H"] ? 
			  	error("h > H") : 
			   	h / posterior[i, "H"] * posterior[i, "A_t"] + 
			  	    (1 - h / posterior[i, "H"]) * posterior[i, "A_b"]
		)

		# sample an initial condtion
		h₀_obs = test_data[1, "h [cm]"]
		h₀_distn = Truncated(
			Normal(h₀_obs, posterior[i, "σ"]),
			0.0, posterior[i, "H"]
		)
		h₀ = rand(h₀_distn)

		# simulate trajectory
		sim_data = DataFrame(
			simulate(h₀, params, 1.25 * test_data[end, "t [s]"], cb)
		)
			
		# plot trajectories
		lines!(ax, sim_data[:, "timestamp"], sim_data[:, "value"], 
			label="model", color=(colors["model"], 0.1))
		# h hole
		hlines!(ax, posterior[i, "h_hole"], color=("gray", 0.1), 
			linestyle=:dash)
	end

	# plot dist'n of emptying times
	hist!(ax_stopping, emptying_time, color=Cycled(3))

	scatter!(
		ax,
		test_data[:, "t [s]"], 
		test_data[:, "h [cm]"],
		label="experiment",
		color=colors["data"]
		)
	# @assert 1.25 * test_data[end, "t [s]"] > maximum(emptying_time)
	xlims!(0, 1.01 * maximum(emptying_time))
	ylims!(ax, 0, nothing)
	ylims!(ax_stopping, 0, nothing)

	linkxaxes!(ax, ax_stopping)
	axislegend(ax, unique=true)
	
	if ! isnothing(savename)
		save("$savename.pdf", fig)
	end
	
	return fig 
end	

# ╔═╡ a3ba0c9d-5f81-4023-9ce0-ff29536aa968
viz_test(train_posterior, test_data, savename="test")

# ╔═╡ 67c46219-2183-47fd-bd3c-82facff98d53
md"## visualize prior too"

# ╔═╡ 38a03e22-d596-4227-a9de-3ef54dc7256e
begin
	train_model_prior = forward_model(train_data, tank_measurements, prior_only=true)
	
	train_prior = DataFrame(
		sample(
			train_model_prior, 
			NUTS(0.65), MCMCSerial(), n_MC_sample, 3; progress=true
		)
	)
end

# ╔═╡ 2148cbb2-b41f-4df3-ab5c-55a89eff7bf1
viz_fit(train_prior, train_data, savename="prior_train", only_ic=true)

# ╔═╡ 9533c662-80af-4dd4-bf25-02e894867360
md"""
# Bayesian inference of object shape

## read experimental data
"""

# ╔═╡ b06a1c07-6250-4324-8802-010e5d847edb
begin
	data_w_object_filenames = ["obs_4_18_1.csv", "obs_4_18_2.csv"]
	@bind data_w_object_filename Select(data_w_object_filenames)
end

# ╔═╡ 8b6d766a-8f7b-4b9a-9a15-0f7375087120
all_data_w_object = read_h_time_series(data_w_object_filename)

# ╔═╡ 807222ba-5ff8-4f33-a9a0-7c69b1dccf52
data_w_object = downsample(all_data_w_object, n_data_sample)

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
		save("h_of_t_with_without_object.pdf", fig)
		fig
	end
	
	viz_free_vs_occupied(data_w_object, train_data)
end

# ╔═╡ 580de17a-625d-420e-974c-86766197025e
md"## measured area of object

a ground truth"

# ╔═╡ cb59f55b-c748-4a94-b344-e50a8fa7c690
begin
	object_true_area = CSV.read("obstacle_area.csv", DataFrame)
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
		r_hole = tank_measurements.r_hole,
		# fudge factor
		c = c_opt,
		# height of the hole
		h_hole = tank_measurements.h_hole,
		# area as a function of h"
		A_of_h = h -> A_of_h(h, tank_measurements) - A_of_object(h)
	)

# ╔═╡ b12963ae-bf7d-4ef7-b1a8-e2d1e24f9b4b
sim_object_data = DataFrame(
	simulate(
		data_w_object[1, "h [cm]"], object_params, data_w_object[end, "t [s]"]
	)
)

# ╔═╡ cfbe753d-85a8-445f-9eda-14a376d7e0c6
viz_sim_fit(data_w_object, sim_object_data)

# ╔═╡ c53edeef-324a-418f-907d-aaf557cb8d24
md"## classical method


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
		row -> row["h [cm]"] >= tank_measurements.h_hole + 1.0, data_w_object
	)
	
	h_of_t_object = Spline1D(
		data_w_object_spline_fit[:, "t [s]"], data_w_object_spline_fit[:, "h [cm]"]; k=4, s=20.0, bc="nearest"
	)
end

# ╔═╡ f3f886d6-3010-4dd9-b42a-d5309463beb6
h_of_t_object(30.0)

# ╔═╡ 5003d6c6-fa30-423e-80c1-a0f82e4085b9
derivative(h_of_t_object, 30.0) # WARNING: doesn't extrapolate.

# ╔═╡ 36b1822e-fe08-494b-a57d-5888163a7b54
function viz_spline_fit(data_w_object::DataFrame, h_of_t_object::Spline1D, h_hole::Float64)
	fig = Figure(size=(500, 600))
	axs = [Axis(fig[i, 1]) for i = 1:2]
	linkxaxes!(axs...)
	axs[2].xlabel = "t [s]"
	axs[1].ylabel = "h(t) [cm]"
	axs[2].ylabel = "dh/dt [cm/s]"
	
	ts = range(0.0, data_w_object[end, "t [s]"], length=150)

	# top
	scatter!(axs[1], data_w_object[:, "t [s]"], data_w_object[:, "h [cm]"])
	hs = h_of_t_object.(ts)
	lines!(axs[1], ts, hs, color=Cycled(3))
	ylims!(axs[1], 0, nothing)

	# bottom
	h′ = [derivative(h_of_t_object, tᵢ) for tᵢ in ts]
	# flatten derivative since it doesn't handle extrapolation
	h′[h′ .> 0.0] .= 0.0
	# shut off derivative
	lines!(axs[2], ts, h′)
	
	for ax in axs
		xlims!(ax, 0, nothing)
	end
	# save("block_data_spline_fit.png", fig)
	fig
end

# ╔═╡ 29ccc89b-f76a-44fa-8a89-e3ca10742ba1
viz_spline_fit(data_w_object, h_of_t_object, tank_measurements.h_hole)

# ╔═╡ 29d3cc8f-780b-449a-87ba-8d543ad2473b
function classical_soln_A_object(
	h_of_t_object::Spline1D,
	c_opt::Float64,
	tank_measurements::TankMeasurements
)
	ts = range(0.0, 1000.0, length=150)

	# dh/dt
	h′ = [derivative(h_of_t_object, tᵢ) for tᵢ in ts]
	h′[h′ .> 0.0] .= 0.0 # flatten to handle BC

	# h
	hs = h_of_t_object.(ts)

	# inferred area in tank
	A_tank = hs / tank_measurements.H * tank_measurements.A_t .+ 
		(1 .- hs / tank_measurements.H) * tank_measurements.A_b
	Aₒ = A_tank .+ π * tank_measurements.r_hole ^ 2 * c_opt * 
			sqrt.(2 * g * (hs .- tank_measurements.h_hole)) ./ h′
	
	inferred_object_data = DataFrame("h [cm]" => hs, "Aₒ [cm²]" => Aₒ)
	
	# filter data where (i) dh/dt
	filter!(row -> ! isinf(row["Aₒ [cm²]"]), inferred_object_data)
	
	return inferred_object_data
end

# ╔═╡ 0a48e016-2fba-47cc-a212-47b4a3324b20
classical_Aₒ = classical_soln_A_object(h_of_t_object, c_opt, tank_measurements)

# ╔═╡ 5feb46c0-3888-4586-8b12-f990d4d38912
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="water level, h [cm]", 
		ylabel="cross-sectional area\nof object\nAₒ [cm²]"
	)
	
	lines!(
		classical_Aₒ[:, "h [cm]"], classical_Aₒ[:, "Aₒ [cm²]"], 
		label="predicted", color=Cycled(4)
	)
	scatter!(
		object_true_area[:, "h [cm]"], object_true_area[:, "area [cm²]"], 
		label="measured"
	)
	
	xlims!(0, nothing)
	ylims!(-10, 100)
	hlines!(0.0, color="lightgray")
	axislegend()
	fig
end

# ╔═╡ b23dc763-d91f-4d66-94d2-dcf96cb07f54
md"## Bayesian inference

### forward model
"

# ╔═╡ 798d8d16-1c19-400d-8a94-e08c7f991e33
@model function forward_model_object(
	data_w_object::DataFrame, train_posterior::DataFrame, γ::Float64, N::Int;
	prior_only::Bool=false
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
	H ~ Truncated(
		Normal(mean(train_posterior.H), std(train_posterior.H)), 
		mean(train_posterior.H) - 4 * std(train_posterior.H),
		mean(train_posterior.H) + 4 * std(train_posterior.H)
	)

	# radius of the hole
	r_hole ~ Truncated(
		Normal(mean(train_posterior.r_hole), std(train_posterior.r_hole)),
		mean(train_posterior.r_hole) - 4 * std(train_posterior.r_hole),
		mean(train_posterior.r_hole) + 4 * std(train_posterior.r_hole)
	)
		

	# discharge coefficient. 
	c ~ Truncated(
		Normal(mean(train_posterior.c), std(train_posterior.c)),
		mean(train_posterior.c) - 2 * std(train_posterior.c),
		mean(train_posterior.c) + 2 * std(train_posterior.c)
	)

	# height of the hole
	h_hole ~ Truncated(
		Normal(mean(train_posterior.h_hole), std(train_posterior.h_hole)),
		mean(train_posterior.h_hole) - 4 * std(train_posterior.h_hole),
		mean(train_posterior.h_hole) + 4 * std(train_posterior.h_hole)
	)

	# initial liquid level
	h₀_obs = data_w_object[1, "h [cm]"]
	h₀ ~ Truncated(
		Normal(h₀_obs, σ), 
		h₀_obs - 4 * σ, 
		0.999 * H
	)
		
	# bottom, top tank area measurements
	A_b ~ Truncated(
		Normal(mean(train_posterior.A_b), std(train_posterior.A_b)),
		mean(train_posterior.A_b) - 4 * std(train_posterior.A_b),
		mean(train_posterior.A_b) + 4 * std(train_posterior.A_b)
	)
	
	A_t ~ Truncated(
		Normal(mean(train_posterior.A_t), std(train_posterior.A_t)),
		mean(train_posterior.A_t) - 4 * std(train_posterior.A_t),
		mean(train_posterior.A_t) + 4 * std(train_posterior.A_t)
	)

	# Random distribution of obstacle area grid at points N 
	As ~ filldist(Normal(), N) # cm²

	function A_of_tank(h)
		h < 0 ? error("h < 0") : nothing
		h > H ? error("h > H") : nothing
		return h / H * A_t + (1 - h / H) * A_b
	end

	# corresponding h's for the unknown A's
	hs = range(h_hole, h₀, length=N)
	# prior: object could be ANY size
	As[1] ~ Uniform(0.0, 0.99 * A_of_tank(hs[1]))
	for i in 2:N
		As[i] ~ Truncated(
						  As[i - 1] + Normal(0.0, γ), 
						  0.0, 
						  0.99 * A_of_tank(hs[i])
						 )
	end

	# for prior, do not show the algo the data :)
	if prior_only
		return
	end

	#=
	set up dynamic model for h(t)
	=#
	A_of_object = linear_interpolation(hs, As)

	# parameter for ODE solver
	params = (
			  r_hole=r_hole,
			  c=c,
			  h_hole=h_hole,
			  A_of_h=h ->  A_of_tank(h) - A_of_object(h)
			)
	
	# set up, solve ODE
	h_of_t = simulate(h₀, params, 1000.0)

	# observations.
	for i in 2:nrow(data_w_object)
		tᵢ = data_w_object[i, "t [s]"]
		ĥᵢ = h_of_t(tᵢ, continuity=:right)[1]
		data_w_object[i, "h [cm]"] ~ Normal(ĥᵢ, σ)
	end
	
	return nothing
end

# ╔═╡ da44647a-36e4-4116-9698-df1cb059c2b7
md"### posterior"

# ╔═╡ fb3ece76-f85c-41e1-a332-12c71d9d3cc0
begin
	γ = 5.0 # smoothness param
	N = 15  # number of points to infer area on
end

# ╔═╡ 02939a87-e811-4ae4-8b6b-173370029889
begin
	object_tank_model = forward_model_object(data_w_object, train_posterior, γ, N)
	object_posterior = DataFrame(
		sample(object_tank_model, NUTS(0.65), MCMCSerial(), 
			n_MC_sample, 3; progress=true
		)
	)
end

# ╔═╡ 3c9a219f-74ef-45fb-83e7-c497e0bee362
@assert all(object_posterior[:, "h₀"] .< object_posterior[:, "H"])

# ╔═╡ e1264f57-f675-4f37-b4db-313cfc52ab8e
viz_fit(object_posterior, data_w_object, savename="posterior_object")

# ╔═╡ a127225a-5b79-4074-a16b-cecd11030800
function viz_inferred_area(
	object_posterior::DataFrame, 
	object_true_area::DataFrame, 
	γ::Float64, 
	N::Int
)
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="height, h [cm]", 
		ylabel="area, a′ [cm²]"
	)
	
	for i in 1:nrow(object_posterior)
		Aₒs = [object_posterior[i, "As[$n]"] for n in 1:N]
		hs = range(
			object_posterior[i, "h_hole"], object_posterior[1, "h₀"], length=N
		)
		lines!(hs, Aₒs, label="model", color=(theme_colors[8], 0.1))
	end

	scatter!(object_true_area[:, "h [cm]"], object_true_area[:, "area [cm²]"], 
			label="data", color=colors["data"])

	ylims!(0, tank_measurements.A_t)

	axislegend("γ=$γ; N=$N", unique=true, position=:rt, titlefont="normal")

	return fig
end

# ╔═╡ b43f9f58-94fd-4c92-8e91-9a6b86cfc041
viz_inferred_area(object_posterior, object_true_area, γ, N)

# ╔═╡ bd95428d-1077-4417-bfca-0c5da7378af2
md"### prior"

# ╔═╡ 65d81268-9ff2-4a18-b0ce-4b105740dc8b
begin
	object_tank_model_prior = forward_model_object(data_w_object, train_posterior, γ, N, prior_only=true)
	object_prior = DataFrame(
		sample(object_tank_model_prior, NUTS(0.65), MCMCSerial(), 
			n_MC_sample, 3; progress=true
		)
	)
end

# ╔═╡ 8c1d1401-bc6b-4be3-8481-1c9a8f86f63d
viz_inferred_area(object_prior, object_true_area, γ, N)

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
# ╟─6ebe0cb0-ba35-411c-9a7a-a8b6eecf326f
# ╠═385442da-f101-4ddf-8293-46d71d6a48fc
# ╟─ef43e50a-5af8-4733-88a4-cd159d173034
# ╠═9dabad13-cfa4-4e06-950d-f7c7d96c1147
# ╟─e040094c-7511-4831-b94a-1c1185868202
# ╠═23ee0e85-a84b-4b63-b432-5526559efcee
# ╟─078c01f7-e47e-4af0-be1c-ac4527b735fd
# ╠═8b00d2b3-9182-42ab-8393-91707b813f60
# ╟─7899f488-9c48-466f-857d-f5a31b5820ab
# ╟─96f26378-846c-4964-935c-0372e2e86e91
# ╠═2ddf387c-5a61-4490-9746-96e1589c7a74
# ╠═b2228d5c-16b4-4fee-b9b6-1112d7cf391c
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
# ╟─a1a10e2f-1b78-4b93-9295-7c0055e32692
# ╠═58eff13c-44b5-4f19-8a42-cf9907ac9515
# ╟─8a21fa0f-d3c3-4aa2-8b8b-74001d921c4a
# ╠═8f5b8859-6b8c-4f2a-af3a-b13c2d33fe2a
# ╠═8082559e-a5b0-41a8-b8ed-aec3b09e5b2b
# ╟─2ee1ca40-141f-40ad-b4c1-a2e025f69f95
# ╠═c2d877b5-d309-4868-925d-dab8d7d23403
# ╟─c239deed-8291-45aa-95cf-94df26e0136d
# ╠═ccb1f005-567d-47f8-bec1-8db268d878ec
# ╠═5bb0b72a-8c77-4fcb-bbde-d144986d9c1e
# ╠═ded5b462-06dd-43a4-93b0-c52ad87174eb
# ╟─86b56683-c80e-4c0f-8b03-a4869860d04f
# ╠═2ab35999-3615-4f5c-8d89-36d77802fe9b
# ╠═2a01b228-f281-46c4-9764-fac6cc1b4217
# ╟─a5ae695b-bfc0-4425-9b64-bbeeba7da015
# ╠═eaf470e9-2898-41d5-a6d5-4cd846e9c0de
# ╠═a3ba0c9d-5f81-4023-9ce0-ff29536aa968
# ╟─67c46219-2183-47fd-bd3c-82facff98d53
# ╠═38a03e22-d596-4227-a9de-3ef54dc7256e
# ╠═2148cbb2-b41f-4df3-ab5c-55a89eff7bf1
# ╟─9533c662-80af-4dd4-bf25-02e894867360
# ╠═b06a1c07-6250-4324-8802-010e5d847edb
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
# ╟─c53edeef-324a-418f-907d-aaf557cb8d24
# ╠═4831a3be-35d3-420c-8463-bb14a597cc6a
# ╠═f3f886d6-3010-4dd9-b42a-d5309463beb6
# ╠═5003d6c6-fa30-423e-80c1-a0f82e4085b9
# ╠═36b1822e-fe08-494b-a57d-5888163a7b54
# ╠═29ccc89b-f76a-44fa-8a89-e3ca10742ba1
# ╠═29d3cc8f-780b-449a-87ba-8d543ad2473b
# ╠═0a48e016-2fba-47cc-a212-47b4a3324b20
# ╠═5feb46c0-3888-4586-8b12-f990d4d38912
# ╟─b23dc763-d91f-4d66-94d2-dcf96cb07f54
# ╠═798d8d16-1c19-400d-8a94-e08c7f991e33
# ╟─da44647a-36e4-4116-9698-df1cb059c2b7
# ╠═fb3ece76-f85c-41e1-a332-12c71d9d3cc0
# ╠═02939a87-e811-4ae4-8b6b-173370029889
# ╠═3c9a219f-74ef-45fb-83e7-c497e0bee362
# ╠═e1264f57-f675-4f37-b4db-313cfc52ab8e
# ╠═a127225a-5b79-4074-a16b-cecd11030800
# ╠═b43f9f58-94fd-4c92-8e91-9a6b86cfc041
# ╟─bd95428d-1077-4417-bfca-0c5da7378af2
# ╠═65d81268-9ff2-4a18-b0ce-4b105740dc8b
# ╠═8c1d1401-bc6b-4be3-8481-1c9a8f86f63d
