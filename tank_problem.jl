### A Pluto.jl notebook ###
# v0.19.38

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
    using CSV, Interpolations, DataFrames, CairoMakie, DifferentialEquations, Turing, StatsBase, PlutoUI, Distributions
end

# ╔═╡ e1b54f1e-3c1b-4657-ae3d-d68a78617a4d
import AlgebraOfGraphics

# ╔═╡ 4391f124-cbef-46e5-8462-e4e5126f5b38
begin
	AlgebraOfGraphics.set_aog_theme!(
		fonts=[AlgebraOfGraphics.firasans("Light"), 
			   AlgebraOfGraphics.firasans("Light")]
	)
	resolution = (0.9*500, 0.9*380)
	update_theme!(
		fontsize=20, 
		linewidth=4,
		markersize=14,
		titlefont=AlgebraOfGraphics.firasans("Light"),
		# resolution=resolution
	)
end

# ╔═╡ 245836a9-6b44-4639-9209-e7ad9035e293
TableOfContents()

# ╔═╡ 418525b7-c358-41da-b865-5df3feb15855
md"
## Sensor Calibration

read in data characterizing the calibration curve of the liquid level sensor.
"

# ╔═╡ a95e371e-9319-4c7e-b5d9-4c4a50d12cd7
calibration = CSV.read("liquid calibration.csv", DataFrame)

# ╔═╡ 23ee0e85-a84b-4b63-b432-5526559efcee
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="level sensor reading", 
		ylabel="liquid level [cm]"
	)
	scatterlines!(
		calibration[:, "level sensor reading"], 
		calibration[:, "liquid level (cm)"]
	)
	fig
end

# ╔═╡ 7752316d-9dd0-4403-aa08-22c977ff3727
md"""
## Tank Geometry
"""

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
	local δ = (L_t - L_b) / 2
	H = sqrt(H★ ^ 2 - δ ^ 2)

	# solve for the r consistent with a rounded rectangle [cm]
	r_b = (p_b / 2 - (L_b + W_b)) / (π - 4)
	r_t = (p_t / 2 - (L_t + W_t)) / (π - 4)
	@show r_b, r_t

	# finally, calculate areas.
	A_b = (L_b - 2 * r_b) * (W_b - 2 * r_b) + 
		2 * r_b * (L_b + W_b - 4 * r_b) + π * r_b ^ 2
	A_t = (L_t - 2 * r_t) * (W_t - 2 * r_t) + 
		2 * r_t * (L_t + W_t - 4 * r_t) + π * r_t ^ 2
	@show A_b, A_t

	function A_of_h(h)
		if h < 0 || h > H
			error("tank over/under-flow!")
		end
		return h / H * A_t + (1 - h / H) * A_b
	end
end

# ╔═╡ 078c01f7-e47e-4af0-be1c-ac4527b735fd
md"""
## Data Preprocessing

the liquid level sensor is up against the side of the tank hence has a bit of a slant. we correct for this.
"""

# ╔═╡ 9af216b7-4bf2-42fb-bd95-5b2040d019a7
function map_slant_height_to_perpendicular(h★)
	return H * h★ / H★
end

# ╔═╡ d3eeccff-39b3-429d-aec6-e7f1a500b729
map_slant_height_to_perpendicular(1.0) # note, this correction is TINY!

# ╔═╡ 6e922ea4-01d0-4202-baed-3cb693b663bc
function level_strip_reading_to_h(sensor_reading::Array; 
								  calibration::DataFrame=calibration)
	mapping = linear_interpolation(
							calibration[:, "level sensor reading"], 
							calibration[:, "liquid level (cm)"]
							)
	ĥ = mapping(sensor_reading) # convert sensor reading to slant height
	
	h =map_slant_height_to_perpendicular(ĥ) # convert slant height to perpendicular 
	return h
end

# ╔═╡ 759852ee-50e7-4deb-ac7e-c4693103c2a7
begin
	experiments = ["experiment 3- 1-12.csv", "experiment 4- 1-12.csv"]
	@bind experiment Select(experiments)
end

# ╔═╡ b06a1c07-6250-4324-8802-010e5d847edb
begin
	obstruction_data = ["obstruction-1-2-28.csv", "obstruction-2-2-28.csv", 
						"obstruction-3-2-28.csv"]
	@bind obstruction Select(obstruction_data)
end

# ╔═╡ 8b00d2b3-9182-42ab-8393-91707b813f60
function read_h_time_series(file::String)
	data = CSV.read(file, DataFrame) # read in file
	
	rename!(data, ["Time [s]", "liquid_level_reading"])
	
	id = argmax(data[:, "liquid_level_reading"]) # index of highest water level
	data = data[id:end, :] # index of highest water level
	data[:, "Time [s]"] = data[:, "Time [s]"] .- data[1, "Time [s]"]
	
	data[:, "liquid level [cm]"] = level_strip_reading_to_h(
												data[:, "liquid_level_reading"]
												)
	
	# data[:, "p_liquid level [cm]"] = sensor_height_to_perpendicular.(
	# 											data[:, "liquid level [cm]"]
	# 											)
	
	return data[:, ["Time [s]", "liquid level [cm]"]]
end

# ╔═╡ 2a3d5a3f-5d20-4efb-91fd-cc84e77d5922
test_data_filename = symdiff(experiments, [experiment])[1]

# ╔═╡ f0675774-b9ba-40a9-953f-ce5784454378
experiment

# ╔═╡ 5100771a-38b6-437a-91e4-70495391067b
begin
	train_data = read_h_time_series(experiment)
	test_data = read_h_time_series(test_data_filename)
end

# ╔═╡ 9a7e5903-69be-4e0a-8514-3e05feedfed5
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="water level, h [cm]", 
		ylabel="area of water surface at top, A(h) [cm²]"
	)

	water_level = train_data[:, "liquid level [cm]"]
	lines!(water_level, A_of_h.(water_level))
	ylims!(0, nothing)
	fig
end

# ╔═╡ 8b6d766a-8f7b-4b9a-9a15-0f7375087120
block_data = read_h_time_series(obstruction)

# ╔═╡ df589086-cc94-4da9-b982-563cd6ecc643
md"""
## Tank Geometry 
"""

# ╔═╡ 33f889d7-e875-40d8-9d6d-cc87b0fbaf22
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1], xlabel="Time [s]", ylabel="water level [cm]")
	lines!(
			train_data[:, "Time [s]"], 
			train_data[:, "liquid level [cm]"], 
			label="no obstruction"
		)
	
	# lines!(
	# 		block_data[:, "Time [s]"], 
	# 		block_data[:, "liquid level [cm]"], 
	# 		label="obstruction"
	# 	)
	
	axislegend()
	fig
end

# ╔═╡ e379461f-8896-4e2a-a71b-1871a8a37eb5
md"""
## Differential tank model setup
"""

# ╔═╡ 31306e0b-9748-48a8-b9d2-892cb501b7ba
begin
	g = 980.665 # cm/s²
	h0 = [train_data[1, "liquid level [cm]"]]
	tspan = (0, train_data[end, "Time [s]"])
	 
	measurements = (
			  a = 0.079, # d = 0.3175; a = (π * d²) / 4  [cm²] 
			  c = 0.63, 
			  h_hole = 2.3, # cm 
			  A_of_h = A_of_h # A(h)
			)

end

# ╔═╡ c6a263eb-cb45-4ee7-9c02-549c89298652
f(h, p, t) = - p.a .* sqrt.(max.(2 * g * (h .- p.h_hole), 0.0)) .* p.c ./ p.A_of_h(h)

# ╔═╡ a934a309-6ab6-46e0-b448-2572334f5b7a
# function condition(h, t, integrator) # Event when condition(u,t,integrator) == 0
#     h[1] = p.h_hole
# end

# ╔═╡ d3307918-1fdb-4f87-bb92-67330d22e58b
begin
	prob = ODEProblem(f, h0, tspan, measurements)
	condition(h, t, integrator) = h[1] <= measurements.h_hole
	affect!(integrator) = terminate!(integrator)
	cb = ContinuousCallback(condition, affect!)
	sol = solve(prob, Tsit5(), callback=cb)
end

# ╔═╡ 3853c2fb-3665-4c7b-8511-af01694f04b3


# ╔═╡ 444f6d74-273e-486d-905a-1443ec0e98df
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1], xlabel="Time [s]", ylabel="liquid level [cm]")
	scatter!(
			train_data[:, "Time [s]"], 
			train_data[:, "liquid level [cm]"], 
			label="experimental")
	
	ts = range(0, maximum(train_data[:, "Time [s]"]), length=1000)
	lines!(ts, [sol.(t, continuity=:right)[1] for t in ts], label="model", color=:red)
	
	axislegend()
	fig
end

# ╔═╡ a1a10e2f-1b78-4b93-9295-7c0055e32692
md"""
## Bayesian inference approach for measurement estimation  
"""

# ╔═╡ f21dc58e-d4e8-4314-b5dd-abbcb29efe86
function downsample(data::DataFrame, n::Int)
	ids = trunc.(Int, collect(range(1, nrow(data), length=n)))
	return data[ids, :]
end

# ╔═╡ 8f5b8859-6b8c-4f2a-af3a-b13c2d33fe2a
@model function infer_params(data, measurements)

	σ_measurement = 0.1
	
	# Prior distributions.
	A_b ~ Normal(A_bottom, σ_measurement) # cross-sectional area at the base of tank
	
	A_t ~ Normal(A_top, σ_measurement) # cross-sectional area at the top of tank
	
	a ~ TruncatedNormal(
						measurements.a, 
						0.01, 
						measurements.a - 0.02, 
						measurements.a + 0.02
						) # area of the orifice [cm²]
	
	c ~ TruncatedNormal(
						measurements.c, 
						σ_measurement, 
						measurements.c - 2 * σ_measurement, 
						measurements.c + 2 * σ_measurement
						) # discharge coefficient
	
	h_hole ~ TruncatedNormal(
							 measurements.h_hole, 
							 σ_measurement, 
						  	 measurements.h_hole - 2 * σ_measurement, 
						  	 measurements.h_hole + 2 * σ_measurement
							) # height of orifice from the base of tank
	
	σ ~ Uniform(0.0, 1.0) # measurement noise
	
	h0_obs = data[1, "liquid level [cm]"]
	h_0 ~ TruncatedNormal(h0_obs, σ, h0_obs - 2 * σ, h0_obs +  2 * σ) # initial liquid level


	# recalibrate area interpolation
	
	A_of_h = linear_interpolation([0.0, H_tank], [A_b, A_t])
	
	# parameter for ODE
	params = (
			  a = a, # area of the orifice [cm²]
			  c = c, # discharge coefficient
			  h_hole = h_hole, # height of hole to the base of tank [cm] 
			  A_of_h = A_of_h 
			)
	
	# set up ODE
	prob = ODEProblem(f, [h_0], tspan, params)

	# callbacks
	condition(h, t, integrator) = h[1] <= params.h_hole
	affect!(integrator) = terminate!(integrator)
	cb = ContinuousCallback(condition, affect!)
	@show params 
	sol = solve(prob, Tsit5(), callback=cb)
	
	# Observations.
	for i in 2:nrow(data)
		tᵢ = data[i, "Time [s]"]
		data[i, "liquid level [cm]"] ~ Normal(sol(tᵢ, continuity=:right)[1], σ)
	end

	return nothing
end



# ╔═╡ 8082559e-a5b0-41a8-b8ed-aec3b09e5b2b
begin
	inference_data = downsample(train_data, 12)
	model = infer_params(inference_data, measurements)
	
	chain = sample(model, NUTS(0.65), MCMCSerial(), 3, 3; progress=true)
end

# ╔═╡ 7ebe4680-c583-4f92-8bae-dd84c3fb5139
posterior = DataFrame(chain)

# ╔═╡ a2048a65-7b7c-41fc-b7ee-f9199f3e96b5
begin
	posterior_to_true = Dict(
							"A_b" => A_bottom, 
							"A_t" => A_top,
							"a" => measurements.a,
							"c" => measurements.c,
							"h_hole" => measurements.h_hole,
							"h_0" => h0[1]
							)
	
	params_to_title = Dict(
							"A_b" => "Area of the bottom tank", 
							"A_t" => "Area of the top tank",
							"a" => "Area of the oriface",
							"c" => "discharge coefficient",
							"h_hole" => "height of orifice from the base of tank",
							"h_0" => "initial water level"
							)
end

# ╔═╡ 5bb0b72a-8c77-4fcb-bbde-d144986d9c1e
function viz_posterior(posterior)
	fig = Figure()
	ax = [Axis(fig[i, j]) for i in 1:3, j in 1:2]
	j, i = 1, 1
	 
	for p in keys(posterior_to_true)
		# vizualize the distribution
		density!(ax[i, j], posterior[:, p], color = (:blue, 0.0), strokewidth = 3, 
				 strokecolor = :blue)

		# plot real value
		vlines!(ax[i, j], [posterior_to_true[p]], linestyle="--",
				 color=:black, label="true value")

		# plot 95% interval
		lo, hi = quantile(posterior[:, p], [0.025, 0.975])
		lines!(ax[i, j], [lo, hi], [0, 0], linewidth=5, color=:black)

		# hidedecorations!(ax[i, j], ticks=false)
		
		ax[i, j].title = params_to_title[p]
		 i += 1
		if i > 3
			j += 1
			i = 1
		end
	end
	return fig	
end

# ╔═╡ 2ab35999-3615-4f5c-8d89-36d77802fe9b
function viz_fit(posterior, data)
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="Time [s]", ylabel="liquid level [cm]")
	
	ts = range(0, maximum(data[:, "Time [s]"]), length=1000)
	
	
	ids_post = sample(1:nrow(posterior), 100; replace=false)
	for i in ids_post
		areas =  # cm^2
		A_of_h = linear_interpolation(
								[0.0, H_tank], 
								[posterior[i, "A_b"], posterior[i, "A_t"]]
								)
		
		params = (
			  a = posterior[i, "a"], # area of the orifice [cm²]
			  c = posterior[i, "c"], # discharge coefficient
			  h_hole = posterior[i, "h_hole"], # height of hole to the base of tank [cm] 
			  A_of_h = A_of_h # A(h)
			)

		
	
		# set up ODE
		_prob = ODEProblem(f, [posterior[i, "h_0"]], tspan, params)
		sol = solve(_prob, Tsit5())
			
		
		lines!(ts, [sol.(t)[1] for t in ts], label="model", color=(:green, 0.1))
	end

	scatter!(data[:, "Time [s]"], data[:, "liquid level [cm]"], label="experimental")
	
	axislegend(unique=true)
	return fig
end

# ╔═╡ 5bea087e-c241-424d-a526-25eae30bfe15
viz_posterior(posterior)

# ╔═╡ 765ff940-1328-4806-aeaa-de8a41a6f4df


# ╔═╡ 2a01b228-f281-46c4-9764-fac6cc1b4217
viz_fit(posterior, inference_data)

# ╔═╡ 193d0e02-988e-4f58-b57d-7a1d125069a4
@model function infer_test(test_data, posterior)
	
	# Prior distributions.
	A_b ~ Normal(mean(posterior.A_b), std(posterior.A_b)) # cross-sectional area at the base of tank
	
	A_t ~ Normal(mean(posterior.A_t), std(posterior.A_t)) # cross-sectional area at the top of tank
	
	a ~ TruncatedNormal(
						 mean(posterior.a), 
						 std(posterior.a), 
						 mean(posterior.a) - 2 * std(posterior.a), 
						 mean(posterior.a) + 2 * std(posterior.a)
						) # area of the orifice [cm²]
	
	c ~ TruncatedNormal(
						mean(posterior.c), 
						std(posterior.c), 
						mean(posterior.c) - 2 * std(posterior.c), 
						mean(posterior.c) + 2 * std(posterior.c)
						) # discharge coefficient
	
	h_hole ~ TruncatedNormal(
							 mean(posterior.h_hole), 
							 std(posterior.h_hole),
							 mean(posterior.h_hole) - 2 * std(posterior.h_hole), 
						     mean(posterior.h_hole) + 2 * std(posterior.h_hole)
							) # height of orifice from the base of tank
	
	σ ~ Uniform(0.0, 1.0) # measurement noise
	
	h0_obs = test_data[1, "liquid level [cm]"]
	h_0 ~ TruncatedNormal(h0_obs, σ, h0_obs - 2 * σ, h0_obs +  2 * σ) # initial liquid level


	# recalibrate area interpolation
	A_of_h = linear_interpolation([0.0, H_tank], [A_b, A_t])
	
	# parameter for ODE
	params = (
			  a = a, # area of the orifice [cm²]
			  c = c, # discharge coefficient
			  h_hole = h_hole, # height of hole to the base of tank [cm] 
			  A_of_h = A_of_h # A(h)
			)
	tspan = (0, test_data[end, "Time [s]"])
	
	
	# set up ODE
	prob = ODEProblem(f, [h_0], tspan, params)
	
	# callbacks
	condition(h, t, integrator) = h[1] <= params.h_hole
	affect!(integrator) = terminate!(integrator)
	cb = ContinuousCallback(condition, affect!)
	# @show p 
	sol = solve(prob, Tsit5(), callback=cb)
	
	# Observations.
	tᵢ = test_data[1, "Time [s]"]
	test_data[1, "liquid level [cm]"] ~ Normal(sol(tᵢ)[1], σ)
	

	return nothing
end



# ╔═╡ aa31509b-b1e3-4e5e-8c34-54f89b6d6e30
test_data

# ╔═╡ 706a0c20-d42c-4fd6-980c-122b9c66de46
begin
	test_infer = downsample(test_data, 12)
	test_model = infer_test(test_infer, posterior)
	
	test_chain = sample(test_model, NUTS(0.65), MCMCSerial(), 100, 3; progress=false)
end

# ╔═╡ 4ceeb508-413e-4a3c-9790-bfec1900ef4d
test_posterior = DataFrame(test_chain)

# ╔═╡ 319ecc01-f2d6-46e8-87dd-9cb7992d544c
viz_fit(test_posterior, test_data, test_data)

# ╔═╡ e95b3140-7658-483d-bfc3-c86399dfd6a0
viz_posterior(test_posterior)

# ╔═╡ 2086cdce-516e-424a-a010-9433197e0699
A_top

# ╔═╡ 5c53d8b4-4929-47d9-ab18-aee0ec8e9efc
std(posterior.A_b)

# ╔═╡ 798d8d16-1c19-400d-8a94-e08c7f991e33
@model function infer_area(data_infer, posterior)
	N = nrow(data_infer)

	# get start area from posterior
	h0_obs = data_infer[1, "p_liquid level [cm]"]
	prior_area = linear_interpolation(slices, [mean(posterior.A_b), 
											  mean(posterior.A_t)]) 
	
	# Prior distributions.
	p_A ~ Normal(prior_area(h0_obs), std(posterior.A_t))
		# Normal(prior_area(h0_obs), std(posterior.A_t))
	
	W ~ filldist(Normal(0, 1), N)
	
	γ ~ Uniform(0.0, 10000.0)
	
	
	_a ~ TruncatedNormal(mean(posterior._a), std(posterior._a), 0.0, Inf) # area of the orifice [cm²]
	
	_c ~ TruncatedNormal(mean(posterior._c), std(posterior._c), 0.0, 1.0) # discharge coefficient
	
	_h_hole ~ TruncatedNormal(mean(posterior._h_hole), std(posterior._h_hole), 2.0, 3.0) # height of orifice from the base of tank
	
	σ ~ Uniform(0.0, 1.0) # measurement noise
	
	h_0 ~ TruncatedNormal(h0_obs, σ, h0_obs - 2 * σ, h0_obs +  2 * σ) # initial liquid level
	
	# for j in 2:N
	# 	p_A[j] = γ * W[j]  + p_A[j-1]
	# end

	_areas = [p_A]
	
	[append!(_areas, _areas[i-1] + γ * W[i]) for i in 2:N]

	idx = sortperm(data_infer[:, "p_liquid level [cm]"])
	
	_A = linear_interpolation(data_infer[idx, "p_liquid level [cm]"], _areas[idx], 
							  extrapolation_bc=Line())
	# print(W)
	# parameter for ODE
	p = (
			  a = _a, # area of the orifice [cm²]
			  g =  980.665, # cm/s², 
			  c = _c, # discharge coefficient
			  h_hole = _h_hole, # height of hole to the base of tank [cm] 
			  A = _A # A(h)
			)

	
	
	# # set up ODE
	_prob = ODEProblem(f, [h_0], tspan, p)

	sol = solve(_prob, Tsit5())
	
	# Observations.
	for i in 2:N
		tᵢ = data_infer[i, "Time [s]"]
		data_infer[i, "p_liquid level [cm]"] ~ Normal(sol(tᵢ)[1], σ)
	end
	return nothing
end



# ╔═╡ b4c62168-24e4-4cd3-8358-7599813af45d
tspan

# ╔═╡ bcd9167a-8fe3-4458-9afe-42d750719d35
A_top

# ╔═╡ 02939a87-e811-4ae4-8b6b-173370029889
begin
	area_model = infer_area(infer_data, posterior)
	area_posterior = DataFrame(sample(area_model, NUTS(0.65), MCMCSerial(), 100, 3; 
									  progress=false))
end

# ╔═╡ a127225a-5b79-4074-a16b-cecd11030800
function viz_area(original_post, posterior, data)
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="height [cm]", ylabel="Area [cm²]")
	
	ids_post = sample(1:nrow(posterior), 100; replace=false)

	areas = [mean(original_post[:, "A_b"]), mean(original_post[:, "A_t"])] # cm^2
	A = linear_interpolation(slices, areas)
	
	for j in ids_post
		
		_areas = [posterior[j, "p_A"]]
		[append!(_areas, _areas[i-1] + posterior[j, "γ"] * posterior[j, "W[$(i)]"]) for i in 2:nrow(data)]
		
		
		lines!(data[:, "p_liquid level [cm]"], _areas, label="model", color=(:green, 0.1))
	end

	scatter!(data[:, "p_liquid level [cm]"], A.(data[:, "p_liquid level [cm]"]), label="experimental")
	ylims!(A_bottom, A_top)
	axislegend(unique=true, position=:rb)
	return fig
end

# ╔═╡ b43f9f58-94fd-4c92-8e91-9a6b86cfc041
viz_area(posterior, area_posterior, infer_data)

# ╔═╡ 9c7e9485-2b7c-4847-bc7a-2049fbddf2cc
infer_data

# ╔═╡ Cell order:
# ╠═e1b54f1e-3c1b-4657-ae3d-d68a78617a4d
# ╠═faf59350-8d67-11ee-0bdd-2510e986118b
# ╠═4391f124-cbef-46e5-8462-e4e5126f5b38
# ╠═245836a9-6b44-4639-9209-e7ad9035e293
# ╟─418525b7-c358-41da-b865-5df3feb15855
# ╠═a95e371e-9319-4c7e-b5d9-4c4a50d12cd7
# ╠═23ee0e85-a84b-4b63-b432-5526559efcee
# ╟─7752316d-9dd0-4403-aa08-22c977ff3727
# ╠═20fa4266-be80-4d8e-b1d0-155a40a1241f
# ╠═9a7e5903-69be-4e0a-8514-3e05feedfed5
# ╟─078c01f7-e47e-4af0-be1c-ac4527b735fd
# ╠═9af216b7-4bf2-42fb-bd95-5b2040d019a7
# ╠═d3eeccff-39b3-429d-aec6-e7f1a500b729
# ╠═6e922ea4-01d0-4202-baed-3cb693b663bc
# ╠═759852ee-50e7-4deb-ac7e-c4693103c2a7
# ╠═b06a1c07-6250-4324-8802-010e5d847edb
# ╠═8b00d2b3-9182-42ab-8393-91707b813f60
# ╠═2a3d5a3f-5d20-4efb-91fd-cc84e77d5922
# ╠═f0675774-b9ba-40a9-953f-ce5784454378
# ╠═5100771a-38b6-437a-91e4-70495391067b
# ╠═8b6d766a-8f7b-4b9a-9a15-0f7375087120
# ╟─df589086-cc94-4da9-b982-563cd6ecc643
# ╠═33f889d7-e875-40d8-9d6d-cc87b0fbaf22
# ╠═e379461f-8896-4e2a-a71b-1871a8a37eb5
# ╠═c6a263eb-cb45-4ee7-9c02-549c89298652
# ╠═31306e0b-9748-48a8-b9d2-892cb501b7ba
# ╠═a934a309-6ab6-46e0-b448-2572334f5b7a
# ╠═d3307918-1fdb-4f87-bb92-67330d22e58b
# ╠═3853c2fb-3665-4c7b-8511-af01694f04b3
# ╠═444f6d74-273e-486d-905a-1443ec0e98df
# ╠═a1a10e2f-1b78-4b93-9295-7c0055e32692
# ╠═f21dc58e-d4e8-4314-b5dd-abbcb29efe86
# ╠═8f5b8859-6b8c-4f2a-af3a-b13c2d33fe2a
# ╠═8082559e-a5b0-41a8-b8ed-aec3b09e5b2b
# ╠═7ebe4680-c583-4f92-8bae-dd84c3fb5139
# ╠═a2048a65-7b7c-41fc-b7ee-f9199f3e96b5
# ╠═5bb0b72a-8c77-4fcb-bbde-d144986d9c1e
# ╠═2ab35999-3615-4f5c-8d89-36d77802fe9b
# ╠═5bea087e-c241-424d-a526-25eae30bfe15
# ╠═765ff940-1328-4806-aeaa-de8a41a6f4df
# ╠═2a01b228-f281-46c4-9764-fac6cc1b4217
# ╠═193d0e02-988e-4f58-b57d-7a1d125069a4
# ╠═aa31509b-b1e3-4e5e-8c34-54f89b6d6e30
# ╠═706a0c20-d42c-4fd6-980c-122b9c66de46
# ╠═4ceeb508-413e-4a3c-9790-bfec1900ef4d
# ╠═319ecc01-f2d6-46e8-87dd-9cb7992d544c
# ╠═e95b3140-7658-483d-bfc3-c86399dfd6a0
# ╠═2086cdce-516e-424a-a010-9433197e0699
# ╠═5c53d8b4-4929-47d9-ab18-aee0ec8e9efc
# ╠═798d8d16-1c19-400d-8a94-e08c7f991e33
# ╠═b4c62168-24e4-4cd3-8358-7599813af45d
# ╠═bcd9167a-8fe3-4458-9afe-42d750719d35
# ╠═02939a87-e811-4ae4-8b6b-173370029889
# ╠═a127225a-5b79-4074-a16b-cecd11030800
# ╠═b43f9f58-94fd-4c92-8e91-9a6b86cfc041
# ╠═9c7e9485-2b7c-4847-bc7a-2049fbddf2cc
