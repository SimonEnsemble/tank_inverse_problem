### A Pluto.jl notebook ###
# v0.20.1

using Markdown
using InteractiveUtils

# ╔═╡ 1f8b0daa-ff3f-4bf3-b2f4-caba7e5b56cd
begin
	import Pkg; Pkg.activate()
	using CSV, DataFrames, CairoMakie, OrdinaryDiffEq, PlutoUI, Interpolations, Optim
end

# ╔═╡ acee1559-0ac7-483e-9745-a7c55c9425af
begin
	import AlgebraOfGraphics
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
		resolution=resolution
	)
end

# ╔═╡ 39c4cec8-270a-4157-923f-663219e568d1
TableOfContents()

# ╔═╡ b5192de4-2149-4a47-8ab1-29fcf3420591
md"## tank geometry"

# ╔═╡ 26438937-e8a8-4523-8fa3-cf0637de6720
A_t = 128.95 # cm²

# ╔═╡ 037639b6-81e3-4e9b-bca0-9559049a38e6
A_b = 102.97 # cm²

# ╔═╡ 30a0346e-9d92-4632-95cd-556e661275c3
H = 28.6 # cm

# ╔═╡ ad79ab7a-74a6-4043-a1c1-7d4738e2d1ed
begin
	r_hole = 5 / 64 * 2.54 / 2 # cm
	a_hole = π * r_hole ^ 2 # cm²

	h_hole = 0.9 # cm
end

# ╔═╡ 633af570-2e83-41b8-8854-fcaad776b1de
function A(h::Float64)
	if h < 0 || h > H
		error("tank over- or under-flow")
	end
	return h / H * A_t + (1 - h / H) * A_b
end

# ╔═╡ e884b0ef-f5b1-44c8-8f55-35f8480cadc1
md"## calibration of liquid sensor

read in data for calibration of the liquid level sensor.
"

# ╔═╡ e602eab3-b4a1-41d5-9181-4f9c253a5677
calibration_data = CSV.read("../calibration_curve.csv", DataFrame)

# ╔═╡ 03a52d5f-2931-42be-9270-e15ac322057a
md"create a linear interpolator to obtain liquid level at level sensor readings between those recorded in the calibration data set."

# ╔═╡ ad0d77c2-2c9d-4d6d-849a-52a95343d8ab
sensor_output_to_h = linear_interpolation(
	calibration_data[:, "level sensor reading"], 
	calibration_data[:, "h [cm]"]
)

# ╔═╡ e73ae703-5814-44b1-b546-f983299a591e
sensor_output_to_h(522.0)

# ╔═╡ da2df80d-313f-4a8c-ac2b-483eea98c142
md"visualize the calibration curve and check the linear interpolator."

# ╔═╡ 3a2a1f79-bdae-411b-8f7e-62b3edcacb2e
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="level sensor reading", 
		ylabel="liquid level (cm)"
	)
	
	scatter!(
		calibration_data[:, "level sensor reading"], 
		calibration_data[:, "h [cm]"]
	)
	
	lr = range(
		minimum(calibration_data[:, "level sensor reading"]),
		maximum(calibration_data[:, "level sensor reading"]),
		length=150
	)
	lines!(lr, sensor_output_to_h.(lr))
	
	fig
end

# ╔═╡ 180dfe22-6948-431f-a7eb-8e85e4b0f6df
md"## process data

read in and process time series data collected during an experiment where the tank is emptying.
"

# ╔═╡ 27c461b1-9ded-417a-83a8-d08643e72e32
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
	data[:, "t [min]"] = data[:, "t [s]"] / 60

	#=
	apply calibration curve to get h
	=#
	data[:, "h [cm]"] = sensor_output_to_h.(data[:, "liquid_level_reading"])
	
	return select(data, ["t [min]", "h [cm]"])
end

# ╔═╡ b5edda70-4771-48f8-a11a-afcb2b4bf33d
function downsample(data::DataFrame, n::Int)
	# filter out sensor data that is below 508
	# this is the liquid level where the jet shooting outside the tank
	# dies, and it just runs down the tank (Toricelli's law invalid)
	data = filter("h [cm]" => (x -> x >= sensor_output_to_h(508)), data)
	
	ids = trunc.(Int, collect(range(1, nrow(data), length=n)))
	return data[ids, :]
end

# ╔═╡ bd4f9b9f-7af4-44cf-8b54-65222d64de13
train_experiment = "../no_obs_4_18_2.csv"

# ╔═╡ 176a95b8-c2e9-4713-ad51-8997a7b65c09
test_experiment  = "../no_obs_4_18_3.csv"

# ╔═╡ b2757732-d93b-479c-a20b-a8e8aef75efb
data = downsample(read_h_time_series(train_experiment), 20)

# ╔═╡ 4d603ba3-2b2d-4c37-a7af-88e333c32a6c
function viz_data(data::DataFrame)
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="time [min]", 
		ylabel="h [cm]"
	)
	scatter!(
		data[:, "t [min]"], 
		data[:, "h [cm]"],
		label="data"
	)
	xlims!(0, nothing)
	ylims!(0, nothing)
	fig
end

# ╔═╡ b375842f-f762-4ac4-a645-19d0f479b0c3
viz_data(data)

# ╔═╡ bb175715-e92c-46db-83a6-1941ff7b5815
md"## ODE model
set up ODE model governing $h(t;c)$ where $c$ is a paramter.
"

# ╔═╡ e8de8919-4147-4255-afbf-de9f13e5e37a
g = 980.665 * 60^2 # cm/s² -> cm/min²

# ╔═╡ afaa06c2-da85-11ee-1846-cb19e3e8c03b
# right-hand side of ODE
function f(h, params, t)
	c = params[1]
	
	if h < h_hole
		return 0.0
	end
	return - c * a_hole * sqrt(2 * g * (h - h_hole)) / A(h)
end

# ╔═╡ 9bdca4d1-9c1d-4f7f-a30d-bf5bfc571028
h₀ = data[1, "h [cm]"] # initial liquid level

# ╔═╡ 7afa6f66-4e6e-4d45-b038-19a26ffb605d
tspan = (0.0, 20.0) # min

# ╔═╡ dce788c9-ed23-453e-afe4-76b384799ee1
prob = ODEProblem(f, h₀, tspan, [1.0], saveat=1.0)

# ╔═╡ d3e4cf56-22de-44c3-9d65-4c46cc2f6392
sim_data = DataFrame(solve(prob, Tsit5()))

# ╔═╡ db2a6ec1-c424-4839-93c3-967ff3fff973
md"assess quality of fit with $c=1$"

# ╔═╡ eb2afb65-f7f1-49e4-b5d0-b41366c4e350
begin
	fig = viz_data(data)
	ax = current_axis(fig)
	lines!(ax, sim_data[:, "timestamp"], sim_data[:, "value"], label="sim", color=Cycled(2))
	axislegend()
	fig
end

# ╔═╡ 66f9b070-a677-45ce-ae2e-c67b4b816243
md"## tune c

let's use the time series data to find the optimal $c$ parameter.

first, we code up a loss function
"

# ╔═╡ b846feb2-f072-4b5d-bd0c-ced6f5e91849
function loss(c)
	# use ODE solver to create h(t; c)
	prob = ODEProblem(f, h₀, tspan, [c], saveat=1.0)
	h = solve(prob, Tsit5())

	# compute loss
	ℓ = 0.0
	for row in eachrow(data[1:end-3, :])
		# extract data point
		tᵢ = row["t [min]"]
		hᵢ = row["h [cm]"]
		
		# get predicted h
		ĥᵢ = h(tᵢ)

		# increment loss
		ℓ += (hᵢ - ĥᵢ) ^ 2
	end
	return ℓ
end

# ╔═╡ c8e6165e-7bca-4784-bed9-e655fb33644d
cs = range(0.75, 1.0, length=15)

# ╔═╡ f9860b59-9943-4eac-b7aa-d7a06340ad79
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1], xlabel="c", ylabel="loss")
	lines!(cs, loss.(cs))
	fig
end

# ╔═╡ a321075d-4da7-4123-a03f-82b28cd56046
md"find optimal $c$ by minimizing the loss."

# ╔═╡ 0de278a0-0f4d-416c-800a-f7790c02e343
c_opt = optimize(loss, 0.4, 1.0).minimizer

# ╔═╡ 92068db5-faa1-4d15-bdba-e1d0bf9ddc93
md"## assess fit of optimal model

re-solve ODE with optimal $c$. plot against data to assess fit.
"

# ╔═╡ f30f8862-81fc-4ae6-a218-fd9df8e60093
prob_opt = ODEProblem(f, h₀, tspan, [c_opt], saveat=1.0)

# ╔═╡ 06804f34-1bc3-43e9-8c52-bd5953ce1628
sim_data_opt = DataFrame(solve(prob_opt, Tsit5()))

# ╔═╡ 17199b2d-fd12-4460-8581-219681ad6553
begin
	local fig = viz_data(data)
	local ax = current_axis(fig)
	lines!(
		ax, 
		sim_data_opt[:, "timestamp"], 
		sim_data_opt[:, "value"], 
		label="simulation",
		color=Cycled(2)
	)
	xlims!(0, nothing)
	ylims!(0, nothing)
	axislegend()
	save("train_fit.pdf", fig)
	fig
end

# ╔═╡ 947b8285-36c6-4949-9d3b-7b8da0f19e2c
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="time [min]", 
		ylabel="h [cm]"
	)
	lines!(
		ax, 
		sim_data_opt[:, "timestamp"], 
		sim_data_opt[:, "value"],
		color=Cycled(2)
	)
	scatter!([0], [h₀])
	xlims!(0, tspan[2])
	ylims!(0, nothing)
	save("just_sim.pdf", fig)
	fig
end

# ╔═╡ 5acc0622-68ed-48f8-ad5d-56225fa51afe
md"### test model

the ultimate judge of the quality of the model is how well it predicts the outcome of a different experiment.

the test data contains time series data for $h(t)$ when the tank was draining in a second experiment.
"

# ╔═╡ d8fe78f3-ab8b-4716-bea2-47edcad2d912
data_test = downsample(read_h_time_series(test_experiment), 20)

# ╔═╡ 2ea68858-158b-48e9-b8eb-5230542f9f04
viz_data(data_test)

# ╔═╡ ddc8d1ad-a5a6-46a2-b52f-acfeeb128823
h₀_test = data_test[1, "h [cm]"]

# ╔═╡ c21066d3-6fd7-45f7-ba6a-ccbad32c1b29
prob_test = ODEProblem(f, h₀_test, tspan, [c_opt], saveat=1.0)

# ╔═╡ d795458a-10ab-4ea5-9403-8fbe6804511d
sim_data_test = DataFrame(solve(prob_test, Tsit5()))

# ╔═╡ e07dfbad-5f2b-4e5a-aa25-17f224c01c9b
begin
	local fig = viz_data(data_test)
	local ax = current_axis(fig)
	lines!(
		ax, 
		sim_data_test[:, "timestamp"], 
		sim_data_test[:, "value"], 
		label="simulation",
		color=Cycled(2)
	)
	axislegend()
	fig
end

# ╔═╡ Cell order:
# ╠═1f8b0daa-ff3f-4bf3-b2f4-caba7e5b56cd
# ╠═acee1559-0ac7-483e-9745-a7c55c9425af
# ╠═39c4cec8-270a-4157-923f-663219e568d1
# ╟─b5192de4-2149-4a47-8ab1-29fcf3420591
# ╠═26438937-e8a8-4523-8fa3-cf0637de6720
# ╠═037639b6-81e3-4e9b-bca0-9559049a38e6
# ╠═30a0346e-9d92-4632-95cd-556e661275c3
# ╠═ad79ab7a-74a6-4043-a1c1-7d4738e2d1ed
# ╠═633af570-2e83-41b8-8854-fcaad776b1de
# ╟─e884b0ef-f5b1-44c8-8f55-35f8480cadc1
# ╠═e602eab3-b4a1-41d5-9181-4f9c253a5677
# ╟─03a52d5f-2931-42be-9270-e15ac322057a
# ╠═ad0d77c2-2c9d-4d6d-849a-52a95343d8ab
# ╠═e73ae703-5814-44b1-b546-f983299a591e
# ╟─da2df80d-313f-4a8c-ac2b-483eea98c142
# ╠═3a2a1f79-bdae-411b-8f7e-62b3edcacb2e
# ╟─180dfe22-6948-431f-a7eb-8e85e4b0f6df
# ╠═27c461b1-9ded-417a-83a8-d08643e72e32
# ╠═b5edda70-4771-48f8-a11a-afcb2b4bf33d
# ╠═bd4f9b9f-7af4-44cf-8b54-65222d64de13
# ╠═176a95b8-c2e9-4713-ad51-8997a7b65c09
# ╠═b2757732-d93b-479c-a20b-a8e8aef75efb
# ╠═4d603ba3-2b2d-4c37-a7af-88e333c32a6c
# ╠═b375842f-f762-4ac4-a645-19d0f479b0c3
# ╟─bb175715-e92c-46db-83a6-1941ff7b5815
# ╠═e8de8919-4147-4255-afbf-de9f13e5e37a
# ╠═afaa06c2-da85-11ee-1846-cb19e3e8c03b
# ╠═9bdca4d1-9c1d-4f7f-a30d-bf5bfc571028
# ╠═7afa6f66-4e6e-4d45-b038-19a26ffb605d
# ╠═dce788c9-ed23-453e-afe4-76b384799ee1
# ╠═d3e4cf56-22de-44c3-9d65-4c46cc2f6392
# ╟─db2a6ec1-c424-4839-93c3-967ff3fff973
# ╠═eb2afb65-f7f1-49e4-b5d0-b41366c4e350
# ╟─66f9b070-a677-45ce-ae2e-c67b4b816243
# ╠═b846feb2-f072-4b5d-bd0c-ced6f5e91849
# ╠═c8e6165e-7bca-4784-bed9-e655fb33644d
# ╠═f9860b59-9943-4eac-b7aa-d7a06340ad79
# ╟─a321075d-4da7-4123-a03f-82b28cd56046
# ╠═0de278a0-0f4d-416c-800a-f7790c02e343
# ╟─92068db5-faa1-4d15-bdba-e1d0bf9ddc93
# ╠═f30f8862-81fc-4ae6-a218-fd9df8e60093
# ╠═06804f34-1bc3-43e9-8c52-bd5953ce1628
# ╠═17199b2d-fd12-4460-8581-219681ad6553
# ╠═947b8285-36c6-4949-9d3b-7b8da0f19e2c
# ╟─5acc0622-68ed-48f8-ad5d-56225fa51afe
# ╠═d8fe78f3-ab8b-4716-bea2-47edcad2d912
# ╠═2ea68858-158b-48e9-b8eb-5230542f9f04
# ╠═ddc8d1ad-a5a6-46a2-b52f-acfeeb128823
# ╠═c21066d3-6fd7-45f7-ba6a-ccbad32c1b29
# ╠═d795458a-10ab-4ea5-9403-8fbe6804511d
# ╠═e07dfbad-5f2b-4e5a-aa25-17f224c01c9b
