### A Pluto.jl notebook ###
# v0.19.38

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
	r_hole = 1 / 8 * 2.54 / 2 # cm
	a_hole = π * r_hole ^ 2 # cm²

	h_hole = 2.5 # cm
end

# ╔═╡ 633af570-2e83-41b8-8854-fcaad776b1de
function A(h::Float64)
	if h < 0 || h > H
		error("tank over- or under-flow")
	end
	return h / H * A_t + (1 - h / H) * A_b
end

# ╔═╡ e884b0ef-f5b1-44c8-8f55-35f8480cadc1
md"## calibration of liquid sensor"

# ╔═╡ e602eab3-b4a1-41d5-9181-4f9c253a5677
calibration_data = CSV.read("../liquid calibration.csv", DataFrame)

# ╔═╡ 3a2a1f79-bdae-411b-8f7e-62b3edcacb2e
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="level sensor reading", 
		ylabel="liquid level (cm)"
	)
	lines!(
		calibration_data[:, "level sensor reading"], 
		calibration_data[:, "liquid level (cm)"]
	)
	fig
end

# ╔═╡ ad0d77c2-2c9d-4d6d-849a-52a95343d8ab
sensor_output_to_h = linear_interpolation(
	calibration_data[:, "level sensor reading"], 
	calibration_data[:, "liquid level (cm)"]
)

# ╔═╡ 7c804847-d4fa-4ad3-8063-59f0c793150c
sensor_output_to_h(522.0)

# ╔═╡ 180dfe22-6948-431f-a7eb-8e85e4b0f6df
md"## process data"

# ╔═╡ 27c461b1-9ded-417a-83a8-d08643e72e32
function read_process_data(expt_no::Int; t_begin::Float64=0.0)
	# read in
	data = CSV.read(
		"../experiment $expt_no- 1-12.csv", 
		DataFrame, 
		types=[Float64, Float64]
	)
	# downsample
	data = data[1:100:nrow(data), :]

	# shift time
	filter!(row -> row["Time [s]"] > t_begin, data)
	data[:, "Time [s]"] = data[:, "Time [s]"] .- t_begin

	# apply calibration curve
	data[:, "h [cm]"] = map(sensor_output_to_h, data[:, " liquid_level_reading"])
	
	return data
end

# ╔═╡ 9e760ae6-bd5b-4fdc-a385-133b4869e163
data = read_process_data(3, t_begin=50.0)

# ╔═╡ 4d603ba3-2b2d-4c37-a7af-88e333c32a6c
function viz_data(data::DataFrame)
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		xlabel="time [s]", 
		ylabel="h [cm]"
	)
	scatter!(
		data[:, "Time [s]"], 
		data[:, "h [cm]"],
		label="data"
	)
	fig
end

# ╔═╡ b375842f-f762-4ac4-a645-19d0f479b0c3
viz_data(data)

# ╔═╡ bb175715-e92c-46db-83a6-1941ff7b5815
md"## ODE model"

# ╔═╡ e8de8919-4147-4255-afbf-de9f13e5e37a
g = 980.0 # cm / s²

# ╔═╡ afaa06c2-da85-11ee-1846-cb19e3e8c03b
function f(h, params, t)
	c = params[1]
	
	if h < h_hole
		return 0.0
	end
	return - c * a_hole * sqrt(2 * g * (h - h_hole)) / A(h)
end

# ╔═╡ 9bdca4d1-9c1d-4f7f-a30d-bf5bfc571028
h₀ = data[1, "h [cm]"]

# ╔═╡ 7afa6f66-4e6e-4d45-b038-19a26ffb605d
tspan = (0.0, 410.0) # s

# ╔═╡ dce788c9-ed23-453e-afe4-76b384799ee1
prob = ODEProblem(f, h₀, tspan, [1.0], saveat=1.0)

# ╔═╡ d3e4cf56-22de-44c3-9d65-4c46cc2f6392
sim_data = DataFrame(solve(prob, Tsit5()))

# ╔═╡ eb2afb65-f7f1-49e4-b5d0-b41366c4e350
begin
	fig = viz_data(data)
	ax = current_axis(fig)
	lines!(ax, sim_data[:, "timestamp"], sim_data[:, "value"], label="sim", color=Cycled(2))
	axislegend()
	fig
end

# ╔═╡ 66f9b070-a677-45ce-ae2e-c67b4b816243
md"## tune c"

# ╔═╡ b846feb2-f072-4b5d-bd0c-ced6f5e91849
function loss(c)
	# define ODE problem
	prob = ODEProblem(f, h₀, tspan, [c], saveat=1.0)
	sim_h = solve(prob, Tsit5())

	# compute loss
	ℓ = 0.0
	for row in eachrow(data)
		# extract data point
		tᵢ = row["Time [s]"]
		hᵢ = row["h [cm]"]
		
		# get predicted h
		ĥᵢ = sim_h(tᵢ)

		# increment loss
		ℓ += (hᵢ - ĥᵢ) ^ 2
	end
	return ℓ
end

# ╔═╡ c8e6165e-7bca-4784-bed9-e655fb33644d
cs = range(0.7, 1.1, length=15)

# ╔═╡ f9860b59-9943-4eac-b7aa-d7a06340ad79
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1], xlabel="c", ylabel="loss")
	lines!(cs, loss.(cs))
	fig
end

# ╔═╡ 0de278a0-0f4d-416c-800a-f7790c02e343
c_opt = optimize(loss, 0.5, 1.0).minimizer

# ╔═╡ 92068db5-faa1-4d15-bdba-e1d0bf9ddc93
md"## plot optimal model"

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
		label="sim",
		color=Cycled(2)
	)
	axislegend()
	fig
end

# ╔═╡ 5acc0622-68ed-48f8-ad5d-56225fa51afe
md"### test data"

# ╔═╡ d8fe78f3-ab8b-4716-bea2-47edcad2d912
data_test = read_process_data(4, t_begin=50.0)

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
		label="sim",
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
# ╠═3a2a1f79-bdae-411b-8f7e-62b3edcacb2e
# ╠═ad0d77c2-2c9d-4d6d-849a-52a95343d8ab
# ╠═7c804847-d4fa-4ad3-8063-59f0c793150c
# ╟─180dfe22-6948-431f-a7eb-8e85e4b0f6df
# ╠═27c461b1-9ded-417a-83a8-d08643e72e32
# ╠═9e760ae6-bd5b-4fdc-a385-133b4869e163
# ╠═4d603ba3-2b2d-4c37-a7af-88e333c32a6c
# ╠═b375842f-f762-4ac4-a645-19d0f479b0c3
# ╟─bb175715-e92c-46db-83a6-1941ff7b5815
# ╠═e8de8919-4147-4255-afbf-de9f13e5e37a
# ╠═afaa06c2-da85-11ee-1846-cb19e3e8c03b
# ╠═9bdca4d1-9c1d-4f7f-a30d-bf5bfc571028
# ╠═7afa6f66-4e6e-4d45-b038-19a26ffb605d
# ╠═dce788c9-ed23-453e-afe4-76b384799ee1
# ╠═d3e4cf56-22de-44c3-9d65-4c46cc2f6392
# ╠═eb2afb65-f7f1-49e4-b5d0-b41366c4e350
# ╟─66f9b070-a677-45ce-ae2e-c67b4b816243
# ╠═b846feb2-f072-4b5d-bd0c-ced6f5e91849
# ╠═c8e6165e-7bca-4784-bed9-e655fb33644d
# ╠═f9860b59-9943-4eac-b7aa-d7a06340ad79
# ╠═0de278a0-0f4d-416c-800a-f7790c02e343
# ╟─92068db5-faa1-4d15-bdba-e1d0bf9ddc93
# ╠═f30f8862-81fc-4ae6-a218-fd9df8e60093
# ╠═06804f34-1bc3-43e9-8c52-bd5953ce1628
# ╠═17199b2d-fd12-4460-8581-219681ad6553
# ╟─5acc0622-68ed-48f8-ad5d-56225fa51afe
# ╠═d8fe78f3-ab8b-4716-bea2-47edcad2d912
# ╠═2ea68858-158b-48e9-b8eb-5230542f9f04
# ╠═ddc8d1ad-a5a6-46a2-b52f-acfeeb128823
# ╠═c21066d3-6fd7-45f7-ba6a-ccbad32c1b29
# ╠═d795458a-10ab-4ea5-9403-8fbe6804511d
# ╠═e07dfbad-5f2b-4e5a-aa25-17f224c01c9b
