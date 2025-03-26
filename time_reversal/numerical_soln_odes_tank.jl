### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a855718e-9078-11ed-2fd1-73a706483e34
begin
	import Pkg; Pkg.activate()
	using DifferentialEquations, CairoMakie, DataFrames, MakieThemes, Statistics, Interpolations, Turing
	
	# modifying the plot scheme
	# see here for other themes
	#  https://makieorg.github.io/MakieThemes.jl/dev/themes/ggthemr/
	set_theme!(ggthemr(:earth))
	update_theme!(fontsize=20, linewidth=4)
end

# ╔═╡ 1febf009-e353-4fba-90fc-016d1f14e7ba
colors = MakieThemes.GGThemr.ColorTheme[:earth][:swatch]

# ╔═╡ 27efd7e0-24a9-4c14-b515-42c9cd7f521d
md"🦫 take length-measurements to characterize the geometry of the tank."

# ╔═╡ 444dd8bb-9c2a-43a1-9b01-3d4465085ef3
md"height of tank"

# ╔═╡ 42b4f36a-d164-411b-867d-849bf0b8da26
H_obs = 14.0 # cm

# ╔═╡ e6bbe642-f899-43ad-8015-097a65415fe4
md"top and bottom perimeter"

# ╔═╡ 29f37729-12c1-45de-a13b-a71793683c0b
P_b_obs = 33.5 # cm

# ╔═╡ aba01fa8-57dc-4373-b674-7cafd02ef22d
P_t_obs = 34.8 # cm

# ╔═╡ dd438ccc-343c-4421-95af-a17796dc95e8
md"top and bottom dimensions (approximate as a square)"

# ╔═╡ 29f5ce09-b4ab-4c46-bdfc-80f146ca3474
L_b_obs = P_b_obs / 4 # cm

# ╔═╡ 151c5229-12bb-439e-8a35-1506184deaf3
L_t_obs = P_t_obs / 4 # cm

# ╔═╡ ee537c0a-7cab-4bc2-89bd-3db590b51172
md"area of liquid, from a helicopter view, as a function of liquid level"

# ╔═╡ 99dcb519-7ccc-489f-8c53-97a1cfe5d90f
md"radius of orifice (drill bit size)"

# ╔═╡ dde9e79e-0f30-4698-9d8e-35f565931e05
rₒ_obs = 3 / 32 * 2.54 / 2 # cm

# ╔═╡ bcba0331-ae6b-493b-8a9c-763eb4b07426
md"
# simulator 

predict the trajectory of the liquid level $h(t)$ using the model.

for `DifferentialEquations.jl`, we view the model as:

```math
\begin{equation}
\frac{dh}{dt}=f(h, \mathbf{p}, t) = -c \sqrt{h} / A(h)
\end{equation}
```
where $\mathbf{p}$ is an optional vector of parameters we do not use here.


!!! note
	we use `DifferentialEquations.jl` (documentation [here](https://docs.sciml.ai/DiffEqDocs/stable/)) to numerically solve (well, approximate the solution to) nonlinear differential equations.
"

# ╔═╡ 10c9de41-4bf6-447f-89a8-6b438390dd6d
g = 9.8 * 100 # cm / s²

# ╔═╡ 3424875e-c3e9-4a33-a8f2-ec6d3ab1bedd
function sim_h(
	# initial height
	h₀,
	# tank geometry
	P_t, P_b, H, # cm
	# discharge coefficient
	c,
	# orifice radius
	rₒ,
	# simulation time
	tf=625.0 # s
)
	# assume square, calculate lengths from perimeter
	L_t = P_t / 4.0
	L_b = P_b / 4.0

	# area of tank as a function of h
	function A(h)
		θ = h / H
		# length of square here.
		L = θ * L_t + (1 - θ) * L_b
		return L ^ 2 # square
	end

	# RHS of ODE
	function f!(dh, h, 🐸, t)
		if h[1] < 0.0
			dh[1] = 0.0
			return 0.0
		else
			return dh[1] = - c * π * rₒ ^ 2 * sqrt(2 * g * h[1]) / A(h[1])
		end
	end

	prob = ODEProblem(f!, [h₀], (0.0, tf))

	sim_data = DataFrame(
		solve(
			prob, 
			Tsit5(), saveat=0.5, 
			reltol=1e-5, abstol=1e-5
		)
	)

	h_of_t = linear_interpolation(sim_data[:, "timestamp"], sim_data[:, "value1"])
	
	return sim_data, h_of_t
end

# ╔═╡ 43389dd2-9237-4812-bfe5-a6570f552073
md"🦫 viz the solution.
"

# ╔═╡ 0124ddb7-1bb6-44cd-a496-a02ada7b7131
md"## comparison to experiment

🦫 measure time series data $\{(t_i, h_i)\}$ and plot against the simulation.

> The Unreasonable Effectiveness of Mathematics in the Natural Sciences

this was done in Cory's kitchen.
"

# ╔═╡ 5e6d8bd0-3a2a-4b3c-b44e-7c6a2f8b9f4a
begin
	Δts = Dict{Int64, Vector{Float64}}()
	# run #1 Cory's kitchen 3/25
	Δts[1] = [
		21.01,
		24.82,
		24.82,
		26.42,
		27.98,
		29.82,
		32.40,
		34.57,
		38.60,
		41.44,
		52.14,
		60+3.99,
	]

	Δts[2] = [
		21.81,
		23.12, 
		25.15,
		25.33,
		27.77,
		29.25,
		31.98,
		33.57,
		37.52,
		40.27,
		54.12,
		60+2.51
	]

	Δts[3] = [
		21.87,
		23.19,
		25.25,
		25.92,
		26.72,
		28.73,
		31.42,
		32.33,
		37.93,
		41.17,
		48.92,
		60.69,
	]
end

# ╔═╡ 05d6b493-7d6a-4414-94ac-9c80a5be220b
n_runs = length(Δts)

# ╔═╡ 3b57a44c-3c8c-4f8e-8d8e-e6288b8dc073
[length(Δts[r]) for r = 1:n_runs]

# ╔═╡ 22c9bf73-68fa-4b2c-8d7a-e3c63e084e10
h₀_obs = 14.0 # cm

# ╔═╡ 931a34aa-25cd-48e0-968c-d8830c2d241e
sim_data, h_of_t = sim_h(
	h₀_obs,
	P_t_obs, P_b_obs, H_obs,
	0.45,
	rₒ_obs
)

# ╔═╡ 5c1d708e-1442-4a53-941c-5b2013dc34c7
begin
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel="time, t [s]",
		ylabel="liquid height, h(t) [m]"
	)
	lines!(
		sim_data[:, "timestamp"], sim_data[:, "value1"], 
		label="theory", color=Cycled(1)
	)
	
	axislegend()
	fig
end

# ╔═╡ 74f20bc5-e6ea-4a2b-8f4c-9b07a9256f84
function Δts_to_data(Δts::Vector{Float64})
	n = length(Δts) + 1
	return DataFrame(
	    "h [cm]" => [h₀_obs - i for i = 0:n-1],
	    "t [s]" => [sum(Δts[1:i]) for i = 0:n-1]
	)
end

# ╔═╡ 27e8ac1b-88bd-4d48-835f-f01f269a4fe9
function data_to_c(exp_data::DataFrame)
	id = 8
	Δh = exp_data[id, "h [cm]"] - exp_data[1, "h [cm]"] # cm
	Δt = exp_data[id, "t [s]"] - exp_data[1, "t [s]"] # s
	h̄ = mean(exp_data[1:id, "h [cm]"])
	c = - A(h̄) * Δh / Δt / (π * rₒ^2 * sqrt(2 * g * h̄))
	return c
end

# ╔═╡ ca687109-e558-4329-9ba9-5045fc5ee455
exp_data = [Δts_to_data(Δts[r]) for r = 1:n_runs]

# ╔═╡ ee3cf773-0615-4d8e-a9b1-30c22efb0bad
begin
	for r = 1:n_runs
		scatter!(
			ax,
			exp_data[r][:, "t [s]"], exp_data[r][:, "h [cm]"],
			marker=:rect, markersize=18, strokewidth=3,
			color=("white", 0.0), label="run $r", strokecolor=colors[r]
		)
	end
	axislegend(ax)
	fig
end

# ╔═╡ 8a7dfe16-d595-4c00-bfa6-3094749c7bed
begin
	local fig = Figure()
	local ax = Axis(fig[1, 1], xlabel="h [m]", ylabel="residual")
	for r = 1:n_runs
		scatter!(
			ax,
			exp_data[r][:, "h [cm]"], exp_data[r][:, "h [cm]"] .- h_of_t.(exp_data[r][:, "t [s]"]),
			marker=:rect, markersize=18, strokewidth=3,
			color=("white", 0.0), label="run $r", strokecolor=colors[r]
		)
	end
	fig
end

# ╔═╡ ab5ac22a-97ae-4014-ad13-01faa608c82b
md"# Bayesian inference of $c$"

# ╔═╡ 8a45adb0-c3a8-472b-89d2-72ad42e5e92c
@model function forward_model_ID_c(
	data::DataFrame, prior_only::Bool=false
)
	#=
	prior distributions
	=#
	# uncertainty for length measurements via tape
	σ_ℓ = 0.2 # cm
	# uncertainty for liquid level measurements
	σ_h ~ Uniform(0.0, 0.5) # cm

	# tank geo
	P_t ~ Normal(P_t_obs, σ_ℓ)
	P_b ~ Normal(P_b_obs, σ_ℓ)
	H ~ Normal(H_obs, σ_ℓ)

	# discharge coefficient
	c ~ Normal(0.65, 0.4)

	# orifice radius
	rₒ ~ Normal(rₒ_obs, rₒ_obs * 0.1)
	
	# initial liquid level
	h₀ ~ Normal(h₀_obs, σ_h)

	# for prior, do not show the algo the data :)
	if prior_only
		return nothing
	end
	
	# set up, solve ODE
	sim_data, h_of_t = sim_h(
		h₀,
		P_t, P_b, H,
		c,
		rₒ
	)

	# observations.
	for i in 2:nrow(data)
		tᵢ = data[i, "t [s]"]
		ĥᵢ = h_of_t(tᵢ)
		data[i, "h [cm]"] ~ Normal(ĥᵢ, σ_h)
	end
	
	return nothing
end

# ╔═╡ 5de14c78-6c0a-4da9-8df0-e5859d666b8e
vcat(exp_data[1], exp_data[2])

# ╔═╡ b45eddf1-8fdc-424d-b655-70764d13515c
begin
	train_model = forward_model_ID_c(
		# exp_data[1]
		# vcat(exp_data[1], exp_data[2])
		vcat(exp_data[1], exp_data[2], exp_data[3])
	)
	
	train_posterior = DataFrame(
		sample(
			train_model, 
			NUTS(0.65), MCMCSerial(), 50, 2; progress=true
		)
	)
end

# ╔═╡ af0c1f1c-7732-494d-8780-98a7652d0c75
hist(
	train_posterior[:, "c"],
	axis=(; xlabel="c", ylabel="# samples")
)

# ╔═╡ c8ee4a00-89a1-4480-a447-0a5759f398f2
mean(train_posterior[:, "σ_h"])

# ╔═╡ c5737e97-33c2-401c-941e-63d4bf784542
mean(train_posterior[:, "c"])

# ╔═╡ 81373b01-bb57-4641-839c-6808e4feee02
std(train_posterior[:, "c"])

# ╔═╡ 43e439e8-2657-44ea-8641-b2ca75445cbe
train_posterior

# ╔═╡ a03e6a2c-d0fc-4e4f-9477-677ee967921a
function posterior_predictive_check(posterior::DataFrame, exp_data::DataFrame)
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel="time, t [s]",
		ylabel="liquid height, h(t) [m]"
	)
	
	n_sample = 25
	for s = sample(1:nrow(posterior), n_sample)
		sim_data, h_of_t = sim_h(
			posterior[s, :h₀],
			posterior[s, :P_t], posterior[s, :P_b], posterior[s, :H],
			posterior[s, :c],
			posterior[s, :rₒ]
		)

		lines!(
			sim_data[:, "timestamp"], sim_data[:, "value1"], 
			color=(colors[2], 0.2)
		)
	end
	scatter!(
		exp_data[:, "t [s]"], exp_data[:, "h [cm]"],
		marker=:rect, markersize=18, strokewidth=3,
		color=("white", 0.0), strokecolor="white"
	)
	fig
end

# ╔═╡ ecc138dd-fc9a-465f-bf38-0ec681966c09
posterior_predictive_check(train_posterior, exp_data[1])

# ╔═╡ 6eab9bd4-1a04-403a-ba61-fc9c511305bf
posterior_predictive_check(train_posterior, exp_data[2])

# ╔═╡ 78a7aec1-b75a-48c7-bf11-0124556540e1
posterior_predictive_check(train_posterior, exp_data[3])

# ╔═╡ 72ec118e-02b4-4b2d-9498-9d0db08e0556
md"# time reversal I: what was $h_0$?"

# ╔═╡ b60e64f9-e8d8-49f4-8381-fdeb4c3060c2


# ╔═╡ ea4b93dd-4eee-4d49-8a04-6da7cbc30d78
md"# time reversal I: what was $t_0$?"

# ╔═╡ Cell order:
# ╠═a855718e-9078-11ed-2fd1-73a706483e34
# ╠═1febf009-e353-4fba-90fc-016d1f14e7ba
# ╟─27efd7e0-24a9-4c14-b515-42c9cd7f521d
# ╟─444dd8bb-9c2a-43a1-9b01-3d4465085ef3
# ╠═42b4f36a-d164-411b-867d-849bf0b8da26
# ╟─e6bbe642-f899-43ad-8015-097a65415fe4
# ╠═29f37729-12c1-45de-a13b-a71793683c0b
# ╠═aba01fa8-57dc-4373-b674-7cafd02ef22d
# ╟─dd438ccc-343c-4421-95af-a17796dc95e8
# ╠═29f5ce09-b4ab-4c46-bdfc-80f146ca3474
# ╠═151c5229-12bb-439e-8a35-1506184deaf3
# ╟─ee537c0a-7cab-4bc2-89bd-3db590b51172
# ╟─99dcb519-7ccc-489f-8c53-97a1cfe5d90f
# ╠═dde9e79e-0f30-4698-9d8e-35f565931e05
# ╟─bcba0331-ae6b-493b-8a9c-763eb4b07426
# ╠═10c9de41-4bf6-447f-89a8-6b438390dd6d
# ╠═3424875e-c3e9-4a33-a8f2-ec6d3ab1bedd
# ╠═931a34aa-25cd-48e0-968c-d8830c2d241e
# ╟─43389dd2-9237-4812-bfe5-a6570f552073
# ╠═5c1d708e-1442-4a53-941c-5b2013dc34c7
# ╟─0124ddb7-1bb6-44cd-a496-a02ada7b7131
# ╠═5e6d8bd0-3a2a-4b3c-b44e-7c6a2f8b9f4a
# ╠═05d6b493-7d6a-4414-94ac-9c80a5be220b
# ╠═3b57a44c-3c8c-4f8e-8d8e-e6288b8dc073
# ╠═22c9bf73-68fa-4b2c-8d7a-e3c63e084e10
# ╠═74f20bc5-e6ea-4a2b-8f4c-9b07a9256f84
# ╠═27e8ac1b-88bd-4d48-835f-f01f269a4fe9
# ╠═ca687109-e558-4329-9ba9-5045fc5ee455
# ╠═ee3cf773-0615-4d8e-a9b1-30c22efb0bad
# ╠═8a7dfe16-d595-4c00-bfa6-3094749c7bed
# ╟─ab5ac22a-97ae-4014-ad13-01faa608c82b
# ╠═8a45adb0-c3a8-472b-89d2-72ad42e5e92c
# ╠═5de14c78-6c0a-4da9-8df0-e5859d666b8e
# ╠═b45eddf1-8fdc-424d-b655-70764d13515c
# ╠═af0c1f1c-7732-494d-8780-98a7652d0c75
# ╠═c8ee4a00-89a1-4480-a447-0a5759f398f2
# ╠═c5737e97-33c2-401c-941e-63d4bf784542
# ╠═81373b01-bb57-4641-839c-6808e4feee02
# ╠═43e439e8-2657-44ea-8641-b2ca75445cbe
# ╠═a03e6a2c-d0fc-4e4f-9477-677ee967921a
# ╠═ecc138dd-fc9a-465f-bf38-0ec681966c09
# ╠═6eab9bd4-1a04-403a-ba61-fc9c511305bf
# ╠═78a7aec1-b75a-48c7-bf11-0124556540e1
# ╟─72ec118e-02b4-4b2d-9498-9d0db08e0556
# ╠═b60e64f9-e8d8-49f4-8381-fdeb4c3060c2
# ╟─ea4b93dd-4eee-4d49-8a04-6da7cbc30d78
