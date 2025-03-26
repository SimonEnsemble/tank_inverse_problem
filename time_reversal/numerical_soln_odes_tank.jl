### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ a855718e-9078-11ed-2fd1-73a706483e34
begin
	import Pkg; Pkg.activate()
	using DifferentialEquations, CairoMakie, DataFrames, MakieThemes, Statistics, Interpolations
	
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
H = 14.0 # cm

# ╔═╡ e6bbe642-f899-43ad-8015-097a65415fe4
md"top and bottom perimeter"

# ╔═╡ 29f37729-12c1-45de-a13b-a71793683c0b
P_b = 33.5 # cm

# ╔═╡ aba01fa8-57dc-4373-b674-7cafd02ef22d
P_t = 34.8 # cm

# ╔═╡ dd438ccc-343c-4421-95af-a17796dc95e8
md"top and bottom dimensions (approximate as a square)"

# ╔═╡ 29f5ce09-b4ab-4c46-bdfc-80f146ca3474
L_b = P_b / 4 # cm

# ╔═╡ 151c5229-12bb-439e-8a35-1506184deaf3
L_t = P_t / 4 # cm

# ╔═╡ ee537c0a-7cab-4bc2-89bd-3db590b51172
md"area of liquid, from a helicopter view, as a function of liquid level"

# ╔═╡ eeffeefd-5d09-487d-8c62-d5c06fc5355c
function A(h)
	θ = h / H
	# length of square here.
	L = θ * L_t + (1 - θ) * L_b
	return L ^ 2 # square
end

# ╔═╡ a7cc0cba-2046-4e12-92f9-5866982092d2
lines(
	range(0, H, 100), A.(range(0, H, 100)),
	axis=(; xlabel="h [cm]", ylabel="A [cm²]")
)

# ╔═╡ 99dcb519-7ccc-489f-8c53-97a1cfe5d90f
md"radius of orifice (drill bit size)"

# ╔═╡ dde9e79e-0f30-4698-9d8e-35f565931e05
rₒ = 3 / 32 * 2.54 / 2 # cm

# ╔═╡ bcba0331-ae6b-493b-8a9c-763eb4b07426
md"
## prediction 

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

# ╔═╡ 2ba54ac2-6b75-4e59-9ad5-5997d1b7555c
time_span = (0.0, 625.0) # s

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
h₀ = 14.0 # cm

# ╔═╡ 27e8ac1b-88bd-4d48-835f-f01f269a4fe9
function data_to_c(exp_data::DataFrame)
	id = 8
	Δh = exp_data[id, "h [cm]"] - exp_data[1, "h [cm]"] # cm
	Δt = exp_data[id, "t [s]"] - exp_data[1, "t [s]"] # s
	h̄ = mean(exp_data[1:id, "h [cm]"])
	c = - A(h̄) * Δh / Δt / (π * rₒ^2 * sqrt(2 * g * h̄))
	return c
end

# ╔═╡ 9e75496a-456c-4c54-9037-3ea0d951ad30
function Δts_to_data(Δts::Vector{Float64})
	n = length(Δts) + 1
	return DataFrame(
	    "h [cm]" => [h₀ - i for i = 0:n-1],
	    "t [s]" => [sum(Δts[1:i]) for i = 0:n-1]
	)
end

# ╔═╡ ca687109-e558-4329-9ba9-5045fc5ee455
exp_data = [Δts_to_data(Δts[r]) for r = 1:n_runs]

# ╔═╡ 8d2e68a5-d5d3-429e-896e-8aef48a7ab9e
cs = [data_to_c(exp_data[r]) for r = 1:n_runs]

# ╔═╡ ad774b05-ada2-4b75-bf9b-d76b9147b09c
c = mean(cs)

# ╔═╡ a1bd034e-b160-4687-8bea-fe7bf8b44c12
# we aren't going to use the parameter vector argument, nor time here explicitly.
function f(h, 🐸, t)
	if h < 0.0
		return 0.0
	else
		return - c * π * rₒ ^ 2 * sqrt(2 * g * h) / A(h)
	end
end

# ╔═╡ 4abecf6e-9404-4672-bdac-234e735ee817
f(h₀, [], 0.0) # initial rate of change

# ╔═╡ bc95d069-698e-4587-a03e-04c5ea27a9d7
# DifferentialEquations.jl syntax
prob = ODEProblem(f, h₀, time_span, saveat=0.1)

# ╔═╡ c5d7269b-ba22-4596-a0e6-e2266ad2ccc3
# solve ODE, return data frame
sim_data = DataFrame(solve(prob))

# ╔═╡ 12ffea5c-bf4a-4030-a44b-9f2bc2c64b6b
h_sim = linear_interpolation(sim_data[:, "timestamp"], sim_data[:, "value"])

# ╔═╡ 5c1d708e-1442-4a53-941c-5b2013dc34c7
begin
	fig = Figure()
	ax = Axis(
		fig[1, 1],
		xlabel="time, t [s]",
		ylabel="liquid height, h(t) [m]"
	)
	lines!(
		sim_data[:, "timestamp"], sim_data[:, "value"], 
		label="theory", color=Cycled(1)
	)
	
	axislegend()
	fig
end

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
			exp_data[r][:, "h [cm]"], exp_data[r][:, "h [cm]"] .- h_sim.(exp_data[r][:, "t [s]"]),
			marker=:rect, markersize=18, strokewidth=3,
			color=("white", 0.0), label="run $r", strokecolor=colors[r]
		)
	end
	fig
end

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
# ╠═eeffeefd-5d09-487d-8c62-d5c06fc5355c
# ╠═a7cc0cba-2046-4e12-92f9-5866982092d2
# ╟─99dcb519-7ccc-489f-8c53-97a1cfe5d90f
# ╠═dde9e79e-0f30-4698-9d8e-35f565931e05
# ╟─bcba0331-ae6b-493b-8a9c-763eb4b07426
# ╠═10c9de41-4bf6-447f-89a8-6b438390dd6d
# ╠═a1bd034e-b160-4687-8bea-fe7bf8b44c12
# ╠═4abecf6e-9404-4672-bdac-234e735ee817
# ╠═2ba54ac2-6b75-4e59-9ad5-5997d1b7555c
# ╠═bc95d069-698e-4587-a03e-04c5ea27a9d7
# ╠═c5d7269b-ba22-4596-a0e6-e2266ad2ccc3
# ╠═12ffea5c-bf4a-4030-a44b-9f2bc2c64b6b
# ╟─43389dd2-9237-4812-bfe5-a6570f552073
# ╠═5c1d708e-1442-4a53-941c-5b2013dc34c7
# ╟─0124ddb7-1bb6-44cd-a496-a02ada7b7131
# ╠═5e6d8bd0-3a2a-4b3c-b44e-7c6a2f8b9f4a
# ╠═05d6b493-7d6a-4414-94ac-9c80a5be220b
# ╠═3b57a44c-3c8c-4f8e-8d8e-e6288b8dc073
# ╠═22c9bf73-68fa-4b2c-8d7a-e3c63e084e10
# ╠═27e8ac1b-88bd-4d48-835f-f01f269a4fe9
# ╠═8d2e68a5-d5d3-429e-896e-8aef48a7ab9e
# ╠═ad774b05-ada2-4b75-bf9b-d76b9147b09c
# ╠═9e75496a-456c-4c54-9037-3ea0d951ad30
# ╠═ca687109-e558-4329-9ba9-5045fc5ee455
# ╠═ee3cf773-0615-4d8e-a9b1-30c22efb0bad
# ╠═8a7dfe16-d595-4c00-bfa6-3094749c7bed
