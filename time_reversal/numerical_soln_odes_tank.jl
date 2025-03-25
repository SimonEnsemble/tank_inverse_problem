### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ a855718e-9078-11ed-2fd1-73a706483e34
begin
	import Pkg; Pkg.activate()
	using DifferentialEquations, CairoMakie, DataFrames, MakieThemes
	
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
rₒ = 5 / 64 * 2.54 / 2 # cm

# ╔═╡ 3815c34f-60b5-4b03-8959-0cc8e822481a
md"## parameter identification
our dynamic model is:

$$A(h)\dfrac{dh}{dt}=-c\sqrt{h(t)}$$

with the initial condition $h(t)=h_0$.

conduct a quick experiment to identify the $c$ parameter.
"

# ╔═╡ 22c9bf73-68fa-4b2c-8d7a-e3c63e084e10
h₀ = 14.0 # cm

# ╔═╡ 988b38f2-56be-4ab7-8e96-6c89ef523c33
h₁ = 13.0 # cm

# ╔═╡ c01e7bd2-9d31-4adc-95bb-75b64608e9df
Δh = h₁ - h₀ # cm

# ╔═╡ 03b5e719-ecff-47fc-8d6e-6c257cce898c
# mean height during this experiment
h̄ = (h₀ + h₁) / 2 # cm

# ╔═╡ 34a53346-101c-4ed8-b0aa-63a5088f8622
# c = - A(h̄) * Δh / Δt / (π * rₒ^2 * sqrt(2 * g * h̄))

# ╔═╡ 10c9de41-4bf6-447f-89a8-6b438390dd6d
g = 9.8 * 100 # cm / s²

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

# ╔═╡ 2ba54ac2-6b75-4e59-9ad5-5997d1b7555c
time_span = (0.0, 625.0) # s

# ╔═╡ bc95d069-698e-4587-a03e-04c5ea27a9d7
# DifferentialEquations.jl syntax
prob = ODEProblem(f, h₀, time_span, saveat=0.1)

# ╔═╡ c5d7269b-ba22-4596-a0e6-e2266ad2ccc3
# solve ODE, return data frame
sim_data = DataFrame(solve(prob))

# ╔═╡ 43389dd2-9237-4812-bfe5-a6570f552073
md"🦫 viz the solution.
"

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

# ╔═╡ 0124ddb7-1bb6-44cd-a496-a02ada7b7131
md"## comparison to experiment

🦫 measure time series data $\{(t_i, h_i)\}$ and plot against the simulation.

> The Unreasonable Effectiveness of Mathematics in the Natural Sciences

this was done in Cory's kitchen.
"

# ╔═╡ 5e6d8bd0-3a2a-4b3c-b44e-7c6a2f8b9f4a
# times to decrease by one cm in height.
# ts = [0, 22.38, 46.14, 70.03,
# 96.84,
# 124.67,
# 154.83,
# 186.57,
# 220.46,
# 261.27,
# 307.74,
# 368.28,
# 453.17,
# 616.05]

# ts=[0,
# 27.86,
# 58.17,
# 88.55,
# 121.33,
# 155.44,
# 192.83,
# 231.5,
# 272.99,
# 318.95,
# 378.12,
# 447.59,
# 537.39,
# 706.54]

ts = [
	0,
27.63,
57.98,
88.81,
121.89,
156.76,
193.84,
236.84,
278.52,
328.46,
383.4,
447.92,
536.81,
696.53
]

# ╔═╡ 7dd762b9-8505-4349-9c21-89b3d6117f60
Δt = ts[2] # s

# ╔═╡ e48cc019-a9e5-4a39-85df-81742d0fbb81
n = length(ts)

# ╔═╡ ca687109-e558-4329-9ba9-5045fc5ee455
exp_data = DataFrame(
	"h [cm]" => [h₀ - i for i = 0:n-1],
	"t [s]" => ts
)

# ╔═╡ ee3cf773-0615-4d8e-a9b1-30c22efb0bad
begin
	scatter!(
		ax,
		exp_data[:, "t [s]"], exp_data[:, "h [cm]"], label="data",
		marker=:rect, markersize=18, strokewidth=3,
		color=("white", 0.0), strokecolor=colors[1]
	)
	fig
end

# ╔═╡ 2152b8ee-a44c-48d4-8914-78ced2bbce83
c = 0.5

# ╔═╡ 8f8c2f8c-ac4e-44f3-ac7b-d28b3bdfd157
# ╠═╡ disabled = true
#=╠═╡
c = 0.5
  ╠═╡ =#

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
# ╟─3815c34f-60b5-4b03-8959-0cc8e822481a
# ╠═22c9bf73-68fa-4b2c-8d7a-e3c63e084e10
# ╠═988b38f2-56be-4ab7-8e96-6c89ef523c33
# ╠═c01e7bd2-9d31-4adc-95bb-75b64608e9df
# ╠═7dd762b9-8505-4349-9c21-89b3d6117f60
# ╠═03b5e719-ecff-47fc-8d6e-6c257cce898c
# ╠═34a53346-101c-4ed8-b0aa-63a5088f8622
# ╠═2152b8ee-a44c-48d4-8914-78ced2bbce83
# ╠═8f8c2f8c-ac4e-44f3-ac7b-d28b3bdfd157
# ╠═10c9de41-4bf6-447f-89a8-6b438390dd6d
# ╟─bcba0331-ae6b-493b-8a9c-763eb4b07426
# ╠═a1bd034e-b160-4687-8bea-fe7bf8b44c12
# ╠═4abecf6e-9404-4672-bdac-234e735ee817
# ╠═2ba54ac2-6b75-4e59-9ad5-5997d1b7555c
# ╠═bc95d069-698e-4587-a03e-04c5ea27a9d7
# ╠═c5d7269b-ba22-4596-a0e6-e2266ad2ccc3
# ╟─43389dd2-9237-4812-bfe5-a6570f552073
# ╠═5c1d708e-1442-4a53-941c-5b2013dc34c7
# ╟─0124ddb7-1bb6-44cd-a496-a02ada7b7131
# ╠═5e6d8bd0-3a2a-4b3c-b44e-7c6a2f8b9f4a
# ╠═e48cc019-a9e5-4a39-85df-81742d0fbb81
# ╠═ca687109-e558-4329-9ba9-5045fc5ee455
# ╠═ee3cf773-0615-4d8e-a9b1-30c22efb0bad
