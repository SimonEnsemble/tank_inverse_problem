### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 3f3c91b4-d6da-4252-ae8d-91795e0aadd0
begin
	# top
	L_t = 14.6 #cm
	W_t = 9.0 # cm
	p_t = 44.3 # cm
	
	# bottom
	L_b = 13.4 # cm
	W_b = 7.8 # cm
	p_b = 40.1 # cm

	# slant
	H★ = 28.6 # cm
	H = sqrt(H★ ^ 2 - ((L_t - L_b) / 2)^2) # negligible really.
end

# ╔═╡ 693534df-c9c8-4b67-81a6-162027cead65
L_t * W_t

# ╔═╡ 7d689656-18eb-43bd-bdf7-b823f544f0b4
L_b * W_b

# ╔═╡ c172ee18-7809-4264-9508-ce12e82584db
r_b = (p_b / 2 - (L_b + W_b)) / (π - 4)

# ╔═╡ ae081a16-93e5-4227-8bb4-e570f08f4b0d
r_t = (p_t / 2 - (L_t + W_t)) / (π - 4)

# ╔═╡ c15c2b78-3a30-45bc-8874-a9df5c33857b
A_b = (L_b - 2 * r_b) * (W_b - 2 * r_b) + 2 * r_b * (L_b + W_b - 4 * r_b) + π * r_b ^ 2

# ╔═╡ 81893df6-7fc9-4c3e-898b-a3dfabfa5d82
A_t = (L_t - 2 * r_t) * (W_t - 2 * r_t) + 2 * r_t * (L_t + W_t - 4 * r_t) + π * r_t ^ 2

# ╔═╡ baa26362-5f31-434c-b0f7-2caf3df8b83c
A(h) = h / H * A_t + (1 - h / H) * A_b

# ╔═╡ Cell order:
# ╠═3f3c91b4-d6da-4252-ae8d-91795e0aadd0
# ╠═693534df-c9c8-4b67-81a6-162027cead65
# ╠═7d689656-18eb-43bd-bdf7-b823f544f0b4
# ╠═baa26362-5f31-434c-b0f7-2caf3df8b83c
# ╠═c172ee18-7809-4264-9508-ce12e82584db
# ╠═ae081a16-93e5-4227-8bb4-e570f08f4b0d
# ╠═c15c2b78-3a30-45bc-8874-a9df5c33857b
# ╠═81893df6-7fc9-4c3e-898b-a3dfabfa5d82
