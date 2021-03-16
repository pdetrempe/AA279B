### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 8d7043ec-13f9-11eb-236c-7d3c970a700b
begin
	using LinearAlgebra
	using Plots
	using DifferentialEquations
	using Roots
	using OrdinaryDiffEq
	using BoundaryValueDiffEq
	using Pkg; Pkg.add("ParameterizedFunctions"); using ParameterizedFunctions
end

# ╔═╡ 96310354-13f9-11eb-031a-45054109c983
begin
	# Problem parameters
	μₑ = 398600 # [km^3/s^2]
	μₘ = 4903
	R = 384400  # Earth-Moon distance [km]

	R_Earth = 6378 # Earth radius [km]
	R_Moon = 1738  # Moon radius  [km]
	h = 200

	μ₁ = μₑ
	μ₂ = μₘ

	# Synodic reference frame parameters
	R₁ = μ₂/(μ₂+μ₁)*R        # Earth-CM distance
	R₂ = μ₁/(μ₂+μ₁)*R        # Moon-CM distance
	r1 = [-R₁; 0; 0]
	r2 = [ R₂; 0; 0]
	ωₛ = sqrt((μ₁ + μ₂)/(R^3))

	params = (μ₁,μ₂,r1,r2,ωₛ)
end

# ╔═╡ a711edbe-13f9-11eb-1e18-2581ad57050d
begin
	# Initial guess
	v_inj = 10.92367104     # [km/s]
	ϕ_inj = deg2rad(47.70061087)

	# Initial conditions (synodic frame)
	r0 = [-cos(ϕ_inj)*(h+R_Earth)-R₁; -sin(ϕ_inj)*(h+R_Earth); 0] # No z-position/velocity
	v0 = [v_inj*sin(ϕ_inj); -v_inj*cos(ϕ_inj); 0]
	x0 = [r0; v0]
end

# ╔═╡ b68789ca-13f9-11eb-07c1-8dd8e2ebae16
begin
	# 1a) Find L₁ point
	f₁(x) = -μ₁*(x-R₂)^2 + μ₂*(x+R₁)^2 + ωₛ^2*abs(x)*(x+R₁)^2*(x-R₂)^2
	x₁ = find_zero(f₁, .5*R)
	L₁ = [x₁; 0; 0]
end

# ╔═╡ d1180148-13f9-11eb-0f40-bf71bf49a8dc
# 1b) Use a fixed-time, single-shoot method to hit L₁
t_final = float(3*24*3600)


# ╔═╡ c4ac1dae-13f9-11eb-1761-dbc5d8dc3fe2
begin

# Let's play around with BVPs and Shooting Methods in DifferentialEquations
function CR3BP_EOM(t,x,dx,p)
	μ₁,μ₂,r1,r2,ωₛ = p
    ω_vec = [0;0;1]*ωₛ
    r3 = x[1:3]
    velocity = x[4:6]

    r13 = r3-r1
    r23 = r3-r2

    acceleration = -μ₁*r13/norm(r13)^3 -μ₂*r23/norm(r23)^3 -    # gravity
                    cross(ω_vec, cross(ω_vec, r3)) -            # centripetal
                    2*cross(ω_vec, velocity)                    # coriolis

    dx[1:3] = velocity
    dx[4:6] = acceleration
end

end

# ╔═╡ d75848d8-13f9-11eb-2b6f-2345809fc6d6
begin
	
function bc_L1(residual, u)
    residual[1] = u[end][1] - L₁[1]
    residual[2] = u[end][2] - L₁[2]
    residual[3] = u[end][3] - L₁[3]
end
	
end

# ╔═╡ ede9d6d4-13f9-11eb-0a9e-1b5b18254328
begin
	t_span = (0.0, t_final)
bvp = BVProblem(CR3BP_EOM, bc_L1, x0, t_span, params)
end

# ╔═╡ 46995714-13ff-11eb-3654-e9f82250cd35
sol = solve(bvp, Shooting(Vern7()))

# ╔═╡ 55e3eacc-13ff-11eb-0a0d-49c87b1b0f26
plot(sol.u[:,1], sol.u[:,2])

# ╔═╡ Cell order:
# ╠═8d7043ec-13f9-11eb-236c-7d3c970a700b
# ╠═96310354-13f9-11eb-031a-45054109c983
# ╠═a711edbe-13f9-11eb-1e18-2581ad57050d
# ╠═b68789ca-13f9-11eb-07c1-8dd8e2ebae16
# ╠═d1180148-13f9-11eb-0f40-bf71bf49a8dc
# ╠═c4ac1dae-13f9-11eb-1761-dbc5d8dc3fe2
# ╠═d75848d8-13f9-11eb-2b6f-2345809fc6d6
# ╠═ede9d6d4-13f9-11eb-0a9e-1b5b18254328
# ╠═46995714-13ff-11eb-3654-e9f82250cd35
# ╠═55e3eacc-13ff-11eb-0a0d-49c87b1b0f26
