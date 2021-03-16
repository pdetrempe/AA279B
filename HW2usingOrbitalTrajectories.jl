### A Pluto.jl notebook ###
# v0.12.4

# using Markdown, InteractiveUtils, LinearAlgebra,
using OrbitalTrajectories, Plots, DifferentialEquations, Roots, Unitful
using ForwardDiff # For calculating jacobian of injection parameter function
# using OrdinaryDiffEq
# using BoundaryValueDiffEq
# using Pkg; Pkg.add("ParameterizedFunctions"); using ParameterizedFunctions

# Problem parameters
μₑ = 398600u"km^3/s^2"# [km^3/s^2]
μₘ = 4903u"km^3/s^2"
R = 384400u"km" # Earth-Moon distance [km]

R_Earth = 6378u"km" # Earth radius [km]
R_Moon = 1738u"km"# Moon radius  [km]
h = 200u"km"

μ₁ = μₑ
μ₂ = μₘ

# Synodic reference frame parameters
R₁ = μ₂/(μ₂+μ₁)*R       # Earth-CM distance
R₂ = μ₁/(μ₂+μ₁)*R       # Moon-CM distance
r1 = [-R₁; 0; 0]
r2 = [ R₂; 0; 0]
ωₛ = sqrt((μ₁ + μ₂)/(R^3))

params = (μ₁,μ₂,r1,r2,ωₛ)
# Initial guess
v_inj = 10.92367104u"km/s"    # [km/s]
ϕ_inj = deg2rad(47.70061087)u"rad"
x_inj = [ustrip(v_inj); ustrip(ϕ_inj)]

# Initial conditions (synodic frame)


function injectionVelandAngleToCartesian(x)
    v_inj = x[1]
    ϕ_inj = x[2]
    r0 = [-cos(ϕ_inj)*ustrip((h+R_Earth)-R₁); -sin(ϕ_inj)*ustrip((h+R_Earth)); 0] # No z-position/velocity
    v0 = [v_inj*sin(ϕ_inj); -v_inj*cos(ϕ_inj); 0]
    x0 = [r0; v0]
end

x0 = injectionVelandAngleToCartesian(x_inj)

# # Symmetric initial conditions
# ϕ_inj = 0
# v_inj = 10.75
# r0 = [-cos(ϕ_inj)*(h+R_Earth)-R₁; -sin(ϕ_inj)*(h+R_Earth); 0] # No z-position/velocity
# v0 = [v_inj*sin(ϕ_inj); -v_inj*cos(ϕ_inj); 0]
# x0 = [r0; v0]

# 1a) Find L₁ point
f₁(x) = -μ₁*(x-R₂)^2 + μ₂*(x+R₁)^2 + ωₛ^2*abs(x)*(x+R₁)^2*(x-R₂)^2
x₁ = find_zero(f₁, .5*R)
L₁ = [x₁; 0; 0]
println(L₁)

# 1b) Use a fixed-time, single-shoot method to hit L₁
t_final = float(3*24*3600)

EarthMoonSystem = CR3BP(:earth, :moon)

# Compute non-dimensional initial state:
# Function for converting initial system state from position + energy to position + velocity
function compute_initial_state(system, x0_in, v0_in; tmax=365u"d")
    r0 = x0_in ./ system.props.L # Use non-dimensional length
    v0 = v0_in ./ system.props.V # Use non-dimensional velocity
    # ^ Probably issue here, figure out how to non-dimensionalize (or just enter Jacobi constant)

    # v0 = √(2*centrifugal_potential(system.props.μ, [x0, 0., 0.]) - J)
    State(system, [r0' v0'], (0, tmax / system.props.T)) # Use non-dimensional time
end

# Convert initial states
IC = compute_initial_state(EarthMoonSystem, x0[1:3]*u"km", x0[4:6]*u"km/s", tmax=t_final*u"s")

# Solve problem and propagate uncorrected solution
traj_uncorrected = solve(IC)
p1 = plot(traj_uncorrected)

# # Create free return traj_uncorrected using axisymmetric correcotr
# symmetric_solution = solve(IC, DiffCorrectAxisymmetric())
# plot!(p1, symmetric_solution; color=:blue)


# Create and use generic fixed-time single-shoot method to hit L₁ point by varying v_ing and ϕ_inj (not vₓ and v\_y)
# Build the initial state and propagate it
# tspan = [0 t_final]./EarthMoonSystem.props.T # Non-dimensionalize time
# state = State(EarthMoonSystem, SynodicFrame(), IC, tspan)
# trajectory = solve(state)
# plot(trajectory; legend=:bottomleft, case.plot_args...)

# Use automatic differentiation to calculate the STM at the final point of the trajectory (w.r.t initial state)
Φ = sensitivity(AD, IC) # STM, dimensionless

# Calculate Jacobian of x, y, ẋ and ẏ with respect to ϕ_inj and v_inj
g = x -> ForwardDiff.jacobian(injectionVelandAngleToCartesian, x)
x_inj_stripped = [ustrip(x_inj[1]); ustrip(x_inj[2])]
dXdInjectionParams = g(x_inj_stripped) # Input sensitivity matrix, HAS DIMENSIONS!!!

dXf_dInjectionParams = Φ * dXdInjectionParams



# # Let's play around with BVPs and Shooting Methods in DifferentialEquations
# function CR3BP_EOM(t,x,dx,p)
# 	μ₁,μ₂,r1,r2,ωₛ = p
#     ω_vec = [0;0;1]*ωₛ
#     r3 = x[1:3]
#     velocity = x[4:6]

#     r13 = r3-r1
#     r23 = r3-r2

#     acceleration = -μ₁*r13/norm(r13)^3 -μ₂*r23/norm(r23)^3 -    # gravity
#                     cross(ω_vec, cross(ω_vec, r3)) -            # centripetal
#                     2*cross(ω_vec, velocity)                    # coriolis

#     dx[1:3] = velocity
#     dx[4:6] = acceleration
# end

# function bc_L1(residual, u)
#     residual[1] = u[end][1] - L₁[1]
#     residual[2] = u[end][2] - L₁[2]
#     residual[3] = u[end][3] - L₁[3]
# end

# t_span = (0.0, t_final)
# bvp = BVProblem(CR3BP_EOM, bc_L1, x0, t_span, params)

# sol = solve(bvp, Shooting(Vern7()))
# plot(sol.u[:,1], sol.u[:,2])


