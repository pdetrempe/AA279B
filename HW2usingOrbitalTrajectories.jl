### A Pluto.jl notebook ###
# v0.12.4

# using Markdown, InteractiveUtils, LinearAlgebra,
using OrbitalTrajectories, Plots, DifferentialEquations, Roots, Unitful
using ForwardDiff # For calculating jacobian of injection parameter function
using LinearAlgebra
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

function injectionVelandAngleToCartesian(x)
    v_inj = x[1]
    ϕ_inj = x[2]
    r0 = [-cos(ϕ_inj)*ustrip((h+R_Earth))-ustrip(R₁); -sin(ϕ_inj)*ustrip((h+R_Earth)); 0] # No z-position/velocity
    v0 = [v_inj*sin(ϕ_inj); -v_inj*cos(ϕ_inj); 0]
    x0 = [r0; v0]
end

# 1a) Find L₁ point
f₁(x) = -μ₁*(x-R₂)^2 + μ₂*(x+R₁)^2 + ωₛ^2*abs(x)*(x+R₁)^2*(x-R₂)^2
x₁ = find_zero(f₁, .5*R)
L₁ = [x₁; 0u"km"; 0u"km"]
println(L₁)

# 1b) Use a fixed-time, single-shoot method to hit L₁
t_final = float(3*24*3600)

EarthMoonSystem = CR3BP(:earth, :moon)

# Compute non-dimensional initial state:
# Function for converting initial system state from position + energy to position + velocity
function compute_initial_state(system, x0_in, v0_in; tmax=365u"d")
    r0 = x0_in ./ system.props.L # Use non-dimensional length
    v0 = v0_in ./ system.props.V # Use non-dimensional velocity
    State(system, [r0' v0'], (0, tmax / system.props.T)) # Use non-dimensional time
end



# Create and use generic fixed-time single-shoot method to hit L₁ point by varying v_ing and ϕ_inj (not vₓ and v\_y)
# Build the initial state and propagate it
# tspan = [0 t_final]./EarthMoonSystem.props.T # Non-dimensionalize time
# state = State(EarthMoonSystem, SynodicFrame(), IC, tspan)
# trajectory = solve(state)
# plot(trajectory; legend=:bottomleft, case.plot_args...)


# F is our final position constraint (using Pavlak's nomenclature)
x_desired = L₁
F = 1000u"km"
ϵ_F = 1u"km" # Miss tolerance
max_iter = 10
num_iters = 0


global p1 = plot()

function fixedTimeSingleShoot(x_inj, x_desired, t_final, traj_uncorrected; system=CR3BP(:earth, :moon))
    F = 1000u"km"
    ϵ_F = 1u"km" # Miss tolerance
    max_iter = 10
    num_iters = 0

    while norm(F) > ϵ_F && num_iters < max_iter
        x0 = injectionVelandAngleToCartesian(x_inj)

        # Convert initial states
        IC = compute_initial_state(system, x0[1:3]*u"km", x0[4:6]*u"km/s", tmax=t_final*u"s")

        # Solve problem and propagate uncorrected solution
        traj_uncorrected = solve(IC)
        # plot!(p1,traj_uncorrected)

        # Use automatic differentiation to calculate the STM at the final point of the trajectory (w.r.t initial state)
        Φ = sensitivity(AD, IC) # STM, dimensionless

        # Calculate Jacobian of x, y, ẋ and ẏ with respect to ϕ_inj and v_inj
        g = x -> ForwardDiff.jacobian(injectionVelandAngleToCartesian, x)
        x_inj_stripped = [ustrip(x_inj[1]); ustrip(x_inj[2])]
        dXdInjectionParams = g(x_inj_stripped) # Input sensitivity matrix, Need to dimensionalize?
        # dXdInjectionParams[:,1] .= ustrip(dXdInjectionParams[:,1].*system.props.V) # Try dimensionalizing in place

        # Dimensionalize STM
        Φ[1:3, 4:6] .= ustrip( Φ[1:3, 4:6] * system.props.L ./ system.props.V )
        Φ[4:6, 1:3] .= ustrip( Φ[4:6, 1:3] * system.props.V ./ system.props.L )

        dXf_dInjectionParams = Φ * dXdInjectionParams

        # Using Pavlak's notation for derivative of final state w.r.t initial conditions
        DF = ustrip( dXf_dInjectionParams[1:3, :] )# variation of final state w.r.t. initial parameters
        F = traj_uncorrected.u[end][1:3]*system.props.L - L₁# Redimensionalize
        # println("F = ", F)
        # println("x_inj = ", x_inj)
        # println("pinv(DF) = ", pinv(DF))

        params_new = x_inj - pinv(DF) * ustrip(F)
        v_inj_new = params_new[1]u"km/s"
        ϕ_inj_new = params_new[2]u"rad"
        x_inj = [ustrip(v_inj_new); ustrip(ϕ_inj_new)]

        num_iters += 1
    end
    return x_inj, traj_uncorrected
end

function variableTimeSingleShoot(x_inj, x_desired, t_final_0, traj_uncorrected; system=CR3BP(:earth, :moon))
    F = 1000u"km"
    ϵ_F = 1u"km" # Miss tolerance
    max_iter = 10
    num_iters = 0

    while norm(F) > ϵ_F && num_iters < max_iter
        x0 = injectionVelandAngleToCartesian(x_inj)

        # Convert initial states
        IC = compute_initial_state(system, x0[1:3]*u"km", x0[4:6]*u"km/s", tmax=t_final_0*u"s")

        # Solve problem and propagate uncorrected solution
        traj_uncorrected = solve(IC)
        # plot!(p1,traj_uncorrected)

        # Use automatic differentiation to calculate the STM at the final point of the trajectory (w.r.t initial state)
        Φ = sensitivity(AD, IC) # STM, dimensionless

        # Calculate Jacobian of x, y, ẋ and ẏ with respect to ϕ_inj and v_inj
        g = x -> ForwardDiff.jacobian(injectionVelandAngleToCartesian, x)
        x_inj_stripped = [ustrip(x_inj[1]); ustrip(x_inj[2])]
        dXdInjectionParams = g(x_inj_stripped) # Input sensitivity matrix, Need to dimensionalize?
        # dXdInjectionParams[:,1] .= ustrip(dXdInjectionParams[:,1].*system.props.V) # Try dimensionalizing in place

        # Dimensionalize STM
        Φ[1:3, 4:6] .= ustrip( Φ[1:3, 4:6] * system.props.L ./ system.props.V )
        Φ[4:6, 1:3] .= ustrip( Φ[4:6, 1:3] * system.props.V ./ system.props.L )

        dXf_dInjectionParams = Φ * dXdInjectionParams

        # Using Pavlak's notation for derivative of final state w.r.t initial conditions
        DF = [ustrip( dXf_dInjectionParams[1:3, :] ) ustrip(traj_uncorrected.u[end][4:6]*system.props.V) ]# variation of final state w.r.t. initial parameters
        F = traj_uncorrected.u[end][1:3]*system.props.L - L₁# Redimensionalize


        params_new = [x_inj; t_final_0] - pinv(DF) * ustrip(F)
        v_inj_new = params_new[1]u"km/s"
        ϕ_inj_new = params_new[2]u"rad"
        t_final_0 = params_new[3]
        x_inj = [ustrip(v_inj_new); ustrip(ϕ_inj_new)]

        num_iters += 1
    end
    return x_inj, t_final_0, traj_uncorrected
end

# Initial conditions
x0 = injectionVelandAngleToCartesian(x_inj)
IC = compute_initial_state(EarthMoonSystem, x0[1:3]*u"km", x0[4:6]*u"km/s", tmax=t_final*u"s")

# Solve problem and propagate uncorrected solution
traj_uncorrected = solve(IC)
# p1 = plot(traj_uncorrected)

x_inj_fixed, traj_fixed = fixedTimeSingleShoot(x_inj, x_desired, t_final, traj_uncorrected; system=EarthMoonSystem)
x_inj_variable, t_flight_variable, traj_variable = variableTimeSingleShoot(x_inj, x_desired, t_final, traj_uncorrected; system=EarthMoonSystem)

println("Fixed state", x_inj_fixed)
println("Variable state/time", x_inj_variable, t_flight_variable)
plot!(p1, traj_fixed)
plot!(p1, traj_variable)

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


