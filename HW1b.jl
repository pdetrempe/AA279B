using LinearAlgebra
using Plots
using DifferentialEquations

μₑ = 398600 # [km^3/s^2]
μₘ = 4903
R = 384400 # Earth-Moon distance [km]

R_Earth = 6378 # Earth radius [km]
R_Moon = 1738  # Moon radius  [km]
h = 200

μ₁ = μₑ
μ₂ = μₘ

# Synodic reference frame parameters
R₁ = μ₂/(μ₂+μ₁)*R   # Earth-CM distance
R₂ = μ₁/(μ₂+μ₁)*R        # Moon-CM distance

# Initial conditions
v_inj = 10.92367104     # [km/s]
ϕ_inj = deg2rad(47.70061087)

# Equations of motion
r1 = [-R₁; 0; 0]
r2 = [ R₂; 0; 0]
ωₛ = [0;0;sqrt((μ₁ + μ₂)/(R^3))]

function CR3BP_EOM(x, p, t)
    μ₁,μ₂,r1,r2,ωₛ = p
    r3 = x[1:3]
    velocity = x[4:6]

    r13 = r3-r1
    r23 = r3-r2

    acceleration = -μ₁*r13/norm(r13)^3 -μ₂*r23/norm(r23)^3 -    # gravity
                    cross(ωₛ, cross(ωₛ, r3)) -                  # centripetal
                    2*cross(ωₛ, velocity)                       # coriolis
                   # ^gravity, centripetal, coriolis

    x_dot = [velocity; acceleration]
end


# Initial conditions (synodic)
r0 = [-cos(ϕ_inj)*(h+R_Earth)-R₁; -sin(ϕ_inj)*(h+R_Earth); 0] # No z-position/velocity
v0 = [v_inj*sin(ϕ_inj); -v_inj*cos(ϕ_inj); 0]
x0 = [r0; v0]

tspan = (0.,5.8*24*3600)

# Integrate trajectory
params = (μ₁,μ₂,r1,r2,ωₛ)
free_return = ODEProblem(CR3BP_EOM,x0,tspan,params)
sol = solve(free_return,alg_hints=[:stiff],reltol=1e-8, VCABM())

# Plot trajectory in synodic frame
gr()
plot(sol,vars=(1,2),denseplot=false)

# 1d) Plot spacecraft trajectory in inertial coordinates
θ = ωₛ[3].*sol.t
pos_synodic = [pos[1:3] for pos in sol.u]
function pos_synodicToInertial(pos,θ)
    [cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1]*pos
end
state_inert = pos_synodicToInertial.(pos_synodic, θ)

# positions = zeros(3, length(θ))
# for (idx, pos) in enumerate(state_inert)
#     positions[:, idx] = pos
# end
# plot(pos[1,:], pos[2,:])
