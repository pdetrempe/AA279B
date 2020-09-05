using Roots
using LinearAlgebra
using Plots
using DifferentialEquations

μₑ = 398600 # [km^3/s^2]
μₘ = 4902
R = 384400 # Earth-Moon distance [km]

μ₁ = μₑ
μ₂ = μₘ

# Synodic reference frame parameters
R₁ = μ₂/(μ₂+μ₁)*R   # Earth-CM distance
R₂ = R - R₁         # Moon-CM distance
ωₛ = sqrt((μ₁ + μ₂)/R^3)

# L₄ and L₅ Lagrange points
L₄ = [-R₁ + R/2; √3/2*R; 0]
L₅ = [-R₁ + R/2; -√3/2*R; 0]
print(L₄)
print(L₅)

# Find L₁, L₂, and  L₃
f₁(x) = -μ₁*(x-R₂)^2 + μ₂*(x+R₁)^2 + ωₛ^2*abs(x)*(x+R₁)^2*(x-R₂)^2
f₂(x) = -μ₁*(x-R₂)^2 - μ₂*(x+R₁)^2 + ωₛ^2*abs(x)*(x+R₁)^2*(x-R₂)^2
f₃(x) =  μ₁*(x-R₂)^2 + μ₂*(x+R₁)^2 - ωₛ^2*abs(x)*(x+R₁)^2*(x-R₂)^2

x₁ = find_zero(f₁, .5*R)
x₂ = find_zero(f₂, 1.5*R)
x₃ = find_zero(f₃, -R)

L₁ = [x₁; 0; 0]
L₂ = [x₂; 0; 0]
L₃ = [x₃; 0; 0]
print(L₁, L₂, L₃)


# Calculate effective potential energy and plot contours
r₁ = [-R₁; 0; 0]
r₂ = [ R₂; 0; 0]
PE(r₃) = -μ₁/norm(r₃-r₁) -μ₂/norm(r₃-r₂) - .5*(ωₛ^2)*(r₃[1]^2 + r₃[2]^2)

num_points = 1000
points = [range(-1.25*R, stop = -.01*R, length = num_points); range(.01*R, stop = 1.25*R, length = num_points)]

PE_eff = [PE([point_x; point_y; 0]) for point_y in points, point_x in points]
p1 = contour(points, points, PE_eff, fill=true, levels=-1.68:.01:-1.56, xlabel="X", ylabel="Y")

# Free return trajectory
