using Roots
using LinearAlgebra
using Plots

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

# Find L₁ and L₂
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
PE(r₁, r₂, r₃) = -μ₁/norm(r₃-r₁) -μ₂/norm(r₃-r₂);# - .5*(ωₛ^2)*(r₃[1]^2 + r₃[2]^2)

r_grid = range(.25*R, stop = 2*R, length = 100)
θ_grid = range(0, stop = 2*pi, length = 100)
x_grid = [r*sin(θ) for r in r_grid, θ in θ_grid]
y_grid = [r*cos(θ) for r in r_grid, θ in θ_grid]

PE_eff = [PE(r₁, r₂, [x; y; 0]) for x in x_grid, y in y_grid]

p1 = contour(x_grid, y_grid, PE_eff)
plot(p1)
