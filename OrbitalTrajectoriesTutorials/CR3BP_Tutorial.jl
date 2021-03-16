using OrbitalTrajectories, Plots, DifferentialEquations, Roots, Unitful


# Set up Sun/Mars system
SunMars = CR3BP(:sun, :mars)


# Orbit parameters from Russel paper
# X in CR3BP space and J is the Jacobi constant
all_params = Dict(
    "g2" => [
        (x0 = 0.172, J = 3.000034),  # updated to get a similar trajectory
        (x0 = 0.432, J = 3.000118),
        (x0 = 0.556, J = 3.000203),
        (x0 = 0.775, J = 3.000203),
        (x0 = 0.841, J = 3.0001804), # highly sensitive, updated to similar
    ],
    "g1" => [
        (x0 = 0.006, J = 3.011712, tmax=1u"d"), # Units in days
        (x0 = 0.198, J = 3.000410),
        (x0 = 0.161, J = 3.000203),
        (x0 = 0.081, J = 3.000203),
        (x0 = 0.027, J = 3.000196), # highly sensitive, updated to similar
    ],
    "DRO" => [
        (x0 = -0.008, J = 3.009513, tmax=1u"d"),
        (x0 = -0.256, J = 3.000255),
        (x0 = -0.944, J = 3.000015),
    ]
)

# Function for converting initial system state from position + energy to position + velocity
function compute_initial_state(system; x0, J, tmax=365u"d")
    x0 = 1 - system.props.μ + (x0*(10^6)u"km" / system.props.L) # Use non-dimensional length
    v0 = √(2*centrifugal_potential(system.props.μ, [x0, 0., 0.]) - J)
    State(system, [x0, 0., 0., 0., v0, 0.], (0, tmax / system.props.T))
end

# Map the compute_initial_state function to all the parameters in each family above
all_states = Dict(family => begin
    [compute_initial_state(SunMars; p...) for p in params]
end for (family, params) in all_params)


# Plot initial (non-periodic solutions)
trajectory = solve(all_states["g2"][1])
p1 = plot(trajectory)

# Use differential corrector to get periodic solutions
all_half_orbits = Dict(family => begin
    [solve(s, DiffCorrectAxisymmetric()) for s in states]
end for (family, states) in all_states)

# Superimpose corrected plots in an overlay
plot!(p1, all_half_orbits["g2"][1]; color=:blue)

# Lastly, plot the half orbits for a full orbit:
plots = Dict(family => begin
    # half_orbits are given in the (0., t[end]) timespan.
    # We convert these to full orbits by propagating them in the (0., t[end]*2) timespan.
    orbits = [solve(remake(State(o), tspan=(0., o.t[end] * 2))) for o in half_orbits]

    # Plot each of the orbits individually.
    plots = [plot(o; color=:blue) for o in orbits]

    # Plot the whole family as a single row.
    plot(plots...; layout=(1, length(plots)))
end for (family, half_orbits) in all_half_orbits)

# Lay out plots prettily
plot(plots["g2"], plots["g1"], plots["DRO"]; layout=(3, 1),
     legend=false, grid=false, minorgrid=false, showaxis=false, size=(500,500))