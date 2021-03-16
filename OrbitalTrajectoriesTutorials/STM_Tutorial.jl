using OrbitalTrajectories
using DifferentialEquations  # To propagate trajectories
using Plots                  # To plot trajectories
using LinearAlgebra          # To compute matrix eigenvalues

# Use cases (see YouTube video & Pellegrini paper if you can get your hands on it)
test_cases = [(
    system = (:Jupiter, :Europa, :Sun),
    u0     = [-0.1017472008677258, 0., 0., 0., -0.01806028472285857, 0.],
    μ      = 2.528009215182033e-5,
    t_p    = 25.13898226959327,
    λ_max  = 2.468621114047195,
    plot_args = (label="7:4 resonant\n(Jupiter-Europa)", xlim=(-1.15, 1.07))
), (
    system = (:Jupiter, :Europa, :Sun),
    u0     = [0.04867586089512202, 0., 0., 0., -0.09354853663949217, 0.],
    μ      = 2.528009215182033e-5,
    t_p    = 70.53945041512506,
    λ_max  = 2.804814246519340e7,
    plot_args = (label="8:11 + flybys\n(Jupiter-Europa)", xlim=(-1.37, 1.35))
), (
    system = (:Earth, :Moon, :Sun),
    u0     = [-0.013059020050628, 0., 0.07129515195874, 0., -0.526306975588415, 0.],
    μ      = 0.01215509906405700,
    t_p    = 2.517727406553485,
    λ_max  = 17.632688124231755,
    plot_args = (label="L1 halo orbit\n(Earth-Moon)", xlim=(0.8, 1.015))
)]

case = test_cases[1]
# Build the system
# NOTE: overriding the default mass parameter μ with the one provided by [Pellegrini2016]
system = CR3BP(case.system[1:2]...; μ=case.μ)

# Propagation timespan
t_init = π/4  # Arbitrary initial time, does not matter for CR3BP
tspan  = (t_init, t_init + case.t_p)

# Initial state vector in [Pellegrini2016] is relative to secondary body
u0 = copy(case.u0)
u0[1] += 1 - case.μ

# Build the initial state and propagate it
state = State(system, SynodicFrame(), u0, tspan)
trajectory = solve(state)
plot(trajectory; legend=:bottomleft, case.plot_args...)

# Use automatic differentiation to calculate the STM at the final point of the trajectory (w.r.t initial state)
STM = sensitivity(AD, state)

# Check if STM eigenvalues matches what [Pellegrini2016] gives
λ_max = maximum(norm.(eigvals(Matrix(STM))))
@assert isapprox(λ_max, case.λ_max; rtol=1e-5)

# Compute the STM trace (i.e. the STM at every point along the trajectory)
STM_trace = sensitivity_trace(AD, state)

# Plot the STM elements
plot(STM_trace; trace=true)