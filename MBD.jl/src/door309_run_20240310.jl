# Import necessary packages
include("./problem/AD_contact_door309.jl")
include("./problem/AD_contact_door310.jl")
include("./solver/solver.jl")

using LinearAlgebra, DifferentialEquations, OrdinaryDiffEq, Sundials, Plots, CSV, DataFrames

# Define a struct to encapsulate parameters for the ODE problem
struct ODEParams
    app::Int
    tspan::Tuple{Float64, Float64}
    solve_kwargs::Dict
end

struct ODERunResults
    solutions::Vector{Any}
end

include("./post/save.jl")
include("./post/draw.jl")

# Define a function to run a specific ODE problem
function run(params::ODEParams, results::ODERunResults)
    # Constants and initial conditions
    h0, hvar, g = 0.001, 0, 9.81

    # Load application data
    #nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd, p_contact = AD(params.app)
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd, p_contact = AD310(params.app)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Correction step
    t = params.tspan[1]
    qdd, Lam, ECond = mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀, du₀ = vcat(q, Lam, qd), vcat(qd, zeros(nc), qdd)

    # Problem setup
    p = [SMDT, STSDAT, SJDT, par, p_contact]
    prob = ODEProblem(odequation, u₀, params.tspan, p)

    # Solve ODE
    sol = solve(prob; params.solve_kwargs...)
    push!(results.solutions, sol)
end

