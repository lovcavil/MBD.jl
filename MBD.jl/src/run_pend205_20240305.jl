# Import necessary packages
include("./problem/AppData_II7.jl")
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
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(params.app)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Correction step
    t = params.tspan[1]
    qdd, Lam, ECond = mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀, du₀ = vcat(q, Lam, qd), vcat(qd, zeros(nc), qdd)

    # Problem setup
    p = [SMDT, STSDAT, SJDT, par]
    prob = ODEProblem(odequation, u₀, params.tspan, p)

    # Solve ODE
    sol = solve(prob; params.solve_kwargs...)
    push!(results.solutions, sol)
end


# using Cassette

# Cassette.@context PrintCtx
# const unique_function_names = Set{Symbol}()

# # Function to log unique function names
# function log_unique_function_name(f, args)
#     func_name = Symbol(f)
#     if !in(func_name, unique_function_names)
#         push!(unique_function_names, func_name)
#         println("Logging unique function name: ", func_name)
#     end
# end

# # Hook to intercept function calls
# Cassette.prehook(::PrintCtx, f, args...) = log_unique_function_name(f, args)
# Cassette.overdub(PrintCtx(), ()->main())

# # Display the unique function names logged
# println("Logged unique function names: ", unique_function_names)

