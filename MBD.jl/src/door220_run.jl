# Import necessary packages
include("./problem/AD_contact_door_base.jl")
include("./solver/solver.jl")
include("./eval/contact.jl")
using LinearAlgebra, DifferentialEquations, OrdinaryDiffEq, Sundials, Plots, CSV, DataFrames
using DiffEqCallbacks

# Define a struct to encapsulate parameters for the ODE problem
struct ODEParams
    app::Int
    tspan::Tuple{Float64, Float64}
    solve_kwargs::Dict
end

struct ODERunResults
    solutions::Vector{Any}
    l_saved_data
end

include("./post/save.jl")
include("./post/draw.jl")

# Function to generate a tuple type for SavedValues
function create_saved_values_type(num_float64::Int)
    # Generate a tuple type with `num_float64` Float64 elements
    return Tuple{fill(Float64, num_float64)...}
end



# Define a function to run a specific ODE problem
function run(params::ODEParams, results::ODERunResults)
    # Constants and initial conditions
    h0, hvar, g = 0.001, 0, 9.81

    # Load application data
    #nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd, p_contact = AD(params.app)
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd, p_contact = appdata(params.app)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Correction step
    t = params.tspan[1]
    qdd, Lam, ECond = mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀, du₀ = vcat(q, Lam, qd), vcat(qd, zeros(nc), qdd)

    # Problem setup
    p = [SMDT, STSDAT, SJDT, par, p_contact]
    prob = ODEProblem(odequation, u₀, params.tspan, p)


    save_func(u, t, integrator) = (0.,0.)

    # Example: Create a tuple type for saved values with 11 Float64 elements
    saved_values = SavedValues(Float64, create_saved_values_type(2))


    cb = SavingCallback(save_func, saved_values)
    # Solve ODE
    sol = solve(prob; params.solve_kwargs..., callback = cb)


    # # Access saved values
    saved_data = saved_values.saveval
    merge_df=merge_sol_result(sol,saved_values)

    push!(results.solutions, sol)
    push!(results.l_saved_data, merge_df)
end

function merge_sol_result(sol,saved_values)
    df = DataFrame(Time=sol.t)

    # Add solution data to the DataFrame
    for i in 1:length(sol.u[1])
        df[!, Symbol("sol_$i")] = getindex.(sol.u, i)
    end

    # Assuming saved_values.saveval contains the saved data,
    # Append the saved data directly to the DataFrame without type-specific initialization
    for i in 1:length(saved_values.saveval[1]) # Assuming the first tuple represents the saved data structure
        # Directly read and append the saved data to the DataFrame
        df[!, Symbol("ex_$i")] = [sd[i] for sd in saved_values.saveval]
    end
    return df
end