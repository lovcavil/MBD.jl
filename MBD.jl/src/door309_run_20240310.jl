# Import necessary packages
include("./problem/AD_contact_door309.jl")
include("./problem/AD_contact_door310.jl")
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

    # save callback
    sec1=1:nb*7
    sec2=nb*7+1:nb*7+nc
    sec3=nb*7+nc+1:nb*14+nc
    index1=(1-1)*7+3
    index2=(2-1)*7+3
    index3=(3-1)*7+3
    index4=(4-1)*7+3
    index5=(5-1)*7+3
    indexmf=(4-1)*7+2
    indexmr=(5-1)*7+2
    indexlf=(4-1)*7+2
    indexlr=(5-1)*7+2
    indexu=(4-1)*7+2
    save_func(u, t, integrator) = (u[1] + u[2], t^2,calculate_F_prepare(p_contact[2][1],u[sec1],u[sec3]),
                                calculate_F_prepare(p_contact[2][2],u[sec1],u[sec3]),
                                u[index4]-p_contact[2][1]["pos"],u[index5]-p_contact[2][2]["pos"],
                                0.,0.,0.,0.,0.
                                )
    saved_values = SavedValues(Float64, Tuple{Float64, Float64,
     Float64, Float64, 
     Float64, Float64, 
     Float64, Float64, Float64, Float64, Float64})

    # save_func(u, t, integrator) = (u[1] + u[2], t^2)
    # saved_values = SavedValues(Float64, Tuple{Float64, Float64})
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