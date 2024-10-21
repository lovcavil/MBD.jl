# Import necessary packages

#cover the extra data output

include("./problem/AD_240_renew.jl")
include("./solver/solver.jl")
include("./eval/contact.jl")
include("./mathfunction_II7.jl")

using LinearAlgebra, DifferentialEquations, OrdinaryDiffEq, Sundials, Plots, CSV, DataFrames
using DiffEqCallbacks

# Define a struct to encapsulate parameters for the ODE problem
struct ODEParams
    app::Int
    contact_json::String
    tspan::Tuple{Float64,Float64}
    solve_kwargs::Dict
end

struct ODERunResults
    solutions::Vector{Any}
    l_saved_data
end

include("./post/save.jl")
include("./post/draw.jl")
include("./post/recordTime.jl")
# include("./post/extra_return.jl")

function create_saved_values_type(num_float64::Int)
    # Generate a tuple type with `num_float64` Float64 elements
    return Tuple{fill(Float64, num_float64)...}
end

function adds!(str, names, count)
    push!(names, str)
    count += 1
end

function run(params::ODEParams, results::ODERunResults)
    # Constants and initial conditions
    h0, hvar, g = 0.001, 0, 9.81

    # Load application data
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd, p_contact = AD240(params.app, params.contact_json)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Correction step
    t = params.tspan[1]
    qdd, Lam, ECond = mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀, du₀ = vcat(q, Lam, qd), vcat(qd, zeros(nc), qdd)

    # Problem setup
    p = [SMDT, STSDAT, SJDT, par, p_contact]
    prob = ODEProblem(odequation, u₀, params.tspan, p)

    # Initialize the integrator

    # Save callback
    save_func(u, t, integrator) = (
        # Collect relevant values here
        0, 0
    )

    saved_values = SavedValues(Float64, create_saved_values_type(2))
    cb1 = SavingCallback(save_func, saved_values)

    # Set up a termination condition
    condition(u, t, integrator) = (u[1] < 999992.14)
    affect!(integrator) = terminate!(integrator)
    cb2 = DiscreteCallback(condition, affect!)
    integrator = init(prob; params.solve_kwargs..., callback=CallbackSet(cb1))
    # Manually step through the ODE solution

    # -------------
    PosConstrMax = 10^-4  # Limit on position constraint error
    VelConstrMax = 10^-4  # Limit on velocity constraint error
    sec1=1:nb*7
    sec2=nb*7+1:nb*7+nc
    sec3=nb*7+nc+1:nb*14+nc

    # -------------
    while integrator.t < params.tspan[2]
        step!(integrator)
        # Collect any additional results here if necessary
        u=integrator.u
        qd=u[sec3]
        q=u[sec1]
        l=u[sec2]
        tn=integrator.t
        Phiq = mathfunction.PhiqEval(tn, q, SJDT, par)
        if norm(Phiq * qd) > VelConstrMax
            mu = pinv(Phiq * Phiq') * (Phiq * qd)
            delqd = -Phiq' * mu
            qd += delqd
            #corvel += 1
        end
        u[sec3]=qd
        u[sec1]=q
        u[sec2]=l
        integrator.u=u
        
    end
    sol=integrator.sol
    # # Access saved values
    saved_data = saved_values.saveval
    merge_df = merge_sol_result(sol, saved_values)
    push!(results.solutions, sol)
    push!(results.l_saved_data, merge_df)
end


function merge_sol_result(sol, saved_values) #@note name
    df = DataFrame(Time=sol.t)

    # Add solution data to the DataFrame
    for i in 1:length(sol.u[1])
        df[!, Symbol("sol_$i")] = getindex.(sol.u, i)
    end

    # Determine the length of data and handle mismatch in array sizes
    n_time_steps = length(df.Time)
    n_saved_steps = length(saved_values.saveval)
    names = []
    push!(names,"p1","p2")

    # Append the saved data directly to the DataFrame, adjusting for length mismatch
    for i in 1:length(saved_values.saveval[1])  # Assuming the first tuple represents the saved data structure
        col_name=Symbol("ex_$(i)_$(names[i])")
        if n_saved_steps == n_time_steps
            # If sizes match, append directly
            df[!, col_name] = [sd[i] for sd in saved_values.saveval]
        elseif n_saved_steps > n_time_steps
            # If saved data is longer, truncate to match df
            df[!, col_name] = [sd[i] for sd in saved_values.saveval[1:n_time_steps]]
        else
            # If saved data is shorter, pad with zeros
            extended_data = [sd[i] for sd in saved_values.saveval]
            append!(extended_data, zeros(n_time_steps - n_saved_steps))
            df[!, col_name] = extended_data
        end
    end
    return df
end
