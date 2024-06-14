# Import necessary packages

#cover the extra data output

include("./problem/AD_contact_door_base.jl")
include("./solver/solver.jl")
include("./eval/contact.jl")
include("./mathfunction_II7.jl")

using LinearAlgebra, DifferentialEquations, OrdinaryDiffEq, Sundials, Plots, CSV, DataFrames
using DiffEqCallbacks

# Define a struct to encapsulate parameters for the ODE problem
struct ODEParams
    app::Int
    contact_json::String
    tspan::Tuple{Float64, Float64}
    solve_kwargs::Dict
end

struct ODERunResults
    solutions::Vector{Any}
    l_saved_data
end

include("./post/save.jl")
include("./post/draw.jl")
include("./post/recordTime.jl")

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
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd, p_contact = appdata(params.app, params.contact_json)
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
    index_x_mf=(6-1)*7+1
    index_y_mf=(6-1)*7+2
    index_x_mr=(7-1)*7+1
    index_y_mr=(7-1)*7+2
    indexlf=(4-1)*7+2
    indexlr=(5-1)*7+2
    indexu=(4-1)*7+2
    ld_contact=p_contact[2]
    contact_mg,contact_lg,contact_mf,contact_mr=ld_contact
    save_func(u, t, integrator) = (
            u[index_y_mf]-contact_mf["guide"](u[index_x_mf]),#diff_y from center to guide
            u[index_y_mr]-contact_mr["guide"](u[index_x_mr]),#diff_y from center to guide,#
            calculate_contact_geo(contact_mf,u[sec1],u[sec3])[1],#fx #q,qd
            calculate_contact_geo(contact_mf,u[sec1],u[sec3])[2],
            calculate_contact_geo(contact_mr,u[sec1],u[sec3])[1],#fy #q,qd
            calculate_contact_geo(contact_mr,u[sec1],u[sec3])[2],
                                u[nb*7+nc+index_y_mf],
                                u[nb*7+nc+index_y_mr],
                                float(size(mathfunction.PhiqEval(t,u[sec1],SJDT,par), 1)),
                                float(size(mathfunction.PhiqEval(t,u[sec1],SJDT,par), 2)),
                                float(size(u[sec2], 1)),
                                float(size((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2], 1)),
                                float(size((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2], 2)),
                                ((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2])[53],
                                ((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2])[54],
                                ((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2])[55],
                                ((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2])[56],
                                ((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2])[57],
                                ((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2])[58],
                                ((mathfunction.PhiqEval(t,u[sec1],SJDT,par))'*u[sec2])[59]
                                )

    # Example: Create a tuple type for saved values with 11 Float64 elements
    saved_values = SavedValues(Float64, create_saved_values_type(13+7))

    cb1 = SavingCallback(save_func, saved_values)

    condition(u,t,integrator) = (u[1] < 2.14) 
    affect!(integrator) = terminate!(integrator)
    cb2 = DiscreteCallback(condition,affect!)

    # Solve ODE
    sol = solve(prob; params.solve_kwargs..., callback = CallbackSet(cb1,cb2))


    # # Access saved values
    saved_data = saved_values.saveval
    merge_df=merge_sol_result(sol,saved_values)

    push!(results.solutions, sol)
    push!(results.l_saved_data, merge_df)
end


function merge_sol_result(sol, saved_values)
    df = DataFrame(Time=sol.t)

    # Add solution data to the DataFrame
    for i in 1:length(sol.u[1])
        df[!, Symbol("sol_$i")] = getindex.(sol.u, i)
    end

    # Determine the length of data and handle mismatch in array sizes
    n_time_steps = length(df.Time)
    n_saved_steps = length(saved_values.saveval)
    
    # Append the saved data directly to the DataFrame, adjusting for length mismatch
    for i in 1:length(saved_values.saveval[1])  # Assuming the first tuple represents the saved data structure
        if n_saved_steps == n_time_steps
            # If sizes match, append directly
            df[!, Symbol("ex_$i")] = [sd[i] for sd in saved_values.saveval]
        elseif n_saved_steps > n_time_steps
            # If saved data is longer, truncate to match df
            df[!, Symbol("ex_$i")] = [sd[i] for sd in saved_values.saveval[1:n_time_steps]]
        else
            # If saved data is shorter, pad with zeros
            extended_data = [sd[i] for sd in saved_values.saveval]
            append!(extended_data, zeros(n_time_steps - n_saved_steps))
            df[!, Symbol("ex_$i")] = extended_data
        end
    end

    return df
end
