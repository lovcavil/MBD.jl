# Import necessary packages

#cover the extra data output

include("../problem/AD_240_renew.jl")
include("../solver/solver.jl")
include("../eval/contact.jl")
include("../mathfunction_II7.jl")

using LinearAlgebra, DifferentialEquations, OrdinaryDiffEq, Sundials, Plots, CSV, DataFrames
using DiffEqCallbacks

include("../post/save.jl")

function create_saved_values_type(num_float64::Int)
    # Generate a tuple type with `num_float64` Float64 elements
    return Tuple{fill(Float64, num_float64)...}
end


function run(tspan, solver;app = 240, abstol=1e-6, reltol=1e-6, timeout=1.0)
    # Constants and initial conditions
    h0, hvar, g = 0.001, 0, 9.81
    contact_json = "contact_simpley"
    # Load application data
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd, p_contact = AD240(app, contact_json)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]
    # Correction step
    t = tspan[1]
    qdd, Lam, ECond = mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀, du₀ = vcat(q, Lam, qd), vcat(qd, zeros(nc), qdd)

    # Problem setup
    p = [SMDT, STSDAT, SJDT, par, p_contact]
    prob = ODEProblem(odequation!, u₀, tspan, p)

    # Save callback
    save_func(u, t, integrator) = (
        # Collect relevant values here
        0, 0
    )
    saved_values = SavedValues(Float64, create_saved_values_type(2))
    cb1 = SavingCallback(save_func, saved_values)
    
    # Initialize the integrator with callbacks
    integrator = init(prob, solver; callback=CallbackSet(cb1, ), abstol=abstol, reltol=reltol, dtmin = 1e-10,dtmax = 1e-3)

    start_time = time()
    timeout_threshold = timeout  # Timeout after 5 seconds
    termination_reason = "completed"  # Default to normal completion

    # Constraint limits
    PosConstrMax = 10^-4  # Limit on position constraint error
    VelConstrMax = 10^-4  # Limit on velocity constraint error
    sec1 = 1:nb*7
    sec2 = nb*7+1:nb*7+nc
    sec3 = nb*7+nc+1:nb*14+nc

    # Manually step through the ODE solution
    while integrator.t < tspan[2]
        step!(integrator)

        # Update state variables
        u = integrator.u
        qd = u[sec3]
        q = u[sec1]
        l = u[sec2]
        tn = integrator.t
        Phiq = mathfunction.PhiqEval(tn, q, SJDT, par)
        if norm(Phiq * qd) > VelConstrMax
            mu = pinv(Phiq * Phiq') * (Phiq * qd)
            delqd = -Phiq' * mu
            qd += delqd
        end
        u[sec3] = qd
        u[sec1] = q
        u[sec2] = l
        integrator.u = u
        if (time() - start_time) > timeout_threshold
            println("Terminating due to timeout at t = $(integrator.t)")
            termination_reason="timeout"
            break
        end
    end

    sol = integrator.sol

    # Access saved values
    # saved_data = saved_values.saveval
    # merge_df = merge_sol_result(sol, saved_values)
    
    return integrator, termination_reason
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
# tspan=(0.0, 1.0)
# solver=Tsit5()
# run(tspan,solver)