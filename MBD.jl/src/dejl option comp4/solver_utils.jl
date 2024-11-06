# solver_utils.jl
using OrdinaryDiffEq
using Printf
using Plots
using DelimitedFiles
using Base.Threads
using FilePathsBase  # for creating directories
using LinearSolve
using Statistics  # For mean calculations on remote workers
include("./door240_recons_run.jl")
using DataFrames,CSV
export E
export solve_HH_system, SolverConfig, run_solvers,combine_groups,run_solver,convert_to_json_compatible
export analyze_HH_system, analyze_HH_system_compare,convert_to_solver_config,generate_param_combinations,convert_from_json_compatible


# Function to create and solve the ODE problem step-by-step using the integrator interface
function solve_HH_system_step(f, initial_state, tspan, solver; abstol=1e-6, reltol=1e-6)
    prob = ODEProblem(f, initial_state, tspan)
    integrator = init(prob, solver; abstol=abstol, reltol=reltol)

    # Step through the integrator manually until completion
    while integrator.t < tspan[end]
        step!(integrator)  # Perform one step of the integration
    end

    return integrator
end

# Define the SolverConfig structure (unchanged)
struct SolverConfig
    name::String               # Name of the configuration
    solver                     # The solver (e.g., Rodas5)
    autodiff
    linsolve
    abstol::Float64            # Absolute tolerance
    reltol::Float64            # Relative tolerance
    initial_state
    tspan
    timeout
end
function convert_to_json_compatible(value)
    if value isa Val
        return string(value)  # Convert Val{:forward} to a string
    elseif value isa Tuple
        return [convert_to_json_compatible(v) for v in value]  # Convert tuples to arrays
    elseif value isa Dict
        return Dict{String, Any}(k => convert_to_json_compatible(v) for (k, v) in value)  # Recursively handle dictionaries
    else
        return value  # Return the value as is if it's already compatible
    end
end
function convert_from_json_compatible(value)
    if value isa String && occursin("Val{", value)
        return eval(Meta.parse(value))  # Convert string like "Val{:forward}" back to Val{:forward}
    elseif value isa Vector
        return [convert_from_json_compatible(v) for v in value]  # Recursively handle arrays
    elseif value isa Dict
        return Dict{Any, Any}(k => convert_from_json_compatible(v) for (k, v) in value)  # Recursively handle dictionaries
    else
        return value  # Return the value as is
    end
end
function generate_param_combinations(groups::Vector{<:Dict})
    param_combinations = Dict()

    # Unpack the groups (solvers, autodiff, linsolve, tol, timespan, initial_state)
    solver_group, autodiff_group, linsolve_group, tol_group, timespan_group, initial_state_group,timeout_group  = groups

    # Use Iterators.product to get the Cartesian product of all groups
    for product_tuple in Iterators.product(solver_group, autodiff_group, linsolve_group, tol_group, timespan_group, initial_state_group,timeout_group)
        # Unpack the tuple values
        ((solver_name, solver), 
         (autodiff_name, autodiff), 
         (linsolve_name, linsolve), 
         (tol_name, (abstol, reltol)), 
         (timespan_name, timespan), 
         (initial_state_name, initial_state),
         (timeout_name, timeout)
         ) = product_tuple

        # Combine the group names for a unique key in the dictionary
        combined_name = "$(solver_name)_$(autodiff_name)_$(linsolve_name)_$(tol_name)_$(timespan_name)_$(initial_state_name)_$(timeout_name)"

        # Store the tuple of parameters in the dictionary
        param_combinations[combined_name] = (solver, autodiff, linsolve, abstol, reltol, initial_state, timespan,timeout)
    end

    return param_combinations
end
function convert_to_solver_config(param_combinations)
    solver_configs = Dict()

    # Iterate through each combination in the dictionary
    for (combined_name, param_tuple) in param_combinations
        # Unpack the tuple values into individual parameters
        solver, autodiff, linsolve, abstol, reltol, initial_state, timespan,timeout = param_tuple

        # Create the SolverConfig
        solver_configs[combined_name] = SolverConfig(combined_name, solver, autodiff, linsolve, abstol, reltol, 
        initial_state, timespan,timeout)
    end

    return solver_configs
end

function run_solver(config)
    solver_name = config.name
    autodiff = config.autodiff
    linsolve = config.linsolve
    linsolve_mapping = Dict(
        "LUFactorization" => LUFactorization(),
        "QRFactorization" => QRFactorization(),
        "" => nothing,  # Handle the case of no linsolve
        "GenericLUFactorization" => GenericLUFactorization(),
        "SVDFactorization" => SVDFactorization(),
        "SimpleLUFactorization" => SimpleLUFactorization(),
        "KrylovJL_GMRES" => KrylovJL_GMRES(),
        "KrylovJL_BICGSTAB" => KrylovJL_BICGSTAB(),
    )
    function_name = config.solver
    fn_symbol = Symbol(function_name)
    solver = @eval $fn_symbol(; autodiff=$(autodiff[1]), diff_type=$(autodiff[2]), linsolve=$(linsolve_mapping[linsolve]))
    abstol = config.abstol
    reltol = config.reltol
    initial_state = config.initial_state
    tspan = config.tspan
    timeout=config.timeout

    elapsed_time = @elapsed begin
        integrator,termination_reason = run(tspan, solver; abstol=abstol, reltol=reltol, timeout=timeout)
    end
    return integrator,elapsed_time,termination_reason

end




function run_solver_multiple_times(config, n)
    times = Float64[]  # Store the elapsed times of each run

    for i in 1:n
        println("Running $(config.name)_$i...")
        # Run the solver with the time limit
        integrator, elapsed_time,termination_reason = run_solver(config)
        push!(times, elapsed_time)  # Store the elapsed time
    end
    avg_time = isempty(times) ? nothing : mean(times)
    return times, avg_time
end

function run_solvers(solver_configs, folder,csv_filename, n=10)
    sols = Dict{String, Any}()
    energies = Dict{String, Vector{Float64}}()
    avg_times = Dict{String, Union{Float64, Nothing}}()  # Allow nothing for failed runs
    all_times = Dict{String, Vector{Float64}}()
    ΔEs = Dict{String, Float64}()

    solver_names = String[]
    ΔE_values = Float64[]
    avg_time_values = Union{Float64, Nothing}[]

    csv_data = [("Method","Solver", "autodiff", "linsolve", "abstol", "reltol", "initial_state", "tspan", "ΔE", "Avg Time (s)","Termination_reason", "All Times")]

    println("Running solvers and collecting results (each run repeated $n times)...\n")
    df=0
    termination_reason="F"
    for i in 1:length(solver_configs)
        config = solver_configs[i]
        solver_name = config.name
        # Run solver once to capture solution and energy with time limit
        integrator, elapsed_time,termination_reason = run_solver(config)
        if integrator == :error
            println("Error: Solver $solver_name failed or exceeded the time limit.")
            continue  # Skip further processing for this solver if there was an error
        end
        df = DataFrame(time = integrator.sol.t)

        # Add each component of the solution `u` as a separate column
        for (i, component) in enumerate(integrator.sol.u[1]) # Assuming sol.u has vectors
            df[!, "u$i"] = [u[i] for u in integrator.sol.u]   # Extract each component of `u` into its own column
        end
        sols[solver_name] = integrator.sol

        # Run solver multiple times with a time limit and get average and all times
        result = run_solver_multiple_times(config, n)
        run_times, avg_time = result

        # Calculate energy and energy difference
        energy, ΔE = calculate_energy_and_deltaE(integrator)

        all_times[solver_name] = run_times
        avg_times[solver_name] = avg_time
        energies[solver_name] = energy
        ΔEs[solver_name] = ΔE

        push!(solver_names, solver_name)
        push!(ΔE_values, ΔE)
        push!(avg_time_values, avg_time)

        # Add config and result data to CSV, including all run times
        push!(csv_data, (
            solver_name, 
            string(config.solver),
            string(config.autodiff),
            string(config.linsolve),
            string(config.abstol), 
            string(config.reltol), 
            string(config.initial_state), 
            string(config.tspan), 
            string(ΔE), 
            string(avg_time),
            termination_reason,
            string(run_times)
        ))

        println("Solver $solver_name completed. ΔE = $ΔE, Avg Time = $avg_time")
    end

    # Save results to CSV
    export_results_to_csv(csv_data,folder, csv_filename)
    export_dfresults_to_csv(df, folder,csv_filename)
    export_string_to_file(termination_reason, folder,csv_filename)
    return sols, energies, avg_times, ΔEs, all_times
end




# Function to calculate energy and ΔE (energy difference)
function calculate_energy_and_deltaE(integrator)
    # Assuming each u is a vector or tuple of 4 elements (x, y, px, py)
    energy = [E(u) for u in integrator.sol.u]  # Explicitly pass all 4 components of u to E
    ΔE = energy[1] - energy[end]
    return energy, ΔE
end


# Plot and comparison functions remain unchanged
# Adjusted print function to reflect average time
function print_solver_comparison(solver_names, ΔE_values, avg_time_values)
    println("Solver Comparison:")
    println(rpad("Solver", 45), lpad("ΔE", 25), lpad("Avg Computation Time (s)", 25))
    println("-"^95)
    for i in 1:length(solver_names)
        ΔE_str = @sprintf("%.3e", ΔE_values[i])
        avg_time_str = @sprintf("%.6f", avg_time_values[i])
        println(rpad(solver_names[i], 45), lpad(ΔE_str, 25), lpad(avg_time_str, 25))
    end
    println("\n")
end




function export_results_to_csv(csv_data, folder, csv_filename="solver_results_with_configs.csv")
    # Create the folder path based on the script directory
    folder_path = joinpath(@__DIR__, folder)
    
    # Check if the folder exists, and create it if it doesn't
    if !isdir(folder_path)
        mkpath(folder_path)  # Create the folder and any missing intermediate directories
    end
    
    # Construct the full file path
    file_path = joinpath(folder_path, csv_filename)
    
    # Write the CSV file
    writedlm(file_path, csv_data, ',')
end
function export_dfresults_to_csv(df, folder,csv_filename)
    # Create the folder path based on the script directory
    folder_path = joinpath(@__DIR__, folder)
    
    # Check if the folder exists, and create it if it doesn't
    if !isdir(folder_path)
        mkpath(folder_path)  # Create the folder and any missing intermediate directories
    end
    file_path = joinpath(folder_path, "$(csv_filename).2")
    CSV.write(file_path, df)
end
function export_string_to_file(content::String, folder::String, filename::String)
    # Create the folder path based on the script directory
    folder_path = joinpath(@__DIR__, folder)
    
    # Check if the folder exists, and create it if it doesn't
    if !isdir(folder_path)
        mkpath(folder_path)  # Create the folder and any missing intermediate directories
    end
    
    # Create the full file path
    file_path = joinpath(folder_path, "$(filename).state")
    
    # Open the file and write the string content
    open(file_path, "w") do file
        write(file, content)
    end
end
# Function to analyze and compare results from multiple solvers
function analyze_HH_system_compare(sols, energies, ΔEs, times)
    solver_names = collect(keys(sols))

    # Prepare data for plotting
    sol_list = [sols[name] for name in solver_names]
    energy_list = [energies[name] for name in solver_names]
    labels = solver_names

    # Plot comparative results using the comparison recipe
    plt3 = plot_HH_comparison(sol_list, energy_list,labels)

    display(plt3)
end
function plot_HH_comparison(sol_list,energy_list, labels)
    num_solvers = length(labels)
    colors = get_colors(num_solvers)


    # 3. Energy Change Comparison Plot
    plt3 = plot(legend=:outerbottom, size=(900, 600), xlabel="Time", ylabel="ΔEnergy", title="Energy Change Comparison")
    for i in 1:num_solvers
        t = sol_list[i].t
        energy = energy_list[i]
        delta_energy = energy .- energy[1]  # ΔEnergy over time
        plot!(plt3, t, delta_energy, label=labels[i], color=colors[i])
    end

    return plt3
end

# Helper function to get colors based on the number of solvers
function get_colors(num_solvers::Int)
    base_colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :yellow, :black]

    # Repeat colors if needed
    if num_solvers > length(base_colors)
        return vcat(base_colors, repeat(base_colors, ceil(Int, num_solvers / length(base_colors))))
    else
        return base_colors[1:num_solvers]
    end
end

# Define the energy function E(x, y, dx, dy)
function E(u)
    # Constants and initial conditions
    h0, hvar, g = 0.001, 0, 9.81
    contact_json = "contact_g65"
    app = 240
    # Load application data
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd, p_contact = AD240(app, contact_json)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]
    sec1=1:nb*7
    sec2=nb*7+1:nb*7+nc
    sec3=nb*7+nc+1:nb*14+nc
    qd=u[sec3]
    q=u[sec1]
    l=u[sec2]
    # Calculate Total Energy
    M = mathfunction.MEval(q, SMDT, par)
    KE= 0.5 * qd' * M * qd

    i = 1
    PE = 0
    while i <= nb
        PE += SMDT[1, i] * g * q[7*(i-1)+3]
        i += 1
    end

    SE= 0
    T = 1
    while T <= NTSDA
        i, j, s1pr, s2pr, K, C, el0, F = STSDATPart(STSDAT, T)
        r1, p1 = qPart(q, i)
        r2 = [0, 0, 0]
        p2 = [1, 0, 0, 0]
        r1d = [0, 0, 0]
        p1d = zeros(4)
        if j >= 1
            r2, p2 = qPart(q, j)
        end
        A1 = ATran(p1)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        el = sqrt(d12' * d12)
        SET = 0.5 * K * (el - el0)^2
        SE += SET
        T += 1
    end
    TE=0.0
    TE = KE + PE + SE

end