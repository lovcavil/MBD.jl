# solver_utils.jl
using OrdinaryDiffEq
using Printf
using Plots
using Statistics
using DelimitedFiles
include("./door240_recons_run.jl")
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
end

# Combine groups of configurations (unchanged)
function combine_groups(groups::Vector{<:Dict})
    configs = Dict()

    # Unpack the groups (solvers, tolerances, timespans, maxsteps)
    solver_group, autodiff_group, linsolve_group, tol_group, timespan_group, initial_state_group = groups

    # Use Iterators.product to get the Cartesian product of all groups
    for product_tuple in Iterators.product(solver_group, autodiff_group, linsolve_group, tol_group, timespan_group, initial_state_group)
        # Unpack the tuple values
        ((solver_name, solver), 
         (autodiff_name, autodiff), 
         (linsolve_name, linsolve), 
         (tol_name, (abstol, reltol)), 
         (timespan_name, timespan), 
         (initial_state_name, initial_state)) = product_tuple

        # Combine the group names for a unique configuration name
        combined_name = "$(solver_name)_$(autodiff_name)_$(linsolve_name)_$(tol_name)_$(timespan_name)_$(initial_state_name)"

        # Create and store the SolverConfig with all the combined settings
        configs[combined_name] = SolverConfig(combined_name, solver, autodiff, linsolve, abstol, reltol, initial_state, timespan)
    end

    return configs
end

# Function to run the solver multiple times and return average time and all times (adjusted for step-based integration)
function run_solver_multiple_times(config, n)
    times = Float64[]

    for i in 1:n
        sol, elapsed_time = run_solver(config)
        push!(times, elapsed_time)
    end

    avg_time = mean(times)
    return times, avg_time
end

# Function to run the solver step-by-step
function run_solver(config)
    solver_name = config.name
    autodiff = config.autodiff
    linsolve = config.linsolve
    solver_constructor = config.solver
    solver = solver_constructor(autodiff=autodiff[1], diff_type=autodiff[2], linsolve=linsolve)
    abstol = config.abstol
    reltol = config.reltol
    initial_state = config.initial_state
    tspan = config.tspan

    println("Running solver: $solver_name with abstol=$abstol, reltol=$reltol")

    # Run solver step-by-step and time it
    elapsed_time = @elapsed begin
        #integrator = solve_HH_system_step(HH_first_order!, initial_state, tspan, solver; abstol=abstol, reltol=reltol)
        integrator = run(tspan, solver; abstol=abstol, reltol=reltol)
    end

    return integrator, elapsed_time
end

# Main function to run all solvers, collect data, and handle multiple runs (unchanged)
function run_solvers(solver_configs,csv_filename, n=10)
    sols = Dict{String, Any}()  # Changed from ODEIntegrator to Any
    energies = Dict{String, Vector{Float64}}()
    avg_times = Dict{String, Float64}()
    all_times = Dict{String, Vector{Float64}}()
    ΔEs = Dict{String, Float64}()

    solver_names = String[]
    ΔE_values = Float64[]
    avg_time_values = Float64[]

    # Collect all configuration options for CSV export
    csv_data = [("Solver", "autodiff", "linsolve", "abstol", "reltol", "initial_state", "tspan", "ΔE", "Avg Time (s)", "All Times")]

    println("Running solvers and collecting results (each run repeated $n times)...\n")

    for config in solver_configs
        solver_name = config.name
        push!(solver_names, solver_name)
        # Run solver once to capture solution and energy (no need to run multiple times for the same solution)
        integrator, _ = run_solver(config)
        sols[solver_name] = integrator.sol
        # Run solver multiple times and get average and all times
        run_times, avg_time = run_solver_multiple_times(config, n)
        all_times[solver_name] = run_times
        avg_times[solver_name] = avg_time
        # Calculate energy and energy difference
        energy, ΔE = calculate_energy_and_deltaE(integrator)
        energies[solver_name] = energy
        ΔEs[solver_name] = ΔE

        push!(ΔE_values, ΔE)
        push!(avg_time_values, avg_time)

        # Add config and result data to CSV, including all run times
        push!(csv_data, (
            solver_name, 
            string(config.autodiff),
            string(config.linsolve),
            string(config.abstol), 
            string(config.reltol), 
            string(config.initial_state), 
            string(config.tspan), 
            string(ΔE), 
            string(avg_time),
            string(run_times) # Storing all run times in the CSV
        ))
    end

    # Print comparison table
    print_solver_comparison(solver_names, ΔE_values, avg_time_values)

    # Save results to CSV
    export_results_to_csv(csv_data,csv_filename)

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


# Function to export results to CSV
function export_results_to_csv(csv_data, csv_filename="solver_results_with_configs.csv")
    writedlm(csv_filename, csv_data, ',')
    println("Results exported to $csv_filename")
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