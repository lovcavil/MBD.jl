# solver_utils.jl
using OrdinaryDiffEq
using Printf
using Plots
using Statistics

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
        integrator = solve_HH_system_step(HH_first_order!, initial_state, tspan, solver; abstol=abstol, reltol=reltol)
    end

    return integrator, elapsed_time
end

# Main function to run all solvers, collect data, and handle multiple runs (unchanged)
function run_solvers(solver_configs, n=10)
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

        # Run solver multiple times and get average and all times
        run_times, avg_time = run_solver_multiple_times(config, n)
        all_times[solver_name] = run_times
        avg_times[solver_name] = avg_time

        # Run solver once to capture solution and energy (no need to run multiple times for the same solution)
        integrator, _ = run_solver(config)
        sols[solver_name] = integrator

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
    export_results_to_csv(csv_data)

    return sols, energies, avg_times, ΔEs, all_times
end

# Function to calculate energy and ΔE (energy difference)
function calculate_energy_and_deltaE(integrator)
    # Assuming each u is a vector or tuple of 4 elements (x, y, px, py)
    println(integrator.u)
    energy = [E(u[1], u[2], u[3], u[4]) for u in integrator.u]  # Explicitly pass all 4 components of u to E
    ΔE = energy[1] - energy[end]
    return energy, ΔE
end


# Plot and comparison functions remain unchanged



# Function to export results to CSV
function export_results_to_csv(csv_data, csv_filename="solver_results_with_configs.csv")
    writedlm(csv_filename, csv_data, ',')
    println("Results exported to $csv_filename")
end

function analyze_HH_system(sols, energies, ΔEs, times)
    for solver_name in keys(sols)
        sol = sols[solver_name]
        energy = energies[solver_name]
        ΔE = ΔEs[solver_name]
        elapsed_time = times[solver_name]

        # Call the hhplot recipe to generate the plot
        plt_s = plot_HH(sol, E)
        display(plt_s)
    end
end

# Function to analyze and compare results from multiple solvers
function analyze_HH_system_compare(sols, energies, ΔEs, times)
    solver_names = collect(keys(sols))

    # Prepare data for plotting
    sol_list = [sols[name] for name in solver_names]
    labels = solver_names

    # Plot comparative results using the comparison recipe
    plt1, plt2, plt3 = plot_HH_comparison(sol_list, E,labels)
    display(plt1)
    display(plt2)
    display(plt3)
end
