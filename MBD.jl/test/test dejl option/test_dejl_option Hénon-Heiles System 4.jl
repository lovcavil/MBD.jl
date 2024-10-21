using OrdinaryDiffEq, Plots, RecipesBase, Printf

# Define the potential and energy functions
function V(x, y)
    return 1//2 * (x^2 + y^2 + 2x^2 * y - 2//3 * y^3)
end

function E(x, y, dx, dy)
    return V(x, y) + 1//2 * (dx^2 + dy^2)
end

# Convert the second-order ODE to a first-order ODE
function HH_first_order!(du, u, p, t)
    x, y, dx, dy = u
    du[1] = dx
    du[2] = dy
    du[3] = -x - 2x * y
    du[4] = -y - x^2 + y^2
end

# Define a flexible function to create and solve the problem
function solve_HH_system(f, initial_state, tspan, solver)
    prob = ODEProblem(f, initial_state, tspan)
    sol = solve(prob, solver)
    return sol
end

# Define the SolverConfig struct with @kwdef to accept keyword arguments
Base.@kwdef struct SolverConfig
    name::String
    solver
end

# Function to run solvers and store solutions
function run_solvers(initial_state, tspan, solver_configs)
    sols = Dict{String, ODESolution}()
    energies = Dict{String, Vector{Float64}}()
    times = Dict{String, Float64}()
    ΔEs = Dict{String, Float64}()

    solver_names = String[]
    ΔE_values = Float64[]
    time_values = Float64[]

    println("Running solvers and collecting results...\n")

    for config in solver_configs
        solver_name = config.name
        push!(solver_names, solver_name)
        solver = config.solver

        println("Running solver: $solver_name")
        # Run solver and time it
        elapsed_time = @elapsed begin
            sol = solve_HH_system(HH_first_order!, initial_state, tspan, solver)
            sols[solver_name] = sol
        end
        energy = [E(u[1], u[2], u[3], u[4]) for u in sols[solver_name].u]
        energies[solver_name] = energy

        # Compute energy difference
        ΔE = energy[1] - energy[end]
        ΔEs[solver_name] = ΔE
        times[solver_name] = elapsed_time

        push!(ΔE_values, ΔE)
        push!(time_values, elapsed_time)
    end

    # Print comparison table
    println("Solver Comparison:")
    println(rpad("Solver", 45), lpad("ΔE", 25), lpad("Computation Time (s)", 25))
    println("-"^95)
    for i in 1:length(solver_names)
        ΔE_str = @sprintf("%.3e", ΔE_values[i])
        time_str = @sprintf("%.6f", time_values[i])
        println(rpad(solver_names[i], 45), lpad(ΔE_str, 25), lpad(time_str, 25))
    end
    println("\n")

    return sols, energies, times, ΔEs
end

# Function to analyze and plot results for each solver
function analyze_HH_system(sols, energies, ΔEs, times)
    for solver_name in keys(sols)
        sol = sols[solver_name]
        energy = energies[solver_name]
        ΔE = ΔEs[solver_name]
        elapsed_time = times[solver_name]

        println("Analysis for solver: $solver_name")
        println("ΔE = $ΔE")
        println("Computation Time: $(elapsed_time) seconds\n")

        # Plot results using the single solver recipe
        plot_result = hhplot(sol, E)
        display(plot_result)
    end
end

# Define a plot recipe for single solver analysis
@userplot HHPlot

@recipe function f(hhplot::HHPlot)
    # Access arguments passed to hhplot
    sol = hhplot.args[1]
    energy_func = hhplot.args[2]

    # Define the layout
    layout := @layout [a b; c]

    # Extract variables
    x = [u[1] for u in sol.u]
    y = [u[2] for u in sol.u]
    dx = [u[3] for u in sol.u]
    dy = [u[4] for u in sol.u]
    t = sol.t
    energy = [energy_func(u[1], u[2], u[3], u[4]) for u in sol.u]
    delta_energy = energy .- energy[1]

    # Orbit plot
    @series begin
        subplot := 1
        xlabel := "x"
        ylabel := "y"
        title := "Orbit of the Hénon-Heiles system"
        x, y
    end

    # Phase space plots
    @series begin
        subplot := 2
        xlabel := "x"
        ylabel := "dx"
        title := "Phase Space x vs dx"
        x, dx
    end
    @series begin
        subplot := 2
        xlabel := "y"
        ylabel := "dy"
        title := "Phase Space y vs dy"
        y, dy
    end

    # Energy change plot
    @series begin
        subplot := 3
        xlabel := "Time"
        ylabel := "ΔEnergy"
        title := "Change in Energy over Time"
        t, delta_energy
    end
end

# Function to analyze and compare results from multiple solvers
function analyze_HH_system_compare(sols, energies, ΔEs, times)
    solver_names = collect(keys(sols))

    # Prepare data for plotting
    sol_list = [sols[name] for name in solver_names]
    labels = solver_names

    # Plot comparative results using the comparison recipe
    plot_result = hhplotcompare(sol_list, E, labels)
    display(plot_result)
end

# Define a plot recipe for solver comparison
@userplot HHPlotCompare

@recipe function f(hhplotcompare::HHPlotCompare)
    # Access arguments passed to hhplotcompare
    sols = hhplotcompare.args[1]          # List of solutions
    energy_func = hhplotcompare.args[2]   # Energy function
    labels = hhplotcompare.args[3]        # Labels for each solver

    # Define the layout
    layout := @layout [a b; c]

    # Colors for plotting
    colors = [:blue, :red, :green, :orange, :purple, :cyan, :magenta, :yellow, :black]

    # Orbit plot
    for i in 1:length(sols)
        @series begin
            subplot := 1
            xlabel := "x"
            ylabel := "y"
            title := "Orbit Comparison"
            color := colors[i]
            label := labels[i]
            x = [u[1] for u in sols[i].u]
            y = [u[2] for u in sols[i].u]
            x, y
        end
    end

    # Phase space plots
    for i in 1:length(sols)
        # x vs dx
        @series begin
            subplot := 2
            xlabel := "Position"
            ylabel := "Velocity"
            title := "Phase Space Comparison"
            color := colors[i]
            label := "x vs dx - $(labels[i])"
            x = [u[1] for u in sols[i].u]
            dx = [u[3] for u in sols[i].u]
            x, dx
        end
        # y vs dy
        @series begin
            subplot := 2
            xlabel := "Position"
            ylabel := "Velocity"
            color := colors[i]
            linestyle := :dash
            label := "y vs dy - $(labels[i])"
            y = [u[2] for u in sols[i].u]
            dy = [u[4] for u in sols[i].u]
            y, dy
        end
    end

    # Energy change plot
    for i in 1:length(sols)
        @series begin
            subplot := 3
            xlabel := "Time"
            ylabel := "ΔEnergy"
            title := "Energy Change Comparison"
            color := colors[i]
            label := labels[i]
            sol = sols[i]
            t = sol.t
            energy = [energy_func(u[1], u[2], u[3], u[4]) for u in sol.u]
            delta_energy = energy .- energy[1]
            t, delta_energy
        end
    end
end

# Main function to run and analyze solvers
function main()
    # Initial conditions and time span
    initial_state = [0.0, 0.1, 0.5, 0.0]  # [x0, y0, dx0, dy0]

    # Define solver configurations with various Jacobian options
    solver_configs = [
        SolverConfig(
            name="Rodas5_autodiff_true",
            solver=Rodas5(autodiff=true)
        ),
        SolverConfig(
            name="Rodas5_autodiff_false_forward",
            solver=Rodas5(autodiff=false, diff_type=Val{:forward})
        ),
        SolverConfig(
            name="Rodas5_autodiff_false_central",
            solver=Rodas5(autodiff=false, diff_type=Val{:central})
        ),
        SolverConfig(
            name="Rodas5_autodiff_false_complex",
            solver=Rodas5(autodiff=false, diff_type=Val{:complex})
        ),
        # Add more configurations with other Jacobian options if desired
    ]

    # Run solvers and store solutions
    # tspan = (0.0, 5.0)
    # sols, energies, times, ΔEs = run_solvers(initial_state, tspan, solver_configs)
    tspan = (0.0, 1000.0)
    sols, energies, times, ΔEs = run_solvers(initial_state, tspan, solver_configs)
    # Analyze and plot results for each solver
    #analyze_HH_system(sols, energies, ΔEs, times)

    # Analyze and compare solvers
    analyze_HH_system_compare(sols, energies, ΔEs, times)
end

# Execute the main function
main()