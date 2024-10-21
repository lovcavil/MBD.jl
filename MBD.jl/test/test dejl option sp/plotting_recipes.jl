using Plots

# Function to plot the comparison of solvers into three separate figures
function plot_HH_comparison(sols, energy_func, labels)
    num_solvers = length(labels)
    colors = get_colors(num_solvers)

    # 1. Orbit Comparison Plot
    plt1 = plot(legend=:outerbottom, size=(900, 600), xlabel="x", ylabel="y", title="Orbit Comparison")
    for i in 1:num_solvers
        x = [u[1] for u in sols[i].u]
        y = [u[2] for u in sols[i].u]
        plot!(plt1, x, y, label=labels[i], color=colors[i])
    end

    # 2. Phase Space Comparison Plot
    plt2 = plot(legend=:outerbottom, size=(900, 600), xlabel="Position", ylabel="Velocity", title="Phase Space Comparison")
    for i in 1:num_solvers
        # x vs dx
        x = [u[1] for u in sols[i].u]
        dx = [u[3] for u in sols[i].u]
        plot!(plt2, x, dx, label="x vs dx - $(labels[i])", color=colors[i])

        # y vs dy (use dashed lines to differentiate)
        y = [u[2] for u in sols[i].u]
        dy = [u[4] for u in sols[i].u]
        plot!(plt2, y, dy, label="y vs dy - $(labels[i])", color=colors[i], linestyle=:dash)
    end

    # 3. Energy Change Comparison Plot
    plt3 = plot(legend=:outerbottom, size=(900, 600), xlabel="Time", ylabel="ΔEnergy", title="Energy Change Comparison")
    for i in 1:num_solvers
        t = sols[i].t
        energy = [energy_func(u[1], u[2], u[3], u[4]) for u in sols[i].u]
        delta_energy = energy .- energy[1]  # ΔEnergy over time
        plot!(plt3, t, delta_energy, label=labels[i], color=colors[i])
    end

    return plt1, plt2, plt3
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

function plot_HH(sol, energy_func)
    # Extract variables from solution
    x = [u[1] for u in sol.u]
    y = [u[2] for u in sol.u]
    dx = [u[3] for u in sol.u]
    dy = [u[4] for u in sol.u]
    t = sol.t
    energy = [energy_func(u[1], u[2], u[3], u[4]) for u in sol.u]
    delta_energy = energy .- energy[1]

    # Create a layout with 3 subplots
    layout = @layout [a b; c]

    # Initialize the plot
    plt = plot(layout=layout, legend=:outerbottom)

    # Orbit plot (subplot 1)
    plot!(plt, x, y, xlabel="x", ylabel="y", title="Orbit of the Hénon-Heiles system", subplot=1)

    # Phase space plot (subplot 2)
    plot!(plt, x, dx, xlabel="x", ylabel="dx", title="Phase Space x vs dx", subplot=2)
    plot!(plt, y, dy, xlabel="y", ylabel="dy", title="Phase Space y vs dy", subplot=2)

    # Energy change plot (subplot 3)
    plot!(plt, t, delta_energy, xlabel="Time", ylabel="ΔEnergy", title="Change in Energy over Time", subplot=3)

    return plt
end
