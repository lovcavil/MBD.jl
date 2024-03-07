# Adjust the draw function to handle multiple solutions and dynamically set figure limits
function draw(results::ODERunResults; groups, plot3d::Bool=true)
    plots = []
    for group in groups
        if plot3d && length(group) == 3
            # Initialize limits for 3D plot
            limits = Dict("x" => [Inf, -Inf], "y" => [Inf, -Inf], "z" => [Inf, -Inf])
            
            # 3D plot creation
            p = plot(layout=(1,1), legend=:topright)

            for sol in results.solutions
                x_vals = getindex.(sol.u, group[1])
                y_vals = getindex.(sol.u, group[2])
                z_vals = getindex.(sol.u, group[3])

                # Update limits
                limits["x"] = [min(minimum(x_vals), limits["x"][1]), max(maximum(x_vals), limits["x"][2])]
                limits["y"] = [min(minimum(y_vals), limits["y"][1]), max(maximum(y_vals), limits["y"][2])]
                limits["z"] = [min(minimum(z_vals), limits["z"][1]), max(maximum(z_vals), limits["z"][2])]

                plot!(p, x_vals, y_vals, z_vals, label="Run $(findfirst(isequal(sol), results.solutions))")
            end

            # Apply limits with padding
            xlims!(p, (limits["x"][1], limits["x"][2]))
            ylims!(p, (limits["y"][1], limits["y"][2]))
            zlims!(p, (limits["z"][1], limits["z"][2]))

            push!(plots, p)
        else
            # Initialize limits for 2D plot
            x_limits = [Inf, -Inf]
            y_limits = [Inf, -Inf]

            p = plot(legend=:topright)

            for sol in results.solutions
                x_vals = sol.t # Assuming sol.t holds the x-axis values (time)
                x_limits = [min(minimum(x_vals), x_limits[1]), max(maximum(x_vals), x_limits[2])]

                for var in group
                    y_vals = getindex.(sol.u, var)
                    y_limits = [min(minimum(y_vals), y_limits[1]), max(maximum(y_vals), y_limits[2])]
                end

                plot!(p, sol, vars=group, label="Run $(findfirst(isequal(sol), results.solutions))")
            end

            # Apply dynamic limits with padding for both x and y axes in 2D plot
            xlims!(p, (x_limits[1] - 0.05 * (x_limits[2] - x_limits[1]), x_limits[2] + 0.05 * (x_limits[2] - x_limits[1])))
            ylims!(p, (y_limits[1] - 0.1 * (y_limits[2] - y_limits[1]), y_limits[2] + 0.1 * (y_limits[2] - y_limits[1])))

            push!(plots, p)
        end
    end
    
    for p in plots
        display(p)
    end
end