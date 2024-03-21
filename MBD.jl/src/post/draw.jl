# Adjust the draw function to handle multiple solutions and dynamically set figure limits
function draw_0(results::ODERunResults; groups, plot3d::Bool=true)
    colors = ["#012a4a", "#013a63", "#01497c", "#014f86", "#2a6f97", "#2c7da0", "#468faf", "#61a5c2", "#89c2d9", "#a9d6e5"]
    colors1 = ["#582f0e", "#7f4f24", "#936639", "#a68a64", "#b6ad90", "#c2c5aa", "#a4ac86", "#656d4a", "#414833", "#333d29"]
    colors1 = ["#5f0f40", "#9a031e", "#fb8b24", "#e36414", "#0f4c5c"]
    plots = []
    line_styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]

    # Map each variable 'i' to a line style
    var_to_linestyle = Dict()
    for group in groups
        unique_vars = unique(group)
        for (i, var) in enumerate(unique_vars)
            var_to_linestyle[var] = line_styles[i%length(line_styles)]
        end
        if plot3d && length(group) == 3
            # Initialize limits for 3D plot
            limits = Dict("x" => [Inf, -Inf], "y" => [Inf, -Inf], "z" => [Inf, -Inf])

            # 3D plot creation
            p = plot(layout=(1, 1), legend=:topright)

            for (index, sol) in enumerate(results.solutions)
                x_vals = getindex.(sol.u, group[1])
                y_vals = getindex.(sol.u, group[2])
                z_vals = getindex.(sol.u, group[3])

                # Update limits
                limits["x"] = [min(minimum(x_vals), limits["x"][1]), max(maximum(x_vals), limits["x"][2])]
                limits["y"] = [min(minimum(y_vals), limits["y"][1]), max(maximum(y_vals), limits["y"][2])]
                limits["z"] = [min(minimum(z_vals), limits["z"][1]), max(maximum(z_vals), limits["z"][2])]

                plot!(p, x_vals, y_vals, z_vals, label="Run$index", xlabel="X-axis label", ylabel="Y-axis label")
            end

            # Apply limits with padding
            xlims!(p, (limits["x"][1], limits["x"][2]))
            ylims!(p, (limits["y"][1], limits["y"][2]))
            zlims!(p, (limits["z"][1], limits["z"][2]))

            push!(plots, p)
        else
            # Adjusted 2D plot code to include line style changes
            x_limits = [Inf, -Inf]
            y_limits = [Inf, -Inf]

            p = plot(legend=:topright)

            for (index, sol) in enumerate(results.solutions)
                x_vals = sol.t
                x_limits = [min(minimum(x_vals), x_limits[1]), max(maximum(x_vals), x_limits[2])]

                for var in group
                    y_vals = getindex.(sol.u, var)
                    y_limits = [min(minimum(y_vals), y_limits[1]), max(maximum(y_vals), y_limits[2])]

                    line_style = var_to_linestyle[var]
                    plot!(p, sol.t, y_vals, color=colors1[index], linewidth=1.5, label="Run$index, q$var", linestyle=line_style)
                end
            end

            xlims!(p, (x_limits[1] - 0.05 * (x_limits[2] - x_limits[1]), x_limits[2] + 0.05 * (x_limits[2] - x_limits[1])))
            ylims!(p, (y_limits[1] - 0.1 * (y_limits[2] - y_limits[1]), y_limits[2] + 0.1 * (y_limits[2] - y_limits[1])))

            push!(plots, p)
        end
    end

    for p in plots
        display(p)
    end
end


using DataFrames
function draw(dfs; groups, plot3d::Bool=true)
    colors = ["#012a4a", "#013a63", "#01497c", "#014f86", "#2a6f97", "#2c7da0", "#468faf", "#61a5c2", "#89c2d9", "#a9d6e5"]
    colors1 = ["#5f0f40", "#9a031e", "#fb8b24", "#e36414", "#0f4c5c"]
    plots = []
    line_styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    p0 = plot(layout=(1, 1), legend=:topright, projection=:perspective)
    for group in groups
        if plot3d && length(group) == 3
            

            # Initialize limits
            limits = Dict("x" => [Inf, -Inf], "y" => [Inf, -Inf], "z" => [Inf, -Inf])

            for df in dfs
                x_vals, y_vals, z_vals = df[:, Symbol(group[1])], df[:, Symbol(group[2])], df[:, Symbol(group[3])]
                # Update limits
                limits["x"] = [min(minimum(x_vals), limits["x"][1]), max(maximum(x_vals), limits["x"][2])]
                limits["y"] = [min(minimum(y_vals), limits["y"][1]), max(maximum(y_vals), limits["y"][2])]
                limits["z"] = [min(minimum(z_vals), limits["z"][1]), max(maximum(z_vals), limits["z"][2])]
            end

            # Plot with adjusted limits
            for (run_index, df) in enumerate(dfs)
                plot!(p0, df[:, Symbol(group[1])], df[:, Symbol(group[2])], df[:, Symbol(group[3])], 
                      label="Run$run_index", color=colors1[run_index % length(colors1) + 1], xlims=(limits["x"][1], limits["x"][2]), 
                      ylims=(limits["y"][1], limits["y"][2]), zlims=(limits["z"][1], limits["z"][2]))
            end

            push!(plots, p0)
        else
            # Adjusted 2D plot code with dynamic limit adjustments
            x_limits, y_limits = [Inf, -Inf], [Inf, -Inf]

            p = plot(legend=:topright)

            for df in dfs
                for var in group
                    y_vals = df[:, Symbol(var)]
                    # Update limits
                    y_limits = [min(minimum(y_vals), y_limits[1]), max(maximum(y_vals), y_limits[2])]
                end
            end

            # Plot with adjusted limits
            for (run_index, df) in enumerate(dfs)
                for (var_index, var) in enumerate(group)
                    plot!(p, df.Time, df[:, Symbol(var)], label="Run$run_index, $var", 
                          color=colors1[run_index % length(colors1) + 1], linestyle=line_styles[var_index % length(line_styles) + 1], 
                          xlims=(minimum(df.Time), maximum(df.Time)), ylims=(y_limits[1] - 0.1 * diff(y_limits)[1], y_limits[2] + 0.1 * diff(y_limits)[1]))
                end
            end

            push!(plots, p)
        end
    end

# Assuming 'plots' is your array of plot objects
for (i, p) in enumerate(plots)
    display(p)  # Display the plot in the REPL or a notebook
    
    # Construct a file name for each plot
    filename = "plot_$(i).png"  # or .svg, .pdf, etc., depending on your preferred format
    
    # Save the plot to disk
    savefig(p, filename)
end
end
function draw3d(dfs; groups, plot3d::Bool=true, limits = Dict("x" => [4000, -4000], "y" => [1000, -1000], "z" => [3000, -1000]))
    colors = ["#012a4a", "#013a63", "#01497c", "#014f86", "#2a6f97", "#2c7da0", "#468faf", "#61a5c2", "#89c2d9", "#a9d6e5"]
    colors1 = ["#5f0f40", "#9a031e", "#fb8b24", "#e36414", "#0f4c5c"]
    plots = []
    line_styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    p0 = plot(layout=(1, 1), legend=:topright, projection=:perspective)
    for group in groups
            # Initialize limits
            #limits = Dict("x" => [4000, -4000], "y" => [1000, -1000], "z" => [3000, -1000])

            # for df in dfs
            #     x_vals, y_vals, z_vals = df[:, Symbol(group[1])], df[:, Symbol(group[2])], df[:, Symbol(group[3])]
            #     # Update limits
            #     limits["x"] = [min(minimum(x_vals), limits["x"][1]), max(maximum(x_vals), limits["x"][2])]
            #     limits["y"] = [min(minimum(y_vals), limits["y"][1]), max(maximum(y_vals), limits["y"][2])]
            #     limits["z"] = [min(minimum(z_vals), limits["z"][1]), max(maximum(z_vals), limits["z"][2])]
            # end

            # Plot with adjusted limits
            for (run_index, df) in enumerate(dfs)
                plot!(p0, df[:, Symbol(group[1])], df[:, Symbol(group[2])], df[:, Symbol(group[3])], 
                      label="Run$run_index", color=colors1[run_index % length(colors1) + 1], xlims=(limits["x"][2], limits["x"][1]), 
                      ylims=(limits["y"][2], limits["y"][1]), zlims=(limits["z"][2], limits["z"][1]))
            end

            

    end
    push!(plots, p0)
    for p in plots
        display(p)
    end
end