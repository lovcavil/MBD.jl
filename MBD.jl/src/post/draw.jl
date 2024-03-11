# Adjust the draw function to handle multiple solutions and dynamically set figure limits
function draw(results::ODERunResults; groups, plot3d::Bool=true)
    colors = ["#012a4a", "#013a63", "#01497c", "#014f86", "#2a6f97", "#2c7da0", "#468faf", "#61a5c2", "#89c2d9", "#a9d6e5"]
    colors1 = ["#582f0e","#7f4f24","#936639","#a68a64","#b6ad90","#c2c5aa","#a4ac86","#656d4a","#414833","#333d29"]
    colors1 = ["#5f0f40","#9a031e","#fb8b24","#e36414","#0f4c5c"]
    plots = []
    line_styles = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    
    # Map each variable 'i' to a line style
    var_to_linestyle = Dict()
    for group in groups
        unique_vars = unique(group)
        for (i, var) in enumerate(unique_vars)
            var_to_linestyle[var] = line_styles[i % length(line_styles) ]
        end
        if plot3d && length(group) == 3
            # Initialize limits for 3D plot
            limits = Dict("x" => [Inf, -Inf], "y" => [Inf, -Inf], "z" => [Inf, -Inf])
            
            # 3D plot creation
            p = plot(layout=(1,1), legend=:topright)

            for (index, sol) in enumerate(results.solutions)
                x_vals = getindex.(sol.u, group[1])
                y_vals = getindex.(sol.u, group[2])
                z_vals = getindex.(sol.u, group[3])

                # Update limits
                limits["x"] = [min(minimum(x_vals), limits["x"][1]), max(maximum(x_vals), limits["x"][2])]
                limits["y"] = [min(minimum(y_vals), limits["y"][1]), max(maximum(y_vals), limits["y"][2])]
                limits["z"] = [min(minimum(z_vals), limits["z"][1]), max(maximum(z_vals), limits["z"][2])]

                plot!(p, x_vals, y_vals, z_vals, label="Run$index")
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