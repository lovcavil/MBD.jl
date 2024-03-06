# Adjust the draw function to handle multiple solutions and dynamically set figure limits
function draw(results::ODERunResults; groups)
    plots = []
    for group in groups
        p = plot()

        # Determine dynamic limits
        x_min, x_max, y_min, y_max = Inf, -Inf, Inf, -Inf
        for sol in results.solutions
            for var in group
                # Assuming `sol` is a structure where you can access the time (t) and solution (u) vectors
                x_vals = sol.t
                y_vals = getindex.(sol.u, var)

                # Update x and y limits
                x_min = min(x_min, minimum(x_vals))
                x_max = max(x_max, maximum(x_vals))
                y_min = min(y_min, minimum(y_vals))
                y_max = max(y_max, maximum(y_vals))
            end

            plot!(p, sol, vars=group, label="Run $(findfirst(isequal(sol), results.solutions))")
        end

        # Apply dynamic limits with some padding
        xlims!(p, (x_min, x_max))
        ylims!(p, (y_min - 0.1 * (y_max - y_min), y_max + 0.1 * (y_max - y_min)))

        push!(plots, p)
    end
    for p in plots
        display(p)
    end
end