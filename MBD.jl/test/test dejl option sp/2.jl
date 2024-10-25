using DifferentialEquations
using DataFrames
using CSV

# Example solving an ODE
function f(du, u, p, t)
    du .= -u
end

u0 = [1.0, 0.5]  # Example for multiple initial conditions
tspan = (0.0, 10.0)
prob = ODEProblem(f, u0, tspan)
sol = solve(prob)

# Convert the solution to a DataFrame
# Each row corresponds to a time step, and each column corresponds to a component of `u`
df = DataFrame(time = sol.t)

# Add each component of the solution `u` as a separate column
for (i, component) in enumerate(sol.u[1]) # Assuming sol.u has vectors
    df[!, "u$i"] = [u[i] for u in sol.u]   # Extract each component of `u` into its own column
end

# Save the DataFrame to CSV
CSV.write("solution.csv", df)
