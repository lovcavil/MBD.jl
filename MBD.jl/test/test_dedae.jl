using DifferentialEquations
using Plots
using OrdinaryDiffEq, ProgressLogging

# function lorenz!(du, u, p, t)
#     du[1] = 10.0(u[2] - u[1])
#     du[2] = u[1] * (28.0 - u[3]) - u[2]
#     du[3] = u[1] * u[2] - (8 / 3) * u[3]
# end
# u0 = [1.0; 0.0; 0.0]
# tspan = (0.0, 100000.0)
# prob = ODEProblem(lorenz!, u0, tspan)
# sol = solve(prob, Tsit5(), progress = true)

function rober(du, u, p, t)
    y₁, y₂, y₃ = u
    k₁, k₂, k₃ = p
    du[1] = -k₁ * y₁ + k₃ * y₂ * y₃
    du[2] = k₁ * y₁ - k₃ * y₂ * y₃ - k₂ * y₂^2
    du[3] = y₁ + y₂ + y₃ - 1
    nothing
    
end
M = [1.0 0 0
    0 1.0 0
    0 0 0]
f = ODEFunction(rober, mass_matrix = M)
prob_mm = ODEProblem(f, [1.0, 0.0, 0.0], (0.0, 1e6), (0.04, 3e7, 1e4))
sol = solve(prob_mm, Rodas5(), reltol = 1e-9, abstol = 1e-9,progress = true)
println("ee")

plot(sol, linewidth = 5, title = "Solution to the linear ODE with a thick line",
    xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!") # legend=false
#plot!(sol.t, t -> 0.5 * exp(1.01t), lw = 3, ls = :dash, label = "True Solution!")

