using DifferentialEquations
using Plots
using OrdinaryDiffEq, ProgressLogging
using Sundials




function f2(out, du, u, p, t)
    out[1] = -0.04u[1] + 1e4 * u[2] * u[3] - du[1]
    out[2] = +0.04u[1] - 3e7 * u[2]^2 - 1e4 * u[2] * u[3] - du[2]
    out[3] = u[1] + u[2] + u[3] - 1.0
end

u₀ = [1.0, 0, 0]
du₀ = [-0.04, 0.04, 0.0]
tspan = (0.0, 100000.0)

differential_vars = [true, true, false]
prob = DAEProblem(f2, du₀, u₀, tspan, differential_vars = differential_vars)



sol = solve(prob, IDA(), reltol = 1e-9, abstol = 1e-9,progress = true)

println("ee")

plot(sol, linewidth = 5, title = "Solution to the linear ODE with a thick line",
    xaxis = "Time (t)", yaxis = "u(t) (in μm)", label = "My Thick Line!") # legend=false
#plot!(sol.t, t -> 0.5 * exp(1.01t), lw = 3, ls = :dash, label = "True Solution!")

