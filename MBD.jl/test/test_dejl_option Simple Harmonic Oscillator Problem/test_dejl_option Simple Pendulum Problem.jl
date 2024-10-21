# Simple Pendulum Problem
using OrdinaryDiffEq, Plots

#Constants
const g = 9.81
L = 1.0

#Initial Conditions
u₀ = [0, π / 2]
tspan = (0.0, 6.3)

#Define the problem
function simplependulum(du, u, p, t)
    θ = u[1]
    dθ = u[2]
    du[1] = dθ
    du[2] = -(g / L) * sin(θ)
end

#Pass to solvers
prob = ODEProblem(simplependulum, u₀, tspan)
sol = solve(prob, Tsit5())

#Plot
plot(sol, linewidth = 2, title = "Simple Pendulum Problem", xaxis = "Time",
    yaxis = "Height", label = ["\\theta" "d\\theta"])

p = plot(sol, vars = (1, 2), xlims = (-9, 9), title = "Phase Space Plot",
    xaxis = "Angular position", yaxis = "Angular velocity", leg = false)
function phase_plot(prob, u0, p, tspan = 2pi)
    _prob = ODEProblem(prob.f, u0, (0.0, tspan))
    sol = solve(_prob, Vern9()) # Use Vern9 solver for higher accuracy
    plot!(p, sol, idxs = (1, 2), xlims = nothing, ylims = nothing)
end
for i in (-4pi):(pi / 2):(4π)
    for j in (-4pi):(pi / 2):(4π)
        phase_plot(prob, [j, i], p)
    end
end
plot(p, xlims = (-9, 9))