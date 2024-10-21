#Double Pendulum Problem
using OrdinaryDiffEq, Plots

#Constants and setup
const m₁, m₂, L₁, L₂ = 1, 2, 1, 2
initial = [0, π / 3, 0, 3pi / 5]
tspan = (0.0, 50.0)

#Convenience function for transforming from polar to Cartesian coordinates
function polar2cart(sol; dt = 0.02, l1 = L₁, l2 = L₂, vars = (2, 4))
    u = sol.t[1]:dt:sol.t[end]

    p1 = l1 * map(x -> x[vars[1]], sol.(u))
    p2 = l2 * map(y -> y[vars[2]], sol.(u))

    x1 = l1 * sin.(p1)
    y1 = l1 * -cos.(p1)
    (u, (x1 + l2 * sin.(p2),
        y1 - l2 * cos.(p2)))
end

#Define the Problem
function double_pendulum(xdot, x, p, t)
    xdot[1] = x[2]
    xdot[2] = -((g * (2 * m₁ + m₂) * sin(x[1]) +
                 m₂ * (g * sin(x[1] - 2 * x[3]) +
                  2 * (L₂ * x[4]^2 + L₁ * x[2]^2 * cos(x[1] - x[3])) * sin(x[1] - x[3]))) /
                (2 * L₁ * (m₁ + m₂ - m₂ * cos(x[1] - x[3])^2)))
    xdot[3] = x[4]
    xdot[4] = (((m₁ + m₂) * (L₁ * x[2]^2 + g * cos(x[1])) +
                L₂ * m₂ * x[4]^2 * cos(x[1] - x[3])) * sin(x[1] - x[3])) /
              (L₂ * (m₁ + m₂ - m₂ * cos(x[1] - x[3])^2))
end

#Pass to Solvers
double_pendulum_problem = ODEProblem(double_pendulum, initial, tspan)
sol = solve(double_pendulum_problem, Vern7(), abstol = 1e-10, dt = 0.05);

ts, ps = polar2cart(sol, l1 = L₁, l2 = L₂, dt = 0.01)
plot(ps...)