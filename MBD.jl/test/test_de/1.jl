using DifferentialEquations

# Define a simple ODE problem
function model!(du, u, p, t)
    du[1] = p[1] * u[1] * (1 - u[1] / p[2])
end

u0 = [1.0]
tspan = (0.0, 10.0)
params = [0.1, 10.0]
prob = ODEProblem(model!, u0, tspan, params)

# Example of solving with VCABM and setting options
sol = solve(prob, VCABM(autodiff=true, linearsolver=:default))
