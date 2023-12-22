using DifferentialEquations
using Sundials
using Plots
# Define the DAE system
function dae_system(out, du, u, p, t)
    x, y, z = u
    dx, dy, dz = du

    # Equations
    out[1] = dx+1  # Corresponding to du[1] = x - y - z
    out[2] = dy+z   # Corresponding to du[2] = x + y + z - 1
    out[3] = z-sin(t)   # Corresponding to du[3] = x^2 + y^2 - z
end

# Initial conditions
u0 = [0.0, 0.0, 0.0]
du0 = [0.0, 0.0, 0]  # Initial derivatives (guessed or calculated)
tspan = (0.0, 10.0)
v=[true,true,false]
# Create a DAE problem
dae_problem = DAEProblem(dae_system, du0, u0, tspan,differential_vars=v)

# Solve the DAE with different solvers
sol1 = solve(dae_problem, IDA())
sol2 = solve(dae_problem, DFBDF())
#sol3 = solve(dae_problem, TRBDF2())
#sol4 = solve(dae_problem, DABDF2())
#sol5 = solve(dae_problem, Rosenbrock23())

# Display the solutions
#println("Solution with IDA: ", sol1)
p1=plot(sol1,title="Solution with IDA")
display(p1)
#println("Solution with DFBDF(): ", sol2)
p2=plot(sol2,title="Solution with DFBDF")
display(p2)
#println("Solution with TRBDF2: ", sol3)
#println("Solution with DABDF2: ", sol4)
#println("Solution with Rosenbrock23: ", sol5)
