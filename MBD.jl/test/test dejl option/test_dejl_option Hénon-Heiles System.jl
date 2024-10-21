#Origin

using OrdinaryDiffEq, Plots

#Setup
initial = [0.0, 0.1, 0.5, 0]
tspan = (0, 100.0)

#Remember, V is the potential of the system and T is the Total Kinetic Energy, thus E will
#the total energy of the system.
V(x, y) = 1 // 2 * (x^2 + y^2 + 2x^2 * y - 2 // 3 * y^3)
E(x, y, dx, dy) = V(x, y) + 1 // 2 * (dx^2 + dy^2);

#Define the function
function Hénon_Heiles(du, u, p, t)
    x = u[1]
    y = u[2]
    dx = u[3]
    dy = u[4]
    du[1] = dx
    du[2] = dy
    du[3] = -x - 2x * y
    du[4] = y^2 - y - x^2
end

#Pass to solvers
prob = ODEProblem(Hénon_Heiles, initial, tspan)
sol = solve(prob, Vern9(), abstol = 1e-16, reltol = 1e-16);
# Plot the orbit
p1=plot(sol, idxs = (1, 2), title = "The orbit of the Hénon-Heiles system", xaxis = "x",
    yaxis = "y", leg = false)
display(p1) 
#Optional Sanity check - what do you think this returns and why?
@show sol.retcode

#Plot -
p2=plot(sol, idxs = (1, 3), title = "Phase space for the Hénon-Heiles system",
    xaxis = "Position", yaxis = "Velocity")
plot!(sol, idxs = (2, 4), leg = false)
display(p2) 

#We map the Total energies during the time intervals of the solution (sol.u here) to a new vector
#pass it to the plotter a bit more conveniently
energy = map(x -> E(x...), sol.u)

#We use @show here to easily spot erratic behavior in our system by seeing if the loss in energy was too great.
@show ΔE = energy[1] - energy[end]

#Plot
p3=plot(sol.t, energy .- energy[1], title = "Change in Energy over Time",
    xaxis = "Time in iterations", yaxis = "Change in Energy")
display(p3) 
