#辛几何
using OrdinaryDiffEq, Plots
V(x, y) = 1 // 2 * (x^2 + y^2 + 2x^2 * y - 2 // 3 * y^3)
E(x, y, dx, dy) = V(x, y) + 1 // 2 * (dx^2 + dy^2);
function HH_acceleration!(dv, v, u, p, t)
  x, y = u
  dx, dy = dv
  dv[1] = -x - 2x * y
  dv[2] = y^2 - y - x^2
end
initial_positions = [0.0, 0.1]
initial_velocities = [0.5, 0.0]
prob = SecondOrderODEProblem(HH_acceleration!, initial_velocities, initial_positions, tspan)
sol2 = solve(prob, KahanLi8(), dt = 1 / 10);
# Plot the orbit
p1=plot(sol2, idxs = (3, 4), title = "The orbit of the Hénon-Heiles system", xaxis = "x",
    yaxis = "y", leg = false)
display(p1) 
# Plot the orbit
p2=plot(sol2, idxs = (3, 1), title = "Phase space for the Hénon-Heiles system",
    xaxis = "Position", yaxis = "Velocity")
plot!(sol2, idxs = (4, 2), leg = false)
display(p2) 
energy = map(x -> E(x[3], x[4], x[1], x[2]), sol2.u)
#We use @show here to easily spot erratic behaviour in our system by seeing if the loss in energy was too great.
@show ΔE = energy[1] - energy[end]

#Plot
p3=plot(sol2.t, energy .- energy[1], title = "Change in Energy over Time",
    xaxis = "Time in iterations", yaxis = "Change in Energy")
display(p3) 

sol3 = solve(prob, DPRKN6());
energy = map(x -> E(x[3], x[4], x[1], x[2]), sol3.u)
@show ΔE = energy[1] - energy[end]
gr()
p4=plot(sol3.t, energy .- energy[1], title = "Change in Energy over Time",
    xaxis = "Time in iterations", yaxis = "Change in Energy")
display(p4) 
    