# HH_system.jl
using OrdinaryDiffEq

# Define the potential function V(x, y)
function V(x, y)
    return 1//2 * (x^2 + y^2 + 2x^2 * y - 2//3 * y^3)
end

# Define the energy function E(x, y, dx, dy)
function E(x, y, dx, dy)
    return V(x, y) + 1//2 * (dx^2 + dy^2)
end

# Convert the second-order ODE to a first-order ODE system
function HH_first_order!(du, u, p, t)
    x, y, dx, dy = u
    du[1] = dx
    du[2] = dy
    du[3] = -x - 2x * y
    du[4] = -y - x^2 + y^2
end
