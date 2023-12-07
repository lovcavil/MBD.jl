using LinearAlgebra
using Plots

# Define a struct for simulation parameters
struct SimulationParams
    m1::Float64
    g::Float64
    l1::Float64
    l2::Float64
    w::Float64
    dt::Float64
    t::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
end

# Function to initialize simulation parameters
function init_params()
    dt=0.000001
    return SimulationParams(
        30.0,      # m1
        9.81,      # g
        50.0,      # l1
        220.0,     # l2
        5 * π,     # w
        dt,     # dt
        0:dt:5  # t
    )
end

# Function to calculate the constraint function phi
function phi(q)
    return (q[1]^2 + q[2]^2 + q[3]^2) / 2
end

# Function to calculate the A matrix
function calculate_A(q, m1)
    return [m1 0 0 q[1] 0;
            0 m1 0 q[2] q[2];
            0 0 m1 q[3] 0;
            q[1] q[2] q[3] 0 0;
            0 q[2] 0 0 0]
end

# Function to calculate B vector
function calculate_B(q_v, Qa)
    garma = [-q_v[1]^2 - q_v[2]^2 - q_v[3]^2, -q_v[2]^2]
    return vcat(Qa, garma)
end

# Function to perform the simulation
function run_simulation(params::SimulationParams)
    n = length(params.t)
    q = zeros(5, n)
    q_v = zeros(5, n)
    q_ac = zeros(5, n)
    PosConstrNorm = zeros(n)
    q[:, 1] = [1, 1, 0, 0, 0]
    q_v[:, 1] = zeros(5)
    PosConstrNorm[1] = phi(q[:, 1])

    for i in 1:n-1
        A = calculate_A(q[:, i], params.m1)
        Qa = [0, 0, -params.m1 * params.g]
        B = calculate_B(q_v[:, i], Qa)
        temp = A \ B
        q_ac[:, i] = temp[1:5]
        if i == 1
            q_v[:, i+1] = q_v[:, i] + params.dt * q_ac[:, i]
            q[:, i+1] = q[:, i] + params.dt * q_v[:, i]
        else
            q_v[:, i+1] = q_v[:, i] + params.dt * (q_ac[:, i-1] + q_ac[:, i]) / 2
            q[:, i+1] = q[:, i] + params.dt * (q_v[:, i-1] + q_v[:, i]) / 2
        end

        PosConstrNorm[i+1] = phi(q[:, i+1])
    end

    return q, PosConstrNorm
end

# Function to plot the results
function plot_results(t, q, PosConstrNorm)
    gr()
    plot1 = plot(t, [q[1, :], q[2, :], q[3, :]], lw=2, ls=:solid, label=["x" "y" "z"])
    display(plot1)
    plot2 = plot(t, PosConstrNorm, lw=3, ls=:solid, color=:blue, label="PosConstrNorm")
    display(plot2)
end

# Main function to run the simulation
function main()
    params = init_params()
    q, PosConstrNorm = run_simulation(params)
    plot_results(params.t, q, PosConstrNorm)
end

main()
