using LinearAlgebra
using Plots
using DifferentialEquations
using OrdinaryDiffEq, ProgressLogging
using Sundials

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
struct Problem
    A::Function
    B::Function
    phi::Function
    q0::Vector{Float64}
    qd0::Vector{Float64}
end


# Function to perform the simulation
function run_simulation(params::SimulationParams,problem::Problem)
    n = length(params.t)
    nm = 3
    nce = 2
    ne = nm,nce
    q = zeros(5, n)
    q_v = zeros(5, n)
    q_ac = zeros(5, n)
    PosConstrNorm = zeros(n)
    q[:, 1] = problem.q0
    q_v[:, 1] = problem.qd0
    PosConstrNorm[1] = problem.phi(q[:, 1])

    for i in 1:n-1

        Qa = [0, 0, -params.m1 * params.g]
        A = problem.A(q[:, i], params.m1)
        B = problem.B(q_v[:, i], Qa)
        temp = A \ B
        q_ac[:, i] = temp[1:5]
        if i == 1
            q_v[:, i+1] = q_v[:, i] + params.dt * q_ac[:, i]
            q[:, i+1] = q[:, i] + params.dt * q_v[:, i]
        else
            q_v[:, i+1] = q_v[:, i] + params.dt * (q_ac[:, i-1] + q_ac[:, i]) / 2
            q[:, i+1] = q[:, i] + params.dt * (q_v[:, i-1] + q_v[:, i]) / 2
        end

        PosConstrNorm[i+1] = problem.phi(q[:, i+1])
    end

    return q, PosConstrNorm
end

# Function to perform the simulation
function run_simulation2(params::SimulationParams,problem::Problem)
    n = length(params.t)
    nm = 3
    nce = 2
    ne = nm,nce
    q = zeros(5, n)
    q_v = zeros(5, n)
    q_ac = zeros(5, n)
    PosConstrNorm = zeros(n)
    q[:, 1] = problem.q0
    q_v[:, 1] = problem.qd0
    PosConstrNorm[1] = problem.phi(q[:, 1])

    for i in 1:n-1

        Qa = [0, 0, -params.m1 * params.g]
        A = problem.A(q[:, i], params.m1)
        B = problem.B(q_v[:, i], Qa)
        temp = A \ B
        q_ac[:, i] = temp[1:5]
        if i == 1
            q_v[:, i+1] = q_v[:, i] + params.dt * q_ac[:, i]
            q[:, i+1] = q[:, i] + params.dt * q_v[:, i]
        else
            q_v[:, i+1] = q_v[:, i] + params.dt * (q_ac[:, i-1] + q_ac[:, i]) / 2
            q[:, i+1] = q[:, i] + params.dt * (q_v[:, i-1] + q_v[:, i]) / 2
        end

        PosConstrNorm[i+1] = problem.phi(q[:, i+1])
    end

    return q, PosConstrNorm
end




# tspan = (0.0, 1.0)
# "₁₂₃₄₅₆₇₈₉₀"
# Function to plot the results
function plot_results(t, q, PosConstrNorm)
    gr()
    plot1 = plot(t, [q[1, :], q[2, :], q[3, :]], lw=2, ls=:solid, label=["x" "y" "z"])
    display(plot1)
    plot2 = plot(t, PosConstrNorm, lw=3, ls=:solid, color=:blue, label="PosConstrNorm")
    display(plot2)
end

function init_params()
    dt=0.001
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
function calculate_phi(q)
    return (q[1]^2 + q[2]^2 + q[3]^2) / 2
end

# Function to calculate the A matrix
function calculate_A(q, m1)
    return [m1 0 0 q[1];
    0 m1 0 q[2] ;
    0 0 m1 q[3];
    q[1] q[2] q[3] 0]
    # return [m1 0 0 q[1] 0;
    #         0 m1 0 q[2] q[2];
    #         0 0 m1 q[3] 0;
    #         q[1] q[2] q[3] 0 0;
    #         0 q[2] 0 0 0]
end

# Function to calculate B vector
function calculate_B(q_v, Qa)
    garma = [-q_v[1]^2 - q_v[2]^2 - q_v[3]^2]
    return vcat(Qa, garma)
    garma = [-q_v[1]^2 - q_v[2]^2 - q_v[3]^2, -q_v[2]^2]
    return vcat(Qa, garma)
end
function init_problem()
    q0 = [1, 1, 0, 0, 0]
    qd0 = zeros(5)
    return Problem(
        calculate_A,
        calculate_B,
        calculate_phi,
        q0,
        qd0
    )
end

# Main function to run the simulation
function main()
    params = init_params()
    problem = init_problem()
    q, PosConstrNorm = run_simulation(params,problem)
    plot_results(params.t, q, PosConstrNorm)
end

#main()


function f2(out, du, u, p, t)
    dq₁, dq₂, dq₃,  dv₁, dv₂, dv₃,dl₁,dl₂ = du
    q₁, q₂, q₃,     v₁, v₂, v₃,l₁,l₂ = u
    params=p[1]
    problem=p[2]
    m1=params.m1 
    Qa = [0, 0, -m1 * params.g]


    ddq=[dv₁; dv₂; dv₃]
    q=[q₁; q₂; q₃]
    dq=[ v₁;v₂; v₃]
    l=[l₁;l₂]

    A = problem.A(q', params.m1)
    B = problem.B(dq', Qa)

    res=A*vcat(ddq,l)-B
    out[1:3]=res[1:3]
    out[7:8] = res[4:5]

    out[4] = v₁ - dq₁
    out[5] = v₂ - dq₂
    out[6] = v₃ - dq₃

    params=0
end

function function_pendulum(out, du, u, p, t)
    dq₁, dq₂, dq₃,dl₁,  dv₁, dv₂, dv₃ = du
    q₁ , q₂ , q₃ , l₁,   v₁,  v₂,  v₃ = u
    params=p[1]
    problem=p[2]
    m1=params.m1 
    Qa = [0, 0, -m1 * params.g]

    ddq=[dv₁; dv₂; dv₃]
    q=[q₁; q₂; q₃]
    dq=[ v₁;v₂; v₃]
    l=[l₁]

    A = problem.A(q', params.m1)
    B = problem.B(dq', Qa)

    res=A*vcat(ddq,l)-B
    out[1:3]=res[1:3]
    out[4] = res[4]

    out[5] = v₁ - dq₁
    out[6] = v₂ - dq₂
    out[7] = v₃ - dq₃

end

function f3(out, du, u, p, t)
    dq₁,dq₂,dv₁,dv₂, dl₁ = du
    q₁,q₂,v₁,v₂, l₁ = u


    out[1] = dv₁-9.81+l₁*q₁
    out[2] = dv₂+9.81
    out[3] = v₁-dq₁
    out[4] = v₂-dq₂
    out[5] = q₁*q₁-1

    params=0
end

function f4(out, du, u, p, t)
    dq₁,dv₁, dl₁ = du
    q₁,v₁, l₁ = u

    out[1] = dv₁+9.81+l₁

    out[2] = v₁-dq₁

    out[3] = q₁-1

end

function f5(out, du, u, p, t)
    dq₁,dv₁ = du
    q₁,v₁ = u

    out[1] = dv₁+9.81

    out[2] = v₁-dq₁


    params=0
end

function test(du₀, u₀, out, p, t, differential_vars, equation)
    #equation(out, du₀, u₀, p, t)
    #println(out)
    tspan = (0.0, 1.0)

    prob = DAEProblem(equation, du₀, u₀, tspan, p, differential_vars=differential_vars)

    sol = solve(prob, IDA(), reltol=1e-9, abstol=1e-9, progress=true)

    # Create a plot
    p = plot(title="Solution to the linear ODE with a thick line",
        xaxis="Time (t)", yaxis="u(t) (in μm)")

    # Plot only the first 5 curves
    for i in 1:min(3, length(sol))
        plot!(p, sol.t, sol[i, :], linewidth=3, label="Curve $i")
    end

    # Display the plot
    display(p)
    return sol
end
#p=Any[init_params(),init_prolbem()]
p = (init_params(), init_problem())
t=0.0


differential_vars = [true, true,true, true,true, true, false]
u₀ = [1.0, 0.0, 0.0,    0.0,    0.0, 0.0,0.0]
du₀ = [0.0, 0.0, 0,     0.0,    0.0,0.0, -9.81]
out=[0.0, 0.0, 0,0.0, 0.0, 0,0.0]
sol=test(du₀, u₀, out, p, t,differential_vars,function_pendulum)

# differential_vars = [true, true,true, true,true, true, false, false]
# u₀ = [1.0, 0.0, 0,0.0, 0.0, 0,0.0, 0.0]
# du₀ = [0.0, 0.0, 0,0.0, 0.0, -9.81,0.0, 0.0]
# out=zeros(8)
# test(du₀, u₀, out, p, t,differential_vars,f2)

# differential_vars = [true, true,true, true, false]
# du₀ = [0.0  , 0.    ,   0,  -9.81,  0]
# u₀ = [1.0   , 1.0    ,   0.0,    0., 9.81]
# out=zeros(5)
# test(du₀, u₀, out, p, t,differential_vars,f3)

# differential_vars = [true, true, false]
# du₀ = [0.0  ,    0.0,    0.0]
# u₀ = [1.0   ,   0.0,     -9.81]
# out=zeros(3)
# test(du₀, u₀, out, p, t,differential_vars,f4)

# differential_vars = [true, true]
# du₀ = [0.0  ,    -9.81]
# u₀ = [1.0   ,   0.0]
# out=zeros(2)
# test(du₀, u₀, out, p, t,differential_vars,f5)