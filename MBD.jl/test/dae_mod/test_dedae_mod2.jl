using LinearAlgebra
using Plots
using DifferentialEquations
using OrdinaryDiffEq, ProgressLogging
using Sundials

# Define a struct for simulation parameters
struct SimulationParams
    m1::Float64
    g::Float64
    dt::Float64
    t::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}
    du₀::Vector{Float64}
    u₀:: Vector{Float64}
    differential_vars::Vector{Bool}
    equation::Function
end
struct Problem
    #phi::Function
end


# tspan = (0.0, 1.0)
# "₁₂₃₄₅₆₇₈₉₀"


# Function to calculate the constraint function phi
function calculate_phi(q)
    return (q[1]^2 + q[2]^2 + q[3]^2) / 2
end

function rhs(x,m1,q,q_v,g)
    A=[m1 0 0 q[1];
       0 m1 0 q[2] ;
       0 0 m1 q[3];
       q[1] q[2] q[3] 0]
    Qa = [0, 0, -m1 * params.g]
    b=vcat(Qa, -q_v[1]^2 - q_v[2]^2 - q_v[3]^2)
    return A*x-b
end




function function_pendulum(out, du, u, params, t)
    dq₁, dq₂, dq₃,dl₁,  dv₁, dv₂, dv₃ = du
    q₁ , q₂ , q₃ , l₁,   v₁,  v₂,  v₃ = u

    m1=params.m1 
    g=params.g
    ddq=[dv₁; dv₂; dv₃]
    dq=[ v₁;v₂; v₃]
    q=[q₁; q₂; q₃]
    l=[l₁]

    x=vcat(ddq,l)
    res=rhs(x,m1,q,dq,g)
    out[1:4]=res[1:4]

    out[5] = (v₁ - dq₁)
    out[6] = (v₂ - dq₂)
    out[7] = (v₃ - dq₃)

end


function run(out, p, t)
    #equation(out, du₀, u₀, p, t)
    #println(out)
    tspan = (0.0, 1.0)
    du₀, u₀=p.du₀,p.u₀ 
    prob = DAEProblem(p.equation, du₀, u₀, tspan, p, differential_vars=p.differential_vars)

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

function run_with_ode_solver(p, t)
    #equation(out, du₀, u₀, p, t)
    #println(out)
    tspan = (0.0, 1.0)
    u₀=p.u₀ 
    prob = ODEProblem(odequation,  u₀, tspan, p)

    sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-9, progress=true)

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

function init_params(du0,u0,differential_vars,function_pendulum)
    dt=0.001
    return SimulationParams(
        30.0,      # m1
        9.81,      # g
        dt,     # dt
        0:dt:5,  # t
        du0,
        u0,
        differential_vars,
        function_pendulum
    )
end

t=0.0
println("START")
u₀  = [1.0, 0.0, 0.0,    0.0,    0.0, 0.0,0.0]
du₀ = [0.0, 0.0, 0,     0.0,    0.0,0.0, -9.81]
out = [0.0, 0.0, 0,0.0, 0.0, 0,0.0]
differential_vars= [true, true,true, true,true, true, false]
params=init_params(du₀, u₀,differential_vars,function_pendulum)
sol=run(out, params, t)


