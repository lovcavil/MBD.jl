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

function ode_rhs(m1,q,q_v,g)
    A=[m1 0 0           q[1] q[1];
       0 m1 0           q[2] 0;
       0 0 m1           q[3] 0;
       q[1] q[2] q[3]   0 0;
       q[1] 0 0         0 0]
    Qa = [0, 0, -m1 * params.g]
    b=vcat(Qa, -q_v[1]^2 - q_v[2]^2 - q_v[3]^2,-q_v[1]^2)
    return A\b
end

function odequation(du, u, params, t)
    q₁ , q₂ , q₃ , l₁,l₂,   v₁,  v₂,  v₃ = u

    m1=params.m1 
    g=params.g
    #ddq=#
    dq=[ v₁;v₂; v₃]
    q=[q₁; q₂; q₃]
    l=[l₁;l₂]

    #x=vcat(ddq,l)

    res=ode_rhs(m1,q,dq,g)
    du[6:8]=res[1:3]
    du[4:5]=res[4:5]
    du[1] = v₁ 
    du[2] = v₂ 
    du[3] = v₃

end

function run_with_ode_solver(p, t)
    #equation(out, du₀, u₀, p, t)
    #println(out)
    tspan = (0.0, 1.0)
    u₀=p.u₀ 
    prob = ODEProblem(odequation,  u₀, tspan, p)

    sol = solve(prob, Tsit5(), reltol=1e-15, abstol=1e-15, progress=true)

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

function init_params(du0,u0)
    dt=0.001
    return SimulationParams(
        30.0,      # m1
        9.81,      # g
        dt,     # dt
        0:dt:5,  # t
        du0,
        u0
    )
end

t=0.0
println("START")
u₀  = [1.0, 0.0, 0.0,    0.0,0.0,    0.0, 0.0,0.0]
du₀ = [0.0, 0.0, 0,     0.0,0.0,    0.0,0.0, -9.81]

differential_vars= [true, true,true, true,true, true, false]

# params=init_params(du₀, u₀,differential_vars,function_pendulum)
# sol=run_with_ode_solver(params, t)

function run_with_ode_solver_stepwise(p, t)
    tspan = (0.0, 1)
    u₀ = p.u₀ 
    prob = ODEProblem(odequation, u₀, tspan, p)

    # Creating an integrator
    integrator = init(prob, Tsit5(), reltol=1e-5, abstol=1e-5)

    # Create a plot
    p = plot(title="Solution to the linear ODE with a thick line",
        xaxis="Time (t)", yaxis="u(t) (in μm)")
    PosConstrNorm=[]
    VelConstrNorm=[]
    t=[]
    while integrator.t < tspan[2]
        step!(integrator)  # Step through the solution

        q=integrator.u[1:3]
        qd=integrator.u[4:7]
        # # Calculate constraint error
        Phi = calculate_phi(q)
        Phiq = qd
        # Gam = mathfunction.GamEval(tn, q, qd, SJDT, par)
        # #println("3norm(Phiq * qd)",norm(Phiq * qd))
        push!(PosConstrNorm, norm(Phi))
        push!(VelConstrNorm, norm(Phiq .* qd))
        push!(t, integrator.t)

        # push!(AccConstrNorm, norm(Phiq * qdd + Gam))

        # Add data to the plot for the first 3 variables
        for i in 1:min(3, length(integrator.u))
            plot!(p, [integrator.t], [integrator.u[i]], seriestype = :scatter, label=false)
        end
    end
    # Display the plot
    display(p)
    f = plot(title="PosConstrNorm",
        xaxis="Time (t)", yaxis="u(t) (in μm)")
    plot!(f, t, PosConstrNorm, linewidth=3, label="Curve")
    display(f)
    f1 = plot(title="VelConstrNorm",
        xaxis="Time (t)", yaxis="u(t) (in μm)")
    plot!(f1, t, VelConstrNorm, linewidth=3, label="Curve")
    display(f1)


    return integrator
end

params=init_params(du₀, u₀)
itg=run_with_ode_solver_stepwise(params, t)