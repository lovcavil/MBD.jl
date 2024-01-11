include("solver.jl")
include("../AppData_II7.jl")
using LinearAlgebra
using DifferentialEquations
using OrdinaryDiffEq
using Sundials
using Plots

function e(out, du, u, p, t)
    out[1] = du[1] - u[1]
    out[2] = u[1]^2 + u[2]^2 - p
end

function test_II1()

    h0 = 0.001       # Initial time step
    hvar=0
    g = 9.81

    # Application Data Function
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(201)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Integration P
    t=0.0
    tspan = (0.0, 1.0)
    println("START")
    qdd, Lam, ECond=mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀  = vcat(q,Lam,qd)
    du₀ = vcat(qd,zeros(nc),qdd)

    differential_vars = create_bool_vector(nb, nc)
    simulation_params_dict = Dict(
        :g => 9.81,       
        :SMDT => SMDT,                   
        :STSDAT => STSDAT,               
        :SJDT => SJDT,                   
        :par => par                     
    )
    println("u₀=",u₀)
    println("du₀=",du₀)
    println("differential_vars=",differential_vars)  # 输出 differential_vars 的值

    prob = DAEProblem(implicitdae, du₀, u₀,tspan, simulation_params_dict, differential_vars=differential_vars)

    # Integration
    default_solve_kwargs = Dict(:alg=>IDA(), :reltol=>1e-6, :abstol=>1e-6, :progress=>true,:dtmin=>0.00001)
    sol = solve(prob; default_solve_kwargs...)

    return sol
end

function test_II3()

    h0 = 0.001       # Initial time step
    hvar=0
    g = 9.81

    # Application Data Function
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(201)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Integration P
    t=0.0
    tspan = (0.0, 1.0)
    println("START")
    qdd, Lam, ECond=mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀  = vcat(q,Lam,qd)
    du₀ = vcat(qd,zeros(nc),qdd)

    differential_vars = create_bool_vector(nb, nc)
    simulation_params_dict = Dict(
        :g => 9.81,       
        :SMDT => SMDT,                   
        :STSDAT => STSDAT,               
        :SJDT => SJDT,                   
        :par => par                     
    )
    println("u₀=",u₀)
    println("du₀=",du₀)
    println("differential_vars=",differential_vars)  # 输出 differential_vars 的值

    prob = DAEProblem(implicit3dae, du₀, u₀,tspan, simulation_params_dict, differential_vars=differential_vars)

    # Integration
    default_solve_kwargs = Dict(:alg=>IDA(), :reltol=>1e-6, :abstol=>1e-6, :progress=>true,:dtmin=>0.00001)
    sol = solve(prob; default_solve_kwargs...)

    return sol
end

function test_EI1()

    h0 = 0.001       # Initial time step
    hvar=0
    g = 9.81

    # Application Data Function
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(208)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Integration P
    t=0.0
    tspan = (0.0, 5)
    println("START")
    qdd, Lam, ECond=mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀  = vcat(q,Lam,qd)
    du₀ = vcat(qd,zeros(nc),qdd)

    differential_vars = create_bool_vector(nb, nc)

    p=Any[SMDT,STSDAT,SJDT,par]
    println("u₀=",u₀)
    println("du₀=",du₀)
    println("differential_vars=",differential_vars)  # 输出 differential_vars 的值

    prob2 = ODEProblem(odequation, u₀,tspan, p)

    # Integration
    default_solve_kwargs = Dict(:alg => Tsit5(), :reltol => 1e-10, :abstol => 1e-10, :progress => true)
    sol = solve(prob2;default_solve_kwargs...)
    return sol
end


#sol=test_II3()
# sol=test_II1()
# f01 = plot(sol)
# display(f01)

sol=test_EI1()

f02 = plot(sol, vars=1:3)
display(f02)
f03 = plot(sol, vars=8:11)
display(f03)
f04 = plot(sol, vars=13:15)
display(f04)