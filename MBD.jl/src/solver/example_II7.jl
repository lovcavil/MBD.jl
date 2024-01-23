include("solver.jl")
include("../problem/AppData_II7.jl")
using LinearAlgebra
using DifferentialEquations
using OrdinaryDiffEq
using Sundials
using Plots
using CSV, DataFrames

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

function test_EI1(;app::Int=208,tspan::Tuple{Float64, Float64}=(0.0, 5.0))
    # CAKD initialization
    # & Application Data Function
    h0 = 0.001       # Initial time step
    hvar=0
    g = 9.81

    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(app)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Initialization
    println("START of Explicit Index 1")

    # Correction
    t=tspan[1]
    qdd, Lam, ECond=mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀  = vcat(q,Lam,qd)
    du₀ = vcat(qd,zeros(nc),qdd)

    p=Any[SMDT,STSDAT,SJDT,par]
    println("u₀=",u₀)
    println("du₀=",du₀)

    prob2 = ODEProblem(odequation, u₀,tspan, p)

    # Integration
    default_solve_kwargs = Dict(:alg => Tsit5(), :reltol => 1e-6, :abstol => 1e-6, :progress => true)
    sol = solve(prob2;default_solve_kwargs...)
    return sol
end

function draw(sol)
    f01 = plot(sol, vars=1:3)
    display(f01)
    f02 = plot(sol, vars=8:11)
    display(f02)
    f03 = plot(sol, vars=13:15)
    display(f03)
end

function save(sol)
    col_names = ["x1", "y1", "z1", "p1_1", "p1_2", "p1_3", "p1_4"]

    # Create a DataFrame using the column names
    df = DataFrame()
    df[!,"t"] = sol.t
    df[!,col_names[1]] = sol[1, :]
    df[!,col_names[2]] = sol[2, :]
    df[!,col_names[3]] = sol[3, :]

    CSV.write("jl_solver.csv", df)
end

sol=test_EI1(app=306,tspan=(0.0, 1.0))

draw(sol)
save(sol)