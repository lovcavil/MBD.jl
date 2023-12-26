include("mathfunction_II7.jl")
include("AppData_II7.jl")
include("ImplicitIndex1.jl")
include("ImplicitIndex3.jl")
include("ExplicitNystrom4.jl")
include("ExplicitRKFN45.jl")
include("solver/ExplicitIndex1Form.jl")
using LinearAlgebra
using CSV, DataFrames
using IterativeRefinement
using Plots

function run(app,my_p)

    h0 = 0.002       # Initial time step
    h = h0
    hmax = 0.01
    hvar = 1          # hvar=1, variable h; hvar=2, constant h

    tfinal = 1

    g = 9.81

    # Application Data Function
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(app)
    
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA,my_p]

    # Integration P
    t=0.0
    println("START")
    qdd, Lam, ECond=mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    #u₀  = vcat(q,0,0,0,0,qd)
    #du₀ = vcat(qd,zeros(nc),0,0,-9.81,zeros(4))
    u₀  = vcat(q,Lam,qd)
    du₀ = vcat(qd,zeros(nc),qdd)
    out = zeros(14*nb+nc)
    differential_vars = vcat(repeat([true], 7*nb), repeat([false], nc), repeat([true], 7*nb))
    println("u₀=",u₀)
    println("du₀=",du₀)
    #params=init_params(du₀, u₀,differential_vars,function_pendulum,SMDT, STSDAT, SJDT, par)
    params=init_params(du₀, u₀,SMDT, STSDAT, SJDT, par)
    # Integration

    sol=sol_O(out, params, t)
    
    Plots.default(show = true)
    f0 = plot(sol, xlabel="Time", ylabel="Value", title="sss")
    display(f0)
    col_names = ["x1", "y1", "z1", "p1_1", "p1_2", "p1_3", "p1_4"]

    # Create a DataFrame using the column names
    df = DataFrame()
    df[!,"t"] = sol.t
    df[!,col_names[1]] = sol[1, :]
    df[!,col_names[2]] = sol[2, :]
    df[!,col_names[3]] = sol[3, :]

    CSV.write("jl_solver.csv", df)
end

#run(101)

my_p = Dict("p1" => 1, "p2" => 1,"p3" => 1, "p4" => 1)

run(205, my_p)  # Call function from file_b.jl