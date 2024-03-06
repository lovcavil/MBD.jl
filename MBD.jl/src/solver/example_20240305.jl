# Import necessary packages
include("solver.jl")
include("../problem/AppData_II7.jl")
using LinearAlgebra, DifferentialEquations, OrdinaryDiffEq, Sundials, Plots, CSV, DataFrames


# Function to initialize and solve an ODE problem for a specific application
function test_EI1(;app::Int=208, tspan::Tuple{Float64, Float64}=(0.0, 5.0))
    # Constants and initial conditions
    h0, hvar, g = 0.001, 0, 9.81

    # Load application data
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(app)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Log start
    println("START of Explicit Index 1")

    # Correction step
    t = tspan[1]
    qdd, Lam, ECond = mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀, du₀ = vcat(q, Lam, qd), vcat(qd, zeros(nc), qdd)

    # Problem setup
    p = [SMDT, STSDAT, SJDT, par]
    println("u₀=", u₀, "\ndu₀=", du₀)
    prob2 = ODEProblem(odequation, u₀, tspan, p)

    # Solve ODE
    solve_kwargs = Dict(:alg => Tsit5(), :reltol => 1e-6, :abstol => 1e-6, :progress => true)
    # sol = solve(prob2; solve_kwargs...)

    integrator = init(prob2; solve_kwargs...)
    step!(integrator)
    step!(integrator)
    step!(integrator)
    for (u, t) in tuples(integrator)
        @show u, t
    end
    # for i in integrator
    # end
    sol=0
    return sol
end

# Function to draw solution plots
function draw(sol)
    plot(sol, vars=1:3) |> display
    plot(sol, vars=8:11) |> display
    plot(sol, vars=13:15) |> display
end

# Function to save solution data
function save(sol; filename="jl_solver.csv")
    col_names = ["x1", "y1", "z1", "p1_1", "p1_2", "p1_3", "p1_4"]
    # Create a DataFrame using the column names
    df = DataFrame()
    df[!,"t"] = sol.t
    df[!,col_names[1]] = sol[1, :]
    df[!,col_names[2]] = sol[2, :]
    df[!,col_names[3]] = sol[3, :]
    CSV.write("jl_solver.csv", df)
    CSV.write(filename, df)
end

# Main execution block
function main()
    sol = test_EI1(app=205, tspan=(0.0,0.1))
#    draw(sol)
#    save(sol)
end

# using Cassette

# Cassette.@context PrintCtx
# const unique_function_names = Set{Symbol}()

# # Function to log unique function names
# function log_unique_function_name(f, args)
#     func_name = Symbol(f)
#     if !in(func_name, unique_function_names)
#         push!(unique_function_names, func_name)
#         println("Logging unique function name: ", func_name)
#     end
# end

# # Hook to intercept function calls
# Cassette.prehook(::PrintCtx, f, args...) = log_unique_function_name(f, args)
# Cassette.overdub(PrintCtx(), ()->main())

# # Display the unique function names logged
# println("Logged unique function names: ", unique_function_names)

main()