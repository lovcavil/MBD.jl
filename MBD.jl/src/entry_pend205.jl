include("./run_pend205_20240305.jl")
# Main entry point to control execution

# Main entry function to control execution and comparison of different runs
function main()
    results = ODERunResults([])

    params1 = ODEParams(205, (0.0, 6.0), Dict(:alg => Tsit5(), :reltol => 1e-7, :abstol => 1e-7, :progress => true))
    run(params1, results)

    params2 = ODEParams(209, (0.0, 6.0), Dict(:alg => Tsit5(), :reltol => 1e-7, :abstol => 1e-7, :progress => true))
    run(params2, results)

    params3 = ODEParams(210, (0.0, 6.0), Dict(:alg => Tsit5(), :reltol => 1e-7, :abstol => 1e-7, :progress => true))
    run(params3, results)

    draw(results, groups=[1:3])
    draw(results, groups=[1:1,2:2,3:3],plot3d=false)
    save(results, filename="jl_solver_comparison.csv", col_names=["x1", "y1", "z1", "p1_1", "p1_2", "p1_3", "p1_4"])
end

main()