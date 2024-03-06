include("./run_pend205_20240305.jl")
# Main entry point to control execution

# Main entry function to control execution and comparison of different runs
function main()
    results = ODERunResults([])

    params1 = ODEParams(205, (0.0, 6.0), Dict(:alg => Tsit5(), :reltol => 1e-6, :abstol => 1e-6, :progress => true))
    run(params1, results)

    params2 = ODEParams(205, (0.0, 3.0), Dict(:alg => Tsit5(), :reltol => 1e-3, :abstol => 1e-3, :progress => true))
    run(params2, results)

    draw(results, groups=[1:3, 8:11, 13:15])
    save(results, filename="jl_solver_comparison.csv", col_names=["x1", "y1", "z1", "p1_1", "p1_2", "p1_3", "p1_4"])
end

main()