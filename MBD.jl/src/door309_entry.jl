include("./door309_run_20240310.jl")
# Main entry point to control execution

# Main entry function to control execution and comparison of different runs
function main()
    results = ODERunResults([])

    params1 = ODEParams(310, (0.0, 2.0), Dict(:alg => Tsit5(), :reltol => 1e-7, :abstol => 1e-7, :progress => true))
    run(params1, results)

    draw(results, groups=[1:3])
    draw(results, groups=[1:3,0*7+3:0*7+3,3*7+3:3*7+3,4*7+3:4*7+3],plot3d=false)
    draw(results, groups=[60+1:60+3,60+0*7+3:60+0*7+3,60+3*7+3:60+3*7+3,60+4*7+3:60+4*7+3],plot3d=false)
    #save(results, filename="jl_solver_comparison.csv", col_names=["x1", "y1", "z1", "p1_1", "p1_2", "p1_3", "p1_4"])
end

main()