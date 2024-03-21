include("./door220_run.jl")
# Main entry point to control execution

# Main entry function to control execution and comparison of different runs
function main()
    results = ODERunResults([],[])

    params1 = ODEParams(220, (0.0, 10.0), Dict(:alg => Tsit5(), :reltol => 1e-10, :abstol => 1e-10, :progress => true))
    run(params1, results)
    params2 = ODEParams(221, (0.0, 10.0), Dict(:alg => Tsit5(), :reltol => 1e-10, :abstol => 1e-10, :progress => true))
    run(params2, results)
    params3 = ODEParams(222, (0.0, 10.0), Dict(:alg => Tsit5(), :reltol => 1e-10, :abstol => 1e-10, :progress => true))
    run(params3, results)
    params4 = ODEParams(223, (0.0, 10.0), Dict(:alg => Tsit5(), :reltol => 1e-10, :abstol => 1e-10, :progress => true))
    run(params4, results)
    params5 = ODEParams(224, (0.0, 10.0), Dict(:alg => Tsit5(), :reltol => 1e-10, :abstol => 1e-10, :progress => true))
    run(params5, results)
    params5 = ODEParams(225, (0.0, 10.0), Dict(:alg => Tsit5(), :reltol => 1e-10, :abstol => 1e-10, :progress => true))
    run(params5, results)
    dfs=results.l_saved_data

    groups=[]
    push!(groups,[Symbol("sol_1"), Symbol("sol_2"), Symbol("sol_3")])
    draw3d(dfs, groups=groups, plot3d=true,limits = Dict("x" => [2, -2], "y" => [2, -2], "z" => [2, -2]))

    groups=[]
    # push!(groups,[Symbol("sol_1")])
    # push!(groups,[Symbol("sol_2")])
    # push!(groups,[Symbol("sol_3")])
    push!(groups,[Symbol("sol_1"), Symbol("sol_2"), Symbol("sol_3")])
     for i in 1:20
        push!(groups,[Symbol("sol_$i")])
     end
    draw(dfs, groups=groups, plot3d=false)

    save(dfs, filename="jl220_solver_comparison.csv")
end

main()