include("./slide_crank240_run.jl")
# Main entry point to control execution

# Main entry function to control execution and comparison of different runs
function main()
    results = ODERunResults([],[])
    t=10.5
    params1 = ODEParams(240, (0.0, t), Dict(:alg => Tsit5(), :reltol => 1e-6, :abstol => 1e-5,:dtmax=>1e-1, :progress => true))
    run(params1, results)
    # params2 = ODEParams(231, (0.0, t), Dict(:alg => Tsit5(), :reltol => 1e-5, :abstol => 1e-5,:dtmax=>1e-1, :progress => true))
    # run(params2, results)

    dfs=results.l_saved_data

    save(dfs, filename="jl240_solver_comparison.csv")

    # groups=[]
    # push!(groups,[Symbol("sol_1"), Symbol("sol_2"), Symbol("sol_3")])
    # draw3d(dfs, groups=groups, plot3d=true,limits = Dict("x" => [2, -2], "y" => [2, -2], "z" => [2, -2]))

    groups=[]

    push!(groups,[Symbol("sol_1"), Symbol("sol_2"), Symbol("sol_3")])
    #  for i in 1:15
    #     push!(groups,[Symbol("sol_$i")])
    #  end
    # push!(groups,[Symbol("sol_2")])
    # push!(groups,[Symbol("sol_10")])
    # push!(groups,[Symbol("sol_2"), Symbol("sol_10")])
    draw(dfs, groups=groups, plot3d=false)


end

main()