include("./door309_run_20240310.jl")
# Main entry point to control execution

# Main entry function to control execution and comparison of different runs
function main()
    results = ODERunResults([],[])

    # params1 = ODEParams(310, (0.0, 0.5), Dict(:alg => Tsit5(), :reltol => 1e-5, :abstol => 1e-5, :progress => true))
    # run(params1, results)
    # params2 = ODEParams(321, (0.0, 0.5), Dict(:alg => Tsit5(), :reltol => 1e-5, :abstol => 1e-5, :progress => true))
    # run(params2, results)
    params3 = ODEParams(330, (0.0, 0.85), Dict(:alg => Tsit5(), :reltol => 1e-5, :abstol => 1e-5, :progress => true))
    run(params3, results)
    # params4 = ODEParams(322, (0.0, 0.88), Dict(:alg => Tsit5(), :reltol => 1e-5, :abstol => 1e-5, :progress => true))
    # run(params4, results)
    dfs=results.l_saved_data

    groups=[]
    push!(groups,[Symbol("sol_1"), Symbol("sol_2"), Symbol("sol_3")])
    push!(groups,[Symbol("sol_$i") for i in 3*7+1:3*7+3])
    push!(groups,[Symbol("sol_$i") for i in 4*7+1:4*7+3])
    push!(groups,[Symbol("sol_$i") for i in 5*7+1:5*7+3])
    push!(groups,[Symbol("sol_$i") for i in 6*7+1:6*7+3])
    push!(groups,[Symbol("sol_$i") for i in 7*7+1:7*7+3])
    push!(groups,[Symbol("sol_$i") for i in 8*7+1:8*7+3])
    push!(groups,[Symbol("sol_$i") for i in 9*7+1:9*7+3])
    draw3d(dfs, groups=groups, plot3d=true)

    groups=[]
    push!(groups,[Symbol("sol_1")])
    push!(groups,[Symbol("sol_2")])
    push!(groups,[Symbol("sol_3")])
    push!(groups,[Symbol("sol_8")])
    push!(groups,[Symbol("sol_9")])
    push!(groups,[Symbol("sol_10")])
    push!(groups,[Symbol("sol_15")])
    push!(groups,[Symbol("sol_16")])
    push!(groups,[Symbol("sol_17")])
    push!(groups,[Symbol("sol_$i") for i in 1*7+1:1*7+3])
    push!(groups,[Symbol("sol_$i") for i in 2*7+1:2*7+3])
    push!(groups,[Symbol("sol_$i") for i in 3*7+1:3*7+3])
    push!(groups,[Symbol("sol_$i") for i in 4*7+1:4*7+3])
    push!(groups,[Symbol("sol_$i") for i in 5*7+1:5*7+3])
    push!(groups,[Symbol("sol_$i") for i in 6*7+1:6*7+3])
    push!(groups,[Symbol("sol_$i") for i in 7*7+1:7*7+3])
    push!(groups,[Symbol("sol_$i") for i in 8*7+1:8*7+3])
    push!(groups,[Symbol("sol_$i") for i in 9*7+1:9*7+3])
    # push!(groups,[Symbol("sol_$i") for i in 0*7+4:0*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 1*7+4:1*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 2*7+4:2*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 3*7+4:3*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 4*7+4:4*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 5*7+4:5*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 6*7+4:6*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 7*7+4:7*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 8*7+4:8*7+7])
    # push!(groups,[Symbol("sol_$i") for i in 9*7+4:9*7+7])
    push!(groups,[Symbol("sol_$i") for i in [3*7+3,5*7+3,6*7+3]])
    push!(groups,[Symbol("sol_$i") for i in [3*7+3]])
    push!(groups,[Symbol("sol_$i") for i in [5*7+3]])
    push!(groups,[Symbol("sol_$i") for i in [6*7+3]])
    sec2=70
    sec3=139
    s=1
    e=3
    push!(groups,[Symbol("sol_$i") for i in sec3+s:sec3+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+1*7+s:sec3+1*7+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+2*7+s:sec3+2*7+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+3*7+s:sec3+3*7+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+4*7+s:sec3+4*7+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+5*7+s:sec3+5*7+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+6*7+s:sec3+6*7+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+7*7+s:sec3+7*7+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+8*7+s:sec3+8*7+e])
    # push!(groups,[Symbol("sol_$i") for i in sec3+9*7+s:sec3+9*7+e])
    push!(groups,[Symbol("ex_$i") for i in 3:4])
    push!(groups,[Symbol("ex_$i") for i in 5:6])
    push!(groups,[Symbol("ex_$i") for i in 7:8])
    push!(groups,[Symbol("ex_$i") for i in 9:10])
    # push!(groups,[Symbol("sol_$i") for i in 71:80])
    #  for i in 1:length(dfs[1][:,1])
    #     push!(groups,[Symbol("sol_$i")])
    #  end
    draw(dfs, groups=groups, plot3d=false)

    # draw(results, groups=[1:3,1*7+1:1*7+3,2*7+1:2*7+3])
    # draw(results, groups=[1:3,1*7+1:1*7+3,2*7+1:2*7+3,3*7+3:3*7+3,4*7+3:4*7+3],plot3d=false)
    #draw(results, groups=[1:3,0*7+3:0*7+3,3*7+3:3*7+3,4*7+3:4*7+3],plot3d=false)
    # draw(results, groups=[65+1:65+3,65+0*7+3:65+0*7+3,65+3*7+3:65+3*7+3,65+4*7+3:65+4*7+3],plot3d=false)
    #save(results, filename="jl_solver_comparison.csv", col_names=["x1", "y1", "z1", "p1_1", "p1_2", "p1_3", "p1_4"])
    save(dfs, filename="jl_solver_comparison.csv")
end

main()