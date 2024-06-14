include("./door309_run_20240310.jl")
using DifferentialEquations
# Main entry point to control execution
using Alert
# Main entry function to control execution and comparison of different runs
function test(testDict,str)
    results= ODERunResults([],[])
    time=2.2 #1.37#1.6#2.2
    params5 = ODEParams(341, (0.0, time),testDict )
    
    run(params5, results)
    
    dfs=results.l_saved_data
    save(dfs, filename="$(str).csv")
    
    alert("Your julia script is finished!")
end

function main()
    default=Dict(:reltol => 1e-4, :abstol => 1e-4, :dtmin=>1e-6, :dtmax=>1e-2, :progress => true)
    #DP5=Dict(:alg => OrdinaryDiffEq.DP5())
    #Tsit5=Dict(:alg => OrdinaryDiffEq.Tsit5())
    #RKO65 =Dict(:alg => OrdinaryDiffEq.RKO65()) #fix
    #TanYam7 =Dict(:alg => OrdinaryDiffEq.TanYam7())
    #DP8=Dict(:alg => OrdinaryDiffEq.DP8())
    #TsitPap8=Dict(:alg => OrdinaryDiffEq.TsitPap8())
    #Feagin10=Dict(:alg => OrdinaryDiffEq.Feagin10())
    Feagin12=Dict(:alg => OrdinaryDiffEq.Feagin12())
    Feagin14=Dict(:alg => OrdinaryDiffEq.Feagin14())
    strs=["Feagin12","Feagin14"]

    tsetset=[
        #DP5 with default;
        #Tsit5 with default;
        #RKO65 with default;
        # TanYam7 with default;
        # DP8 with default;
        # TsitPap8 with default;
        # Feagin10 with default;
        Feagin12 with default;
        Feagin14 with default;
    ]
    for (t,str) in zip(tsetset,strs)
        test(t,str)
    end
end

 main()   