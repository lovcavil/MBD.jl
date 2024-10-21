include("./door309_run_20240310.jl")
using DifferentialEquations
# Main entry point to control execution
using Alert
# Main entry function to control execution and comparison of different runs
function test(testDict,str)
    results= ODERunResults([],[])
    time=2.2 #1.37#1.6#2.2

    params4 = ODEParams(341, (0.0, time), Dict(:alg => DP5(), :reltol => 1e-3, :abstol => 1e-3, :dtmin=>1e-6, :dtmax=>1e-3, :progress => true))
    #run(params4, results)
    #params5 = ODEParams(341, (0.0, time), Dict(:alg => BS3(), :reltol => 1e-3, :abstol => 1e-3, :dtmin=>1e-6, :dtmax=>1e-3, :progress => true))
    #params5 = ODEParams(341, (0.0, time), Dict(:alg => VCABM(), :reltol => 1e-3, :abstol => 1e-3, :dtmin=>1e-6, :dtmax=>1e-3, :progress => true))
    #params5 = ODEParams(341, (0.0, time), Dict(:alg => Vern7(), :reltol => 1e-3, :abstol => 1e-3, :dtmin=>1e-6, :dtmax=>1e-3, :progress => true))
    #params5 = ODEParams(341, (0.0, time), Dict(:alg => AB3(), :reltol => 1e-3, :abstol => 1e-3, :dt=>1e-3 ,:progress => true))
    # fixDict= Dict(:alg => AB3(), :dt=>1e-3 ,:progress => true)
    # testDict=Dict(:alg => VCAB3(), :reltol => 1e-4, :abstol => 1e-4, :dtmin=>1e-6, :dtmax=>1e-2, :progress => true)
    # refDict= Dict(:alg => VCABM(), :reltol => 1e-12, :abstol => 1e-12, :dtmin=>1e-15, :dtmax=>1e-3, :progress => true)
    params5 = ODEParams(341, (0.0, time),testDict)
    
    run(params5, results)
    
    dfs=results.l_saved_data
    save(dfs, filename="$(str).csv")
    
    alert("Your julia script is finished!")
end

function main()
    default=Dict(:reltol => 1e-4, :abstol => 1e-4, :dtmin=>1e-6, :dtmax=>1e-2, :progress => true)
    VCAB3=Dict(:alg => OrdinaryDiffEq.VCAB3())
    VCAB4=Dict(:alg => OrdinaryDiffEq.VCAB4())
    VCAB5=Dict(:alg => OrdinaryDiffEq.VCAB5())
    VCABM3=Dict(:alg => OrdinaryDiffEq.VCABM3())
    VCABM4=Dict(:alg => OrdinaryDiffEq.VCABM4())
    VCABM5=Dict(:alg => OrdinaryDiffEq.VCABM5())
    VCABM=Dict(:alg => OrdinaryDiffEq.VCABM())
    strs=["VCAB3";"VCAB4";"VCAB5";"VCABM3";"VCABM4";"VCABM5";"VCABM"]

    tsetset=[
        VCAB3 with default;
        VCAB4 with default;
        VCAB5 with default;
        VCABM3 with default;
        VCABM4 with default;
        VCABM5 with default;
        VCABM with default
    ]
    for (t,str) in zip(tsetset,strs)
        test(t,str)
    end
end

 main()   