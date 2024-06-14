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
    MSRK5=Dict(:alg => OrdinaryDiffEq.MSRK5())
    MSRK6=Dict(:alg => OrdinaryDiffEq.MSRK6())
    Stepanov5 =Dict(:alg => OrdinaryDiffEq.Stepanov5()) #fix
    SIR54 =Dict(:alg => OrdinaryDiffEq.SIR54())
    Alshina2=Dict(:alg => OrdinaryDiffEq.Alshina2())
    Alshina3=Dict(:alg => OrdinaryDiffEq.Alshina3())
    Alshina6=Dict(:alg => OrdinaryDiffEq.Alshina6())
    #MSRK5","MSRK6",
    strs=["Stepanov5","SIR54","Alshina2","Alshina3","Alshina6"]

    tsetset=[
        # MSRK5 with default;
        # MSRK6 with default;
        Stepanov5 with default;
        SIR54 with default;
        Alshina2 with default;
        Alshina3 with default;
        Alshina6 with default;
    ]
    for (t,str) in zip(tsetset,strs)
        test(t,str)
    end
end

 main()   