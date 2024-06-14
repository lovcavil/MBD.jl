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
    Midpoint=Dict(:alg => OrdinaryDiffEq.Midpoint())
    Heun=Dict(:alg => OrdinaryDiffEq.Heun())
    Ralston =Dict(:alg => OrdinaryDiffEq.Ralston())
    RK4 =Dict(:alg => OrdinaryDiffEq.RK4())
    BS3=Dict(:alg => OrdinaryDiffEq.BS3())
    OwrenZen3=Dict(:alg => OrdinaryDiffEq.OwrenZen3())
    OwrenZen4=Dict(:alg => OrdinaryDiffEq.OwrenZen4())
    OwrenZen5=Dict(:alg => OrdinaryDiffEq.OwrenZen5())
    strs=["Midpoint","Heun","Ralston","RK4","BS3","OwrenZen3","OwrenZen4","OwrenZen5"]

    tsetset=[
        Midpoint with default;
        Heun with default;
        Ralston with default;
        RK4 with default;
        BS3 with default;
        OwrenZen3 with default;
        OwrenZen4 with default;
        OwrenZen5 with default;
    ]
    for (t,str) in zip(tsetset,strs)
        test(t,str)
    end
end

 main()   