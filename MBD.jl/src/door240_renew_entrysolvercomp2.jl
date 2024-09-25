# using diffeq.jl with correct

include("./door240_renew_run_202403102.jl")
using DifferentialEquations
# Main entry point to control execution
using Alert
# Main entry function to control execution and comparison of different runs

# problem 
function test(sol_Dict,par_Dict)
    #sol_name,config_cfg,contact_cfg,app
    results= ODERunResults([],[])
    params = ODEParams(par_Dict[:app],par_Dict[:contact_cfg], (0.0, par_Dict[:time]),sol_Dict )
    run_name="$(par_Dict[:runname])_$(par_Dict[:contact_cfg])_$(par_Dict[:sol_name])"
    bf=par_Dict[:bf]
    measure_and_save_time(run, params, results, run_name,bf)
    dfs=results.l_saved_data
    filename="$(run_name).csv"
    save(dfs, bf, filename)
    alert("Your julia script is finished!")
end

function main()
    # 
    if length(ARGS)>1
        dt = parse(Float64, ARGS[1])
        config_cfg=ARGS[2]
        contact_cfg=ARGS[3]
        app=parse(Int, ARGS[4])
        time=parse(Float64, ARGS[5])
        runname=ARGS[6]
        bf=ARGS[7]
        reltol = parse(Float64, ARGS[8])
        abstol = parse(Float64, ARGS[9])
        dtmax = parse(Float64, ARGS[10])
        dtmin = parse(Float64, ARGS[11])
        sol_name=ARGS[12]
    else
        dt = 0.0001
        config_cfg=""
        contact_cfg="contact_g65"
        app=240
        time=3.5
        time=1.0
        runname="test_240"
        bf="D:/OneDrive/Articles/10.Working/[D21][20211009]ContactMechanics/MBD.jl/"
        reltol = 1e-6
        abstol = 1e-5
        dtmin = 1e-10
        dtmax = 1e-3
        sol_name="Tsit5"
    end
    if sol_name=="Euler"
        sol_Dict=getEulerDict(dt)
    end
    if sol_name=="Tsit5"
        sol_Dict=getTsit5Dict(reltol,abstol,dtmin,dtmax)
    end
    par_Dict = Dict(:sol_name=>sol_name,:config_cfg=>config_cfg,:contact_cfg=>contact_cfg,:app=>app
    ,:time=>time,:runname=>runname,:bf=>bf,
    :reltol=>reltol,:abstol=>abstol,:dtmin=>dtmin,:dtmax=>dtmax)
    test(sol_Dict,par_Dict)

end
function getEulerDict(dt)
    fixDict= Dict(:dt=>dt ,:progress => true)
    Euler = Dict( :alg =>  OrdinaryDiffEq.Euler() ) 
    sol_Dict = merge(Euler, fixDict)
    return sol_Dict
end

function getTsit5Dict(reltol,abstol,dtmin,dtmax)
    default=Dict(:reltol => reltol, :abstol => abstol, :dtmin=>dtmin, :dtmax=>dtmax, :progress => true)
    Euler = Dict( :alg =>  OrdinaryDiffEq.Tsit5() ) 
    sol_Dict = merge(Euler, default)
    return sol_Dict
end

main()
main()