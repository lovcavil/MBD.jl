include("./door309_run_20240310.jl")
# Main entry point to control execution
using Alert
# Main entry function to control execution and comparison of different runs

function test(sol_Dict,sol_name,config_cfg,contact_cfg)
    results= ODERunResults([],[])
    time=2.2
    params = ODEParams(350,"$contact_cfg.json", (0.0, time), Dict(:alg => Tsit5(),
     :reltol => 1e-2, :abstol => 1e-2, :dtmin=>1e-5, :dtmax=>1e-1, :progress => true))
    params = ODEParams(341,contact_cfg, (0.0, time),sol_Dict )
    run_name="$(contact_cfg)_$(sol_name)"
    measure_and_save_time(run, params, results)
    dfs=results.l_saved_data
    save(dfs, filename="$(run_name).csv")
    alert("Your julia script is finished!")
end

function main()
    # default=Dict(:reltol => 1e-4, :abstol => 1e-4, :dtmin=>1e-6, :dtmax=>1e-2, :progress => true)
    #println(ARGS)
    dt = parse(Float64, ARGS[1])
    config_cfg=ARGS[2]
    contact_cfg=ARGS[3]
    fixDict= Dict(:dt=>dt ,:progress => true)
    Euler = Dict( :alg =>  OrdinaryDiffEq.Euler() ) 
    strs=["Euler"]
    tests=[
        merge(Euler, fixDict);
    ]
    for (sol_Dict,sol_name) in zip(tests,strs)
        test(sol_Dict,sol_name,config_cfg,contact_cfg)
    end
end



main()