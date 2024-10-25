
include("solver_utils.jl")

using OrdinaryDiffEq
using LinearSolve
function main()
    # Solver Group (Set 1)

    solver_dict = Dict(
        "Rodas5"=>"Rodas5",
        "Rosenbrock23"=>"Rosenbrock23", 
        #"TRBDF2"=>TRBDF2,
        "Rodas5P"=>"Rodas5P",
    )    
    solver_dict = Dict(

        "Trapezoid"=>"Trapezoid",
        "Kvaerno3" => "Kvaerno3",  # A-L stable stiffly-accurate 3rd order ESDIRK method
        "KenCarp3" => "KenCarp3",  # A-L stable stiffly-accurate 3rd order ESDIRK method with splitting
        "KenCarp47" => "KenCarp47",  # A-L stable stiffly-accurate 4th order seven-stage ESDIRK method with splitting
        "ESDIRK436L2SA2" => "ESDIRK436L2SA2",  # A-L stable stiffly-accurate 4th order six-stage ESDIRK method
        "ESDIRK437L2SA" => "ESDIRK437L2SA",  # A-L stable stiffly-accurate 4th order seven-stage ESDIRK method
        "KenCarp5" => "KenCarp5",  # A-L stable stiffly-accurate 5th order ESDIRK method with splitting
        "ESDIRK54I8L2SA" => "ESDIRK54I8L2SA",  # A-L stable stiffly-accurate 5th order eight-stage ESDIRK method
        "RosShamp4" => "RosShamp4",  # 4th order A-stable Rosenbrock method
        "Veldd4" => "Veldd4",  # 4th order D-stable Rosenbrock method
        "Velds4" => "Velds4",  # 4th order A-stable Rosenbrock method
        "GRK4T" => "GRK4T",  # 4th order efficient Rosenbrock method
        "GRK4A" => "GRK4A",  # 4th order A-stable Rosenbrock method
        "Ros4LStab" => "Ros4LStab",  # 4th order L-stable Rosenbrock method
        "Rodas4" => "Rodas4",  # 4th order A-stable stiffly stable Rosenbrock method with stiff-aware 3rd order interpolant
        "Rodas42" => "Rodas42",  # 4th order A-stable stiffly stable Rosenbrock method with stiff-aware 3rd order interpolant
        "Rodas4P" => "Rodas4P",  # 4th order A-stable stiffly stable Rosenbrock method, accurate on linear parabolic and nonlinear problems
        "Rodas4P2" => "Rodas4P2",  # 4th order L-stable stiffly stable Rosenbrock method, improved for inexact Jacobians
        "Rodas5" => "Rodas5",  # 5th order A-stable stiffly stable Rosenbrock method with stiff-aware 4th order interpolant
        "ROS3" => "ROS3",  # 3rd order L-stable Rosenbrock-Wanner method with 3 internal stages
        "ROS2PR" => "ROS2PR",  # 2nd order stiffly accurate Rosenbrock-Wanner method with 3 internal stages
        "Scholz4_7" => "Scholz4_7",  # 3rd order stiffly accurate Rosenbrock-Wanner method (converges with order 4 for stiff cases)
        "ROS3PRL" => "ROS3PRL",  # 3rd order stiffly accurate Rosenbrock-Wanner method with 4 internal stages
         "FBDF" => "FBDF",    # Fixed-leading coefficient adaptive-order BDF method
    )
    
    autodiff_dict = Dict(
        #"autodiff_true" => (true,Val{:forward}),
        "autodiff_false_forward" => [false, Val{:forward}],
        #"autodiff_false_central" => (false, Val{:central}),
        #"autodiff_false_complex" => (false, Val{:complex}),
    )
    linsolve_dict = Dict(
        ""=> nothing,
        #"QR" => QRFactorization(),
        #"KrylovJLCG" => KrylovJL_CG(),
        "KrylovJLGMRES" => KrylovJL_GMRES(),
        "CholeskyFactorization" => CholeskyFactorization(),
        #"BunchKaufmanFactorization" => BunchKaufmanFactorization(),
        #"CHOLMODFactorization" => CHOLMODFactorization(),
        "NormalCholeskyFactorization" => NormalCholeskyFactorization(),
        "NormalBunchKaufmanFactorization" => NormalBunchKaufmanFactorization(),
    )
    linsolve_dict = Dict(
        ""=> nothing,
        "LUFactorization" => LUFactorization(),
        "GenericLUFactorization" => GenericLUFactorization(),
        "QRFactorization" => QRFactorization(),
        "SVDFactorization" => SVDFactorization(),
        "SimpleLUFactorization" => SimpleLUFactorization(),
        #"KLUFactorization" => KLUFactorization(),
        #"UMFPACKFactorization" => UMFPACKFactorization(),
        #"KrylovJL_MINRES" => KrylovJL_MINRES(),
        "KrylovJL_GMRES" => KrylovJL_GMRES(),
        "KrylovJL_BICGSTAB" => KrylovJL_BICGSTAB(),
        #"KrylovJL_LSMR" => KrylovJL_LSMR(),
        #"KrylovJL_CRAIGMR" =>KrylovJL_CRAIGMR(),
    )
    linsolve_dict = Dict(
        ""=> nothing,)
    # Tolerance Group (Set 2)
    tolerance_dict = Dict(
        #"tol_1e8" => (1e-8, 1e-8),
        "tol_1e3" => [1e-3, 1e-3],
        "tol_1e6" => (1e-6, 1e-6),
        #"tol_1e9" => (1e-9, 1e-9)
    )
    # Timespan Group (Set 3)
    timespan_dict = Dict(
        "short_timespan" => [0.0, 1.0],
        "short_timespan2" => [0.0, 10.0],
        #"long_timespan" => (0.0, 1000.0)
    )
    initial_state_dict = Dict(
        "state1" => [0.0, 0.1, 0.5, 0.0])
    timeout_dict=Dict(
        "tLimit"=>60.0
    )    
    # Combine the configurations
    param_combinations = generate_param_combinations([solver_dict,autodiff_dict,linsolve_dict, tolerance_dict, timespan_dict,initial_state_dict,timeout_dict])
    # Convert to a list for further use
    println(param_combinations)
    json_compatible_dict = convert_to_json_compatible(param_combinations)
    # Save to a JSON file
    file_path = joinpath(@__DIR__, "my_data.json")
    open(file_path, "w") do file
        JSON.print(file, json_compatible_dict)
    end

    println("End")
end
# Execute the main function
println(Threads.nthreads())
main()
