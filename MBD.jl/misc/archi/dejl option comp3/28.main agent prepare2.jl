
include("solver_utils.jl")

#using OrdinaryDiffEq
using LinearSolve
function main()
    solver_dict = Dict(
        "ImplicitEuler"=>"ImplicitEuler",
        #"ImplicitMidpoint"=>"ImplicitMidpoint", 
        #"Trapezoid"=>"Trapezoid",
        "TRBDF2"=>"TRBDF2",
        #"SDIRK2" => "SDIRK2",  # A-B-L stable 2nd order SDIRK method

        # Third-order ESDIRK methods
        # "Kvaerno3" => "Kvaerno3",  # A-L stable stiffly-accurate 3rd order ESDIRK method
        # "KenCarp3" => "KenCarp3",  # A-L stable stiffly-accurate 3rd order ESDIRK method with splitting
    
        # Fourth-order SDIRK methods
        # "Cash4" => "Cash4",  # A-L stable 4th order SDIRK method
        # "Hairer4" => "Hairer4",  # A-L stable 4th order SDIRK method
        # "Hairer42" => "Hairer42",  # A-L stable 4th order SDIRK method
    
        # Fourth-order ESDIRK methods
        # "Kvaerno4" => "Kvaerno4",  # A-L stable stiffly-accurate 4th order ESDIRK method
        # "KenCarp4" => "KenCarp4",  # A-L stable stiffly-accurate 4th order ESDIRK method with splitting
        # "KenCarp47" => "KenCarp47",  # A-L stable stiffly-accurate 4th order seven-stage ESDIRK method with splitting
        "ESDIRK436L2SA2" => "ESDIRK436L2SA2",  # A-L stable stiffly-accurate 4th order six-stage ESDIRK method
        # "ESDIRK437L2SA" => "ESDIRK437L2SA",  # A-L stable stiffly-accurate 4th order seven-stage ESDIRK method
    
        # Fifth-order ESDIRK methods
        # "Kvaerno5" => "Kvaerno5",  # A-L stable stiffly-accurate 5th order ESDIRK method
        "KenCarp5" => "KenCarp5",  # A-L stable stiffly-accurate 5th order ESDIRK method with splitting
        "KenCarp58" => "KenCarp58",  # A-L stable stiffly-accurate 5th order eight-stage ESDIRK method with splitting
        # "ESDIRK54I8L2SA" => "ESDIRK54I8L2SA",  # A-L stable stiffly-accurate 5th order eight-stage ESDIRK method
        "ESDIRK547L2SA2" => "ESDIRK547L2SA2",  # A-L stable stiffly-accurate 5th order seven-stage ESDIRK method

        # Fully-Implicit Runge-Kutta Methods (FIRK)
        # "RadauIIA3" => "RadauIIA3",  # A-B-L stable fully implicit Runge-Kutta method (3rd order)
        "RadauIIA5" => "RadauIIA5",  # A-B-L stable fully implicit Runge-Kutta method (5th order)
    
        # Parallel Diagonally Implicit Runge-Kutta Methods
        # "PDIRK44" => "PDIRK44",  # 2 processor 4th order diagonally non-adaptive implicit method
    
        # Rosenbrock Methods
        # "ROS3P" => "ROS3P",  # 3rd order A-stable and stiffly stable Rosenbrock method
        # "Rodas3" => "Rodas3",  # 3rd order A-stable and stiffly stable Rosenbrock method
        # "Rodas3P" => "Rodas3P",  # 3rd order A-stable stiffly stable Rosenbrock with stiff-aware 3rd order interpolant
        # "RosShamp4" => "RosShamp4",  # 4th order A-stable Rosenbrock method
        # "Veldd4" => "Veldd4",  # 4th order D-stable Rosenbrock method
        # "Velds4" => "Velds4",  # 4th order A-stable Rosenbrock method
        # "GRK4T" => "GRK4T",  # 4th order efficient Rosenbrock method
        # "GRK4A" => "GRK4A",  # 4th order A-stable Rosenbrock method
        # "Ros4LStab" => "Ros4LStab",  # 4th order L-stable Rosenbrock method
        "Rodas4" => "Rodas4",  # 4th order A-stable stiffly stable Rosenbrock method with stiff-aware 3rd order interpolant
        "Rodas42" => "Rodas42",  # 4th order A-stable stiffly stable Rosenbrock method with stiff-aware 3rd order interpolant
        "Rodas4P" => "Rodas4P",  # 4th order A-stable stiffly stable Rosenbrock method, accurate on linear parabolic and nonlinear problems
        "Rodas4P2" => "Rodas4P2",  # 4th order L-stable stiffly stable Rosenbrock method, improved for inexact Jacobians
        "Rodas5" => "Rodas5",  # 5th order A-stable stiffly stable Rosenbrock method with stiff-aware 4th order interpolant
        "Rodas5P" => "Rodas5P",  # 5th order A-stable stiffly stable Rosenbrock method with improved stability for adaptive time stepping
        # "ROS2" => "ROS2",  # 2nd order L-stable Rosenbrock-Wanner method with 2 internal stages
        # "ROS3" => "ROS3",  # 3rd order L-stable Rosenbrock-Wanner method with 3 internal stages
        # "ROS2PR" => "ROS2PR",  # 2nd order stiffly accurate Rosenbrock-Wanner method with 3 internal stages
        # "ROS3PR" => "ROS3PR",  # 3rd order stiffly accurate Rosenbrock-Wanner method with 3 internal stages
        # "Scholz4_7" => "Scholz4_7",  # 3rd order stiffly accurate Rosenbrock-Wanner method (converges with order 4 for stiff cases)
        #"ROS3PRL" => "ROS3PRL",  # 3rd order stiffly accurate Rosenbrock-Wanner method with 4 internal stages
        # "ROS3PRL2" => "ROS3PRL2",  # 3rd order stiffly accurate Rosenbrock-Wanner method with 4 internal stages (improved version)
        # Rosenbrock-W Methods
        # NDF (Numerical Differentiation Function) Methods
        "QNDF1" => "QNDF1",  # Adaptive order 1 L-stable quasi-constant timestep NDF method
        "QNDF2" => "QNDF2",  # Adaptive order 2 L-stable quasi-constant timestep NDF method
        "QNDF" => "QNDF",    # Adaptive order quasi-constant timestep NDF method (similar to ode15s)
        # BDF (Backward Differentiation Formula) Methods
        "QBDF1" => "QBDF1",  # Adaptive order 1 L-stable BDF method (similar to implicit Euler)
        "QBDF2" => "QBDF2",  # Adaptive order 2 L-stable BDF method with quasi-constant timesteps
        "QBDF" => "QBDF",    # Adaptive order quasi-constant timestep BDF method
        "ABDF2" => "ABDF2",  # Adaptive order 2 L-stable fixed leading coefficient multistep BDF method
        "FBDF" => "FBDF",    # Fixed-leading coefficient adaptive-order BDF method
        # Modified BDF Methods
        # "MEBDF2" => "MEBDF2"  # Second order Modified Extended BDF method (fixed timestep)
    )
    
    autodiff_dict = Dict(
        #"autodiff_true" => [true,Val{:forward}],
         "autodiff_false_forward" => [false, Val{:forward}],
        #  "autodiff_false_central" => [false, Val{:central}],
        #  "autodiff_false_complex" => [false, Val{:complex}],
    )

    linsolve_dict = Dict(
        # ""=> nothing,
        # "QRFactorization" => "QRFactorization",
        "LUFactorization" => "LUFactorization",
        "GenericLUFactorization" => "GenericLUFactorization",
        "SVDFactorization" => "SVDFactorization",
        # "SimpleLUFactorization" => "SimpleLUFactorization",
        # "KrylovJL_GMRES" => "KrylovJL_GMRES",
        # "KrylovJL_BICGSTAB" => "KrylovJL_BICGSTAB",
    )
    # Tolerance Group (Set 2)
    tolerance_dict = Dict(
        #"tol_1e12" => [1e-12, 1e-12],
        "tol_1e3" => [1e-3, 1e-3],
        "tol_1e6" => [1e-6, 1e-6],
        #"tol_1e9" => [1e-9, 1e-9]
    )
    tolerance_dict = Dict(
        "tol_1e$i" => [10.0^(-i), 10.0^(-i)] for i in 3:12
    )
    # Timespan Group (Set 3)
    timespan_dict = Dict(
        #"short_timespan" => [0.0, 1.0],
        "short_timespan2" => [0.0, 3.0],
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
