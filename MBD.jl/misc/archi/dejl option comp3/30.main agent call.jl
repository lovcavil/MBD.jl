
include("solver_utils.jl")

using OrdinaryDiffEq
using LinearSolve
using JSON

# Main function
function main(folder,num,json_string)
    # Parse JSON string into a Julia dictionary
    my_dict_back = JSON.parse(json_string)
    
    # Convert JSON-compatible data back to Julia types
    my_dict_back = convert_from_json_compatible(my_dict_back)

    # Convert to solver configurations
    final_configs = convert_to_solver_config(my_dict_back)
    solver_configs = collect(values(final_configs))
    fname="$(string(num)).csv"
    # Run solvers and store solutions
    sols, energies, times, ΔEs = run_solvers(solver_configs, folder,fname, 5)
    
        # Write the final result to a JSON file
    fname="$(string(num)).json"
    output_json_file = joinpath(@__DIR__,folder, fname)
    open(output_json_file, "w") do file
        JSON.print(file, Dict("energies" => energies, "times" => times))
    end
    return 0

    # Print the results to be captured by Python
    #println(JSON.json(Dict("energies" => energies, "times" => times)))
end

# Get the JSON string from command-line arguments
if length(ARGS) > 0
    folder=ARGS[1]
    num = ARGS[2]
    json_string = ARGS[3]
    println(ARGS)
    main(folder,num,json_string)
else
    println("Please provide the JSON string as a command-line argument.")
end