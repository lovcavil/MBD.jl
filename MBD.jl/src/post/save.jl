using CSV
using Dates

function save(dfs,  bf, filename="jl_solver.csv")
    # Get the current timestamp
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    corrected_path = joinpath(bf, "csv", filename)

    base_name = replace(basename(corrected_path), ".csv" => "")
    dir_name = dirname(corrected_path)


    for (index, df) in enumerate(dfs)
        # Construct a unique filename for each DataFrame
        base_name = replace(basename(corrected_path), ".csv" => "")
        dir_name = dirname(corrected_path)
    
        specific_filename = joinpath(dir_name, "$(timestamp)_$(base_name)_run$(index).csv")
        specific_filename2 = joinpath(dir_name, "$(base_name)_run$(index).csv")
        
        # Write the DataFrame to a CSV file
        CSV.write(specific_filename, df)
        CSV.write(specific_filename2, df)
    end
end