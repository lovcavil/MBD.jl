using CSV
using DataFrames
using Dates

# Function to measure execution time and record it to a CSV file
function measure_and_save_time(measured_func::Function, params::ODEParams, 
    results::ODERunResults, runname::String,bf)
    # Measure the execution time of the provided function

    output_csv = "execution_times.csv"
    output_csv = joinpath(bf, "csv", output_csv)
    start_time = now()
    start_ns = time_ns()  # Get the current time in nanoseconds
    measured_func(params, results)
    end_ns = time_ns()  # Get the end time in nanoseconds
    elapsed_time = (end_ns - start_ns) / 1e9  # convert to seconds

    # Prepare a DataFrame to store the timing information
    timing_df = DataFrame(RunName=[runname],
                          Timestamp=[Dates.format(start_time, "yyyy-mm-dd HH:MM:SS")],
                          ExecutionTime=[elapsed_time])

    # Check if the output CSV file already exists
    if isfile(output_csv)
        # Append to the existing CSV file
        existing_df = CSV.read(output_csv, DataFrame)
        combined_df = vcat(existing_df, timing_df)
        CSV.write(output_csv, combined_df)
    else
        # Write a new CSV file
        CSV.write(output_csv, timing_df)
    end
end
