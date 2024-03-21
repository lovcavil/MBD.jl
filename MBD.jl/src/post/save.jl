# Modify the save function to handle multiple solutions
function save(dfs; filename="jl_solver.csv", col_names=String[])
    for (index, df) in enumerate(dfs)
        # df = DataFrame()
        # df[!,"t"] = sol.t
        # for (i, col_name) in enumerate(col_names)
        #     df[!,col_name] = sol[i, :]
        # end
        specific_filename = replace(filename, ".csv" => "_run$(index).csv")
        CSV.write(specific_filename, df)
    end
end