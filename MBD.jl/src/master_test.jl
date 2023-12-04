# file_a.jl

# Define global variable

# Include file_b.jl
include("entry_DAE_II7.jl")

# Function that modifies the global variable and calls a function in file_b.jl
function modifyAndCall()
    values=[1 -1]
    for p1 in values
        for p2 in values
            for p3 in values
                for p4 in values
                    # Use a tuple (p1, p2, p3, p4) as the key
                    my_p = Dict("p1" => p1, "p2" => p2,"p3" => p3, "p4" => p4)
                    run(301, my_p)  # Call function from file_b.jl

                end
            end
        end
    end

end

# Example call
modifyAndCall()
