function process_vector(q::Vector, flags::Vector{Bool}, target::Vector{Vector{Int}}, func::Function)
    # Ensure `q` and `flags` have the same length
    if length(q) != length(flags)
        error("Vectors 'q' and 'flags' must be of the same length.")
    end

    # Process elements of `q` based on `flags`
    for (i, flag) in enumerate(flags)
        if flag && i <= length(target)
            func(target[i], q[i])
        end
    end
end

# Function to append `element` of `q` to the `target_vector`
function append_to_vector(target_vector::Vector{Int}, element)
    push!(target_vector, element)
end

# Example usage
x = [1]
y = [2]
zero = [3]

# Using an array of arrays
target = [x, y, zero]

q = [10, 20, 30]
flags = [true, true, true]

# Apply the function
process_vector(q, flags, target, append_to_vector)

println("x: ", x)
println("y: ", y)
println("zero: ", zero)
