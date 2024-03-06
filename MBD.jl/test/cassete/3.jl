using Cassette

Cassette.@context PrintCtx
const unique_function_names = Set{Symbol}()

# Function to log unique function names
function log_unique_function_name(f, args)
    func_name = Symbol(f)
    if !in(func_name, unique_function_names)
        push!(unique_function_names, func_name)
        println("Logging unique function name: ", func_name)
    end
end

# Hook to intercept function calls
Cassette.prehook(::PrintCtx, f, args...) = log_unique_function_name(f, args)

# User-defined function 'foo'
function foo(x)
    x + 1
end

# User-defined function 'bar' to demonstrate indirect calls
function bar(x)
    foo(x * 2)
end

# Another function to demonstrate variety
function baz(x)
    x / 2
end

# Activate the logging context and execute 'bar' and 'baz'
Cassette.overdub(PrintCtx(), () -> bar(3))
Cassette.overdub(PrintCtx(), () -> baz(6))

# Display the unique function names logged
println("Logged unique function names: ", unique_function_names)
