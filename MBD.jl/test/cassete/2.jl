using Cassette

Cassette.@context PrintCtx
const function_calls = []

# Function to log the calls
function log_call(f, args)
    push!(function_calls, (f, args))
    println("Logging call to function: ", f, args)
end

# Hook to intercept function calls
Cassette.prehook(::PrintCtx, f, args...) = log_call(f, args)

# User-defined function 'foo'
function foo(x)
    x + 1
end

# User-defined function 'bar' to demonstrate indirect calls
function bar(x)
    foo(x * 2)
end

# Activate the logging context and execute 'bar'
Cassette.overdub(PrintCtx(), () -> bar(3))

# Display the logged function calls
println("Logged function calls: ", function_calls)
