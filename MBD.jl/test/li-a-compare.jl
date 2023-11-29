using LinearAlgebra
using BenchmarkTools
using IterativeRefinement
# Define a placeholder for the IterativeRefinement method
# Replace this with your actual implementation or package function


# Generate a random matrix A and a vector b
n = 100  # Size of the matrix
A = randn(n, n)
b = randn(n)

# Basic method using backslash operator
function solve_basic(A, b)
    return A \ b
end

# Direct inversion
function solve_direct_inv(A, b)
    return inv(A) * b
end

# Iterative Refinement method
function solve_iterative_refinement(A, b)
    w, bnorm, bcomp = rfldiv(A,b)
    return w
end
#w = -JJ \ Resid
#w, bnorm, bcomp = rfldiv(-JJ,Resid)
#w = -inv(JJ) * Resid
# Benchmarking
println("Benchmarking basic method:")
@btime solve_basic($A, $b)

println("\nBenchmarking direct inversion:")
@btime solve_direct_inv($A, $b)

println("\nBenchmarking iterative refinement method:")
@btime solve_iterative_refinement($A, $b)
