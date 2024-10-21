using OrdinaryDiffEq, LinearSolve, SparseArrays

# Define your ODE problem
function f!(du, u, p, t)
    du[1] = u[2]
    du[2] = -u[1]
end

u0 = [1.0, 0.0]
tspan = (0.0, 10.0)
prob = ODEProblem(f!, u0, tspan)

# Create a matrix (for example, a Jacobian matrix)
# If it's dense, convert it to sparse
A_dense = rand(10, 10)  # Example dense matrix
A_sparse = sparse(A_dense)  # Convert to sparse matrix

# Define the KLU factorization solver for sparse matrices
linsolver = QRFactorization()

# Set up the Rodas5 solver with autodiff=true and KLUFactorization as the linear solver
solver = Rodas5(autodiff=true, linsolve=linsolver)

# Solve the ODE problem
sol = solve(prob, solver, reltol=1e-6, abstol=1e-6)

# Now you can work with the solution `sol`

    # solver_dict = Dict(
    #     "Rodas5_autodiff_true" => Rodas5(autodiff=true),
    #     "Rodas5_autodiff_false_forward" => Rodas5(autodiff=false, diff_type=Val{:forward}),
    #     "Rodas5_autodiff_false_central" => Rodas5(autodiff=false, diff_type=Val{:central}),
    #     "Rodas5_autodiff_false_complex" => Rodas5(autodiff=false, diff_type=Val{:complex}),
    #     "Rodas5_autodiff_trueQR" => Rodas5(autodiff=true,linsolve = QRFactorization()),
    #     "Rodas5_autodiff_false_forwardQR" => Rodas5(autodiff=false, diff_type=Val{:forward},linsolve = QRFactorization()),
    #     "Rodas5_autodiff_false_centralQR" => Rodas5(autodiff=false, diff_type=Val{:central},linsolve = QRFactorization()),
    #     "Rodas5_autodiff_false_complexQR" => Rodas5(autodiff=false, diff_type=Val{:complex},linsolve = QRFactorization())
    # )