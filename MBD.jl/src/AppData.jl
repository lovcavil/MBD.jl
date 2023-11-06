using LinearAlgebra
include("mathfunction.jl")
export AppData
function AppData(app)

    z3 = zeros(3)
    zer = z3
    z4 = zeros(4)
    I3 = Matrix{Float64}(I, 3, 3)
    ux = [1, 0, 0]
    uy = [0, 1, 0]
    uz = [0, 0, 1]

    # Initialize variables to hold the return values
    nb = ngc = nh = nhc = nc = nd = 0
    SJDT = []
    q0e = []

    # App 1: 2 Bar with 2 rotational drivers
    if app == 1
        nb = 2        # Number of bodies
        ngc = 7 * nb  # Number of generalized coordinates
        nh = 2        # Number of holonomic constraints
        nhc = 10      # Number of holonomic constraint equations
        nc = nhc + nb # Number of constraint equations
        nd = ngc - nc
        
        # Spatial Joint Data Table
        #SJDT = zeros(Float64, 22, nc+nd)
        #SJDT=Any[]
        SJDT = Array{Any}(undef, 22, nc+nd)
        # Create an empty array of type Any to hold the transposed data
        #transposed_SJDT = Any[]

        # Push arrays representing columns into the transposed array
        #push!(transposed_SJDT, Any[4, 1, 0, z3...,z3...,  0,  ux..., uz..., ux..., uz...])
        #push!(transposed_SJDT, Any[4, 1, 2, ux..., z3..., 0, -uz..., ux..., ux..., uz...])
        #push!(transposed_SJDT, Any[10, 1, 0, z3...,z3..., 0,  ux..., uz..., ux..., uz...])
        #push!(transposed_SJDT, Any[10, 1, 2, ux..., z3...,0, -uz..., ux..., ux..., uz...])
        SJDT[:, 1] = Any[4, 1, 0, z3...,z3...,  0,  ux..., uz..., ux..., uz...]
        SJDT[:, 2] = Any[4, 1, 2, ux..., z3..., 0, -uz..., ux..., ux..., uz...]
        SJDT[:, 3] = Any[10, 1, 0, z3...,z3..., 0,  ux..., uz..., ux..., uz...]
        SJDT[:, 4] = Any[10, 1, 2, ux..., z3...,0, -uz..., ux..., ux..., uz...]

        # Now transpose the array to get the original column-wise structure
        #SJDT = permutedims(hcat(transposed_SJDT...), (2, 1))
        #SJDT = permutedims(hcat(SJDT...), (2, 1))
        
        # Body Initial Configuration Data Estimate
        BPDDT = zeros(Float64, 9, 2)
        BPDDT[:,1] = [z3..., ux..., uy...]
        BPDDT[:,2] = [   ux...; (ux-uz)...; (ux+uy)...]
        BPDD=BPDDT'
        # Initial generalized coordinate estimate
        q0e = zeros(7 * nb)
        for j in 1:nb
            rO = BPDDT[1:3, j]
            rP = BPDDT[4:6, j]
            rQ = BPDDT[7:9, j]
            q, p, A = mathfunction.InitConfig(rO, rP, rQ) # Needs to be defined or imported
            q0e[(7*(j-1)+1):(7*j)] = q
        end
    end

    return nb, ngc, nh, nhc, nc, nd, SJDT, q0e
end
