using LinearAlgebra
include("mathfunction_II7.jl")
export AppData_II7

function AppData_II7(app)
    if app == 1  # Pendulum, Spherical to Ground
        nb = 1          # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 1          # Number of holonomic constraints
        nhc = 3         # Number of holonomic constraint equations
        nc = nhc + nb   # Number of constraint equations
        NTSDA = 0       # Number of TSDA force elements

        ux = [1, 0, 0]
        uy = [0, 1, 0]
        uz = [0, 0, 1]
        zer = zeros(3)

        # SJDT(22, nh): Joint Data Table
        # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
        # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0;
        # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
        #SJDT = zeros(22, nh)
        SJDT = []
        SJDT = Array{Any}(undef, 22, nh)
        SJDT[:, 1] = Any[2, 1, 0, -uz..., zer..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground

        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        SMDT = [30, 90, 90, 30]

        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [0, 0, -1]
        p10 = [0, ux...]
        q0 = [r10..., p10...]
        omeg1pr0 = [2, 0, 0]
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * uz
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        qd0 = [r1d0..., p1d0...]

        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end

    if app == 6   # Spatial Slider-Crank
    
        nb = 2       # Number of bodies
        ngc = 7 * nb # Number of generalized coordinates
        nh = 3       # Number of holonomic constraints
        nhc = 11     # Number of holonomic constraint equations
        nc = nhc + nb # Number of constraint equations
        NTSDA = 0    # Number of TSDA force elements
    
        ux = [1; 0; 0]
        uy = [0; 1; 0]
        uz = [0; 0; 1]
        zer = zeros(3, 1)
    
        # SJDT(22, nh): Joint Data Table 
        # SJDT[:, k] = [t; i; j; sipr; sjpr; d; uxipr; uzipr; uxjpr; uzjpr]
        SJDT = Array{Any}(undef, 22, nh)
        SJDT[:, 1] = Any[4, 1, 0, zer..., (0.1 * ux + 0.12 * uy)..., 0, ux..., uz..., ux..., uz...] # Cyl-Body1 to Ground
        SJDT[:, 2] = Any[5, 2, 0, zer..., zer..., 0, ux..., uz..., ux..., uz...] # Tran-Body2 to Ground
        SJDT[:, 3] = Any[1, 1, 2, (0.08 * uy)..., (0.02 * uy)..., 0.23, zer..., zer..., zer..., zer...] # Dist-Body 1 to 2
    
        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        SMDT = hcat(vcat(0.5, 0.2, 0.2, 0.2), vcat(5, 0.2, 0.2, 0.2))
    
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : # STSDAT initialization logic here if NTSDA != 0
    
        # Initial generalized coordinates
        r10 = [0.1, 0.12, 0]
        r1111 = [0.1, 0.12, 0]
        p10 = [1, 0, 0, 0]
        r20 = [0, 0, 0.1027]
        p20 = [1, 0, 0, 0]
        q0 = [r1111..., p10..., r20..., p20...]
        r1d0 = [0, 0, 0]
        p1d0 = 0.5 * mathfunction.EEval(p10)' * 120 * uz
        r2d0 = [0, 0, 9.378]
        p2d0 = [0, 0, 0, 0]
        qd0 = [r1d0..., p1d0..., r2d0..., p2d0...]
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end
    
end
