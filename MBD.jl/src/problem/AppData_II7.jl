using LinearAlgebra
include("../mathfunction_II7.jl")
include("AppData_door.jl")
export AppData_II7,process_vector,AppDataStruct

struct AppDataStruct
    name::String
    nb::Int
    ngc::Int
    nh::Int
    nc::Int
    NTSDA::Int
    SJDT::Array{Any}
    SMDT::Array
    STSDAT
    q0::Array
    qd0::Array
end

function AppData_II7(app)
    if app == 1 || app == 101 || app == 102 || app == 103 || app == 201  # Pendulum, Spherical to Ground
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
        sipr = uy
        sjpr = [0,0 , 0]
        SJDT[:, 1] = Any[2, 1, 0, sipr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground

        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        SMDT = [30, 90, 90, 30]

        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [0, -1, 0]
        if app == 1 || app == 101
            p10 = [0, ux...]
        end
        if app == 102
            p10 = [0, uy...]
        end
        if app == 103
            p10 = [0, uz...]
        end
        if app == 201
            p10 = [1, zer...]
        end
        q0 = [r10..., p10...]
        omeg1pr0 = [0, 0, 0]
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * sipr
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        qd0 = [r1d0..., p1d0...]
        w1 = 2 * mathfunction.EEval(p10) * p1d0
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end

    if app == 2  # Spin Stabilized Top, Spherical to Ground

        nb = 1  # Number of bodies
        ngc = 7 * nb  # Number of generalized coordinates
        nh = 1  # Number of holonomic constraints
        nhc = 3  # Number of holonomic constraint equations
        nc = nhc + nb  # Number of constraint equations
        nv = ngc - nc
        nu = nc
        NTSDA = 0  # Number of TSDA force elements

        ux = [1, 0, 0]
        uy = [0, 1, 0]
        uz = [0, 0, 1]
        zer = zeros(3)
        SJDT = Array{Any}(undef, 22, nh)
        # SJDT(22,nh): Joint Data Table
        # SJDT(:,k)=[t; i; j; sipr; sjpr; d; uxipr; uzipr; uxjpr; uzjpr]; 
        #SJDT = zeros(22, nh)
        SJDT[:, 1] = Any[2, 1, 0, -uz..., zer..., 0, zer..., zer..., zer..., zer...]  # Sph Jt - Body1 and ground

        # SMDT(4,nb): Mass Data Table (With diagonal inertia matrix)
        SMDT = [30, 90, 90, 90]

        # STSDAT(12,1): TSDA Data Table
        if NTSDA == 0
            STSDAT = zeros(12, NTSDA)
        end

        # Initial generalized coordinates
        r10 = [0, 0, 1]
        p10 = [1, 0, 0, 0]
        q0 = [r10..., p10...]
        omeg1pr0 = [1e-12; 1e-12; 13.5]
        r1d0 = mathfunction.atil(omeg1pr0) * uz
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
        STSDAT =
            NTSDA == 0 ? zeros(12, NTSDA) : # STSDAT initialization logic here if NTSDA != 0

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

    if app == 202  # doubleb Pendulum, Spherical to Ground
        nb = 2          # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 2          # Number of holonomic constraints
        nhc = 6         # Number of holonomic constraint equations
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
        si1pr = [-1, 0, 0]
        sjpr = [0, 0, 0]
        SJDT[:, 1] = Any[2, 1, 0, si1pr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        si2pr = [0, 0, 0]
        sjpr = [1, 0, -2]
        SJDT[:, 2] = Any[2, 1, 2, si2pr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground

        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        SMDT = hcat(vcat(30, 90, 90, 30), vcat(30, 90, 90, 30))

        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [1, 0, 0]
        p10 = [1, zer...]
        r20 = [0, 0, 2]
        p20 = [1, zer...]
        q0 = [r10..., p10..., r20..., p20...]
        omeg1pr0 = [0, 0, 0]
        omeg2pr0 = [0, 0, 0]
        #A=mathfunction.ATran(p10)
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        r2d0 = mathfunction.ATran(p20) * mathfunction.atil(omeg2pr0) * si2pr
        p2d0 = 0.5 * mathfunction.GEval(p20)' * omeg2pr0
        qd0 = [r1d0..., p1d0..., r2d0..., p2d0...]
        #w1=2*mathfunction.EEval(p10)*p1d0
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end

    if app == 203  # double rev
        nb = 2          # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 2          # Number of holonomic constraints
        nhc = 10         # Number of holonomic constraint equations
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
        si1pr = [-1, 0, 0]
        sjpr = [0, 0, 0]
        SJDT[:, 1] = Any[4, 1, 0, si1pr..., sjpr..., 0, ux..., uz..., ux..., uz...]  # Spherical Joint - Body 1 and Ground
        si2pr = [0, 0, 0]
        sjpr = [1, 0, -2]
        SJDT[:, 2] = Any[4, 1, 2, si2pr..., sjpr..., 0, uz..., ux..., uz..., ux...]  # Spherical Joint - Body 1 and Ground

        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        SMDT = hcat(vcat(30, 90, 90, 30), vcat(30, 90, 90, 30))

        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [1, 0, 0]
        p10 = [1, zer...]
        r20 = [0, 0, 2]
        p20 = [1, zer...]
        q0 = [r10..., p10..., r20..., p20...]
        omeg1pr0 = [0, 0, 0]
        omeg2pr0 = [0, 0, 0]
        #A=mathfunction.ATran(p10)
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        r2d0 = mathfunction.ATran(p20) * mathfunction.atil(omeg2pr0) * si2pr
        p2d0 = 0.5 * mathfunction.GEval(p20)' * omeg2pr0
        qd0 = [r1d0..., p1d0..., r2d0..., p2d0...]
        #w1=2*mathfunction.EEval(p10)*p1d0
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end

    if app == 204  # doubleb Pendulum,+plainer  Spherical to Ground
        nb = 1         # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 1         # Number of holonomic constraints
        nhc = 3        # Number of holonomic constraint equations
        nc = nhc + nb   # Number of constraint equations
        NTSDA = 0       # Number of TSDA force elements

        ux = [1, 0, 0]
        uy = [0, 1, 0]
        uz = [0, 0, 1]
        zer = zeros(3)

        # SJDT(22, nh): Joint Data Table
        # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
        # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
        # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
        #SJDT = zeros(22, nh)
        SJDT = []
        SJDT = Array{Any}(undef, 22, nh)
        si1pr = [1, 0, 1]
        sjpr = [0, 0, 0]
        SJDT[:, 1] = Any[2, 1, 0, si1pr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        si2pr = [0, 0, 1]
        sjpr = [0, 0, 0]
        #SJDT[:, 2] = Any[2, 1, 0, si2pr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        si3pr = [0, 0, 0]
        sj3pr = [1, 0, 0]
        #SJDT[:, 2] = Any[1010, 1, 0, si3pr..., sj3pr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        #SJDT[:, 2] = Any[1, 1, 0, si3pr..., sj3pr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        #SMDT = hcat(vcat(30, 90, 90, 30), vcat(30, 90, 90, 30))
        SMDT = hcat(vcat(30, 90, 90, 90))
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [-1, 0, -1]
        p10 = [1, zer...]
        r20 = [0, 0, -1]
        p20 = [1, zer...]
        #q0 = [r10..., p10...,r20..., p20...]
        q0 = [r10..., p10...]
        omeg1pr0 = [0, 0, 0]
        omeg2pr0 = [0, 0, 0]
        #A=mathfunction.ATran(p10)
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        r2d0 = mathfunction.ATran(p20) * mathfunction.atil(omeg2pr0) * si2pr
        p2d0 = 0.5 * mathfunction.GEval(p20)' * omeg2pr0
        #qd0 = [r1d0..., p1d0...,r2d0..., p2d0...]
        qd0 = [r1d0..., p1d0...]
        #w1=2*mathfunction.EEval(p10)*p1d0
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end
    if app == 205  # single Pendulum+plainer x=1  Spherical to Ground
        nb = 1         # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 2        # Number of holonomic constraints
        nhc = 4        # Number of holonomic constraint equations
        nc = nhc + nb   # Number of constraint equations
        NTSDA = 0       # Number of TSDA force elements

        zer = zeros(3)
        # SJDT(22, nh): Joint Data Table
        # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
        # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
        # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
        SJDT = Array{Any}(undef, 22, nh)
        si1pr = [-1, -1, 0]
        sjpr = [0, 0, 0]
        SJDT[:, 1] = Any[2, 1, 0, si1pr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        si3pr = [1, 1, 0]
        sj3pr = [0, 0, 0]
        SJDT[:, 2] = Any[1010, 1, 0, si3pr..., sj3pr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        SMDT = hcat(vcat(30, 90, 90, 90))
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [1, 1, 0]
        p10 = [1, zer...]

        q0 = [r10..., p10...]
        omeg1pr0 = [0, 0, 0]
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        qd0 = [r1d0..., p1d0...]
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end
    if app == 206  # single Pendulum+plainer x=0  Spherical to Ground
        nb = 1         # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 2        # Number of holonomic constraints
        nhc = 4       # Number of holonomic constraint equations
        nc = nhc + nb   # Number of constraint equations
        NTSDA = 0       # Number of TSDA force elements

        zer = zeros(3)
        # SJDT(22, nh): Joint Data Table
        # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
        # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
        # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
        SJDT = Array{Any}(undef, 22, nh)
        si1pr = [0, 0, 1]
        sjpr = [0, 0, 0]
        SJDT[:, 1] = Any[2, 1, 0, si1pr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        si3pr = [0, 0, 0]
        sj3pr = [0, 0, 0]
        SJDT[:, 2] = Any[1010, 1, 0, si3pr..., sj3pr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        #SMDT = hcat(vcat(30, 90, 90, 30), vcat(30, 90, 90, 30))
        SMDT = hcat(vcat(30, 90, 90, 90))
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [0, 0, 0]
        p10 = [0,1,0,0]

        #q0 = [r10..., p10...,r20..., p20...]
        q0 = [r10..., p10...]
        omeg1pr0 = [0, 0, 0]
        #A=mathfunction.ATran(p10)
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        #qd0 = [r1d0..., p1d0...,r2d0..., p2d0...]
        qd0 = [r1d0..., p1d0...]
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end
    if app == 207  # single Pendulum+plainer
        nb = 1         # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 2        # Number of holonomic constraints
        nhc = 4       # Number of holonomic constraint equations
        nc = nhc + nb   # Number of constraint equations
        NTSDA = 0       # Number of TSDA force elements

        zer = zeros(3)
        # SJDT(22, nh): Joint Data Table
        # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
        # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
        # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
        SJDT = Array{Any}(undef, 22, nh)
        si1pr = [-1, 0, -1]
        sjpr = [0, 0, 0]
        SJDT[:, 1] = Any[2, 1, 0, si1pr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        si3pr = [0, 0, 0]
        sj3pr = [0, 0, 0]
        SJDT[:, 2] = Any[1010, 1, 0, si3pr..., sj3pr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        #SMDT = hcat(vcat(30, 90, 90, 30), vcat(30, 90, 90, 30))
        SMDT = hcat(vcat(30, 90, 90, 90))
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [1, 0, 1]
        p10 = [0,1,0,0]

        #q0 = [r10..., p10...,r20..., p20...]
        q0 = [r10..., p10...]
        omeg1pr0 = [0, 0, 0]
        #A=mathfunction.ATran(p10)
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        #qd0 = [r1d0..., p1d0...,r2d0..., p2d0...]
        qd0 = [r1d0..., p1d0...]
        println("q0",q0)
        println("qd0",qd0)
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    
    end

    if app == 208  # single Pendulum+plainer x=1  Spherical to Ground *
        nb = 1         # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 2        # Number of holonomic constraints
        nhc = 4        # Number of holonomic constraint equations
        nc = nhc + nb   # Number of constraint equations
        NTSDA = 0       # Number of TSDA force elements

        zer = zeros(3)
        # SJDT(22, nh): Joint Data Table
        # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
        # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
        # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
        SJDT = Array{Any}(undef, 22, nh)
        si1pr = [-1, -1, 0]
        sjpr = [0, 0, 0]
        SJDT[:, 1] = Any[2, 1, 0, si1pr..., sjpr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        si3pr = [0, 0, 0]
        sj3pr = [0, 0, 0]
        SJDT[:, 2] = Any[1060, 1, 0, si3pr..., sj3pr..., [(1, -0.5), (2, 0.5), (3, 1)], zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
        # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        #SMDT = hcat(vcat(30, 90, 90, 30), vcat(30, 90, 90, 30))
        SMDT = hcat(vcat(30, 0.1, 0.1, 0.1))
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial generalized coordinates
        r10 = [1, 1, 0]
        p10 = [1, zer...]

        #q0 = [r10..., p10...,r20..., p20...]
        q0 = [r10..., p10...]
        omeg1pr0 = [0, 0, 0]
        #A=mathfunction.ATran(p10)
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        #qd0 = [r1d0..., p1d0...,r2d0..., p2d0...]
        qd0 = [r1d0..., p1d0...]
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0
    end


    if app == 301  # single Pendulum+plainer
        apps = model_sph_plain()
        println0(apps)
        return apps.nb, apps.ngc, apps.nh, apps.nc, apps.NTSDA,
         apps.SJDT, apps.SMDT, apps.STSDAT, apps.q0, apps.qd0
    end
    if app == 302  # door
        apps = model_door()
        println0(apps)
        return apps.nb, apps.ngc, apps.nh, apps.nc, apps.NTSDA,
         apps.SJDT, apps.SMDT, apps.STSDAT, apps.q0, apps.qd0
    end
    if app == 303  # space 4
        apps = model_cr()
        println0(apps)
        return apps.nb, apps.ngc, apps.nh, apps.nc, apps.NTSDA,
         apps.SJDT, apps.SMDT, apps.STSDAT, apps.q0, apps.qd0
    end
    if app == 304  # door2
        apps = model_door_2()
        println0(apps)
        return apps.nb, apps.ngc, apps.nh, apps.nc, apps.NTSDA,
         apps.SJDT, apps.SMDT, apps.STSDAT, apps.q0, apps.qd0
    end
    if app == 305  # door3 *
        apps = model_door_305()
        println0(apps)
        return apps.nb, apps.ngc, apps.nh, apps.nc, apps.NTSDA,
         apps.SJDT, apps.SMDT, apps.STSDAT, apps.q0, apps.qd0
    end
    if app == 306  # model_door_306_poly
        apps = model_door_306_poly()
        println0(apps)
        return apps.nb, apps.ngc, apps.nh, apps.nc, apps.NTSDA,
         apps.SJDT, apps.SMDT, apps.STSDAT, apps.q0, apps.qd0
    end
end

function println0(apps::AppDataStruct)
    println("App Data")
    println("========")
    println("Name: ", apps.name)
    println("Number of Bodies (nb): ", apps.nb)
    println("Number of Geometric Constraints (ngc): ", apps.ngc)
    println("Number of Holonomic Constraints (nh): ", apps.nh)
    println("Number of Constraints (nc): ", apps.nc)
    println("\nSJDT :")
    for row in eachrow(apps.SJDT')
        # println(row)
    end
    println("\nSMDT :")
    println(apps.SMDT)
    println("\nInitial (q0): ", apps.q0)
    println("Initial Derivative (qd0): ", apps.qd0)
    println("========")
end

function process_vector(q::Vector, flags::Vector{Bool}, target::Vector{Vector{Float64}}, func::Function)
    # Ensure `q` and `flags` have the same length
    # if length(q) != length(flags)
    #     error("Vectors 'q' and 'flags' must be of the same length.")
    # end

    # Process elements of `q` based on `flags`
    for (i, flag) in enumerate(flags)
        if flag && i <= length(target)
            func(target[i], q[i])
        end
    end
end

# Function to append `element` of `q` to the `target_vector`
function append_to_vector(target_vector::Vector{Float64}, element)
    push!(target_vector, element)
end

function plot_vector(target_vector::Vector{Float64}, element)
    push!(target_vector, element)
end


function model_sph_plain()
    nb = 1         # Number of bodies
    ngc = 7 * nb    # Number of generalized coordinates
    nh = 2        # Number of holonomic constraints
    nhc = 4       # Number of holonomic constraint equations
    nc = nhc + nb   # Number of constraint equations
    NTSDA = 0       # Number of TSDA force elements

    zer = zeros(3)
    
    sphroot1=[0,0,0]
    sphroot2=[2,0,0]
    sphroot3=[1,1,0]
    r10_ = [0.5, 0.5, 0]
    r20_ = [1.5, 0.5, 0]
    r00_ = [0  , 0  , 0]
    # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
    # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
    # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
    # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
    SJDT = Array{Any}(undef, 22, nh)
    si1pr1 = ((sphroot1-r10_)'*[1  0  0;0  1  0;0  0  1])'
    sj2pr1 = ((sphroot1-r00_)'*[1  0  0;0  1  0;0  0  1])'
    SJDT[:, 1] = Any[2, 1, 0, si1pr1..., sj2pr1..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground

    #si1pr2 = ((sphroot2-r20_)'*[1  0  0;0  -1  0;0  0  1])'
    #sj2pr2 = ((sphroot2-r00_)'*[1  0  0;0  -1  0;0  0  1])'
    #SJDT[:, 2] = Any[2, 2, 0, si1pr2..., sj2pr2..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground

    #si1pr3 = ((sphroot3-r10_)'*[1  0  0;0  -1  0;0  0  1])'
    #sj2pr3 = ((sphroot3-r20_)'*[1  0  0;0  -1  0;0  0  1])'
    #SJDT[:, 2] = Any[2, 1, 2, si1pr3..., sj2pr3..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground

    #SJDT[:, 4] = Any[1010, 1, 0, si3pr..., sj3pr..., 0, zer..., zer..., zer..., zer...]  # Spherical Joint - Body 1 and Ground
    
    si1pr4 = ((sphroot3-r10_)'*[1  0  0;0  1  0;0  0  1])'
    sj2pr4 = [0,0,0]
    si1pr4= [0.5,0.5,0]
    SJDT[:, 2] = Any[1020, 1, 0, si1pr4..., zer..., 1, zer..., zer..., zer..., zer...]
    # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
    # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
    #SMDT = hcat(vcat(30, 90, 90, 90),vcat(30, 90, 90, 90))
    SMDT = hcat(vcat(30, 90, 90, 90))
    # STSDAT(12, 1): TSDA Data Table
    STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

    # Initial generalized coordinates

    p10 = [1, 0, 0, 0]
    p20 = [1, 0, 0, 0]
    #q0 = [r10_..., p10...,r20_..., p20...]
    q0 = [r10_..., p10...]
    omeg1pr0 = [0, 0, 0]
    r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr1
    p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
    #r2d0 = mathfunction.ATran(p20) * mathfunction.atil(omeg1pr0) * si1pr3
    #p2d0 = 0.5 * mathfunction.GEval(p20)' * omeg1pr0
    #r1d0 = [0, 0, 0]
    #r2d0 = [0, 0, 0]
    #p1d0 = [0, 1, 0, 0]
    #p2d0 = [0, 1, 0, 0] 
    #qd0 = [r1d0..., p1d0...,r2d0..., p2d0...]
    qd0 = [r1d0..., p1d0...]
    println(qd0)
    return AppDataStruct("model_sph_plain",nb,ngc,nh,nc,NTSDA,SJDT,SMDT,STSDAT,q0,qd0)
end

function model_cr()
    nb = 3         # Number of bodies
    ngc = 7 * nb    # Number of generalized coordinates
    nh = 2   +2     # Number of holonomic constraints
    nhc = 6   +10    # Number of holonomic constraint equations
    nc = nhc + nb   # Number of constraint equations
    NTSDA = 0       # Number of TSDA force elements

    ux = [1; 0; 0]
    uy = [0; 1; 0]
    uz = [0; 0; 1]
    zer = zeros(3)
    revroot=[0.0, 0.1, 0.12]
    sphroot1=[0.0, 0.1, 0.2]
    sphroot2=[0.2, 0.0, 0.0]
    tranroot=[0.2, 0.0, 0.0]

    r1 = [0.0, 0.1, 0.12]
    r2 =   [0.1,	0.05,	0.1]
    r3 =    [0.2,0,0]
    r0 = [0  , 0  , 0]
    # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
    # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
    # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
    # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
    SJDT = Array{Any}(undef, 22, nh)
    si1pr1 = (revroot-r1)
    sj2pr1 = (revroot-r0)
    SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, uz..., ux..., uz..., ux...]  

    si1pr2 = (sphroot1-r1)
    sj2pr2 = (sphroot1-r2)
    SJDT[:, 2] = Any[2, 1, 2, si1pr2..., sj2pr2..., 0, ux..., uz..., ux..., uz...]  

    si1pr3 = (sphroot2-r2)
    sj2pr3 = (sphroot2-r3)
    SJDT[:, 3] = Any[2, 2, 3, si1pr3..., sj2pr3..., 0, ux..., uz..., ux..., uz...]  
    si1pr4 = (tranroot-r3)
    sj2pr4 = (tranroot-r0)
    SJDT[:, 4] = Any[5, 3, 0, si1pr4..., sj2pr4..., 0, uz..., ux..., uz..., ux...]  

    # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
    # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
    SMDT = hcat(vcat(30, 90, 90, 90),vcat(30, 90, 90, 90),vcat(30, 90, 90, 90))
    #SMDT = hcat(vcat(30, 90, 90, 90))
    # STSDAT(12, 1): TSDA Data Table
    STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

    # Initial generalized coordinates

    p10 = [1, 0, 0, 0]
    p20 = [1, 0, 0, 0]
    p30 = [1, 0, 0, 0]

    q0 = [r1..., p10...,r2..., p20...,r3..., p30...]

    omeg1pr0 = [3.14/4, 0, 0]
    r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr1
    p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
    #r2d0 = mathfunction.ATran(p20) * mathfunction.atil(omeg1pr0) * si1pr3
    #p2d0 = 0.5 * mathfunction.GEval(p20)' * omeg1pr0
    #r1d0 = [0, 0, 0]
    r2d0 = [0, 0, 0]
    r3d0 = [0.03, 0, 0]
    #p1d0 = [0, 0, 0, 0]
    p2d0 = [0, 0, 0, 0]
    p3d0 = [0, 0, 0, 0] 
    qd0 = [r1d0..., p1d0...,r2d0..., p2d0...,r3d0..., p3d0...]

    println(qd0)
    return AppDataStruct("model_cr",nb,ngc,nh,nc,NTSDA,SJDT,SMDT,STSDAT,q0,qd0)
end

