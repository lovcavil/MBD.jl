using LinearAlgebra
using JSON
using DelimitedFiles
using Dierckx


function AD240(app)
    if app >= 240 && app < 249
        script_dir = @__DIR__

        # Build the path to the JSON file
        json_path = joinpath(script_dir, "config.json")
        json_data = JSON.parsefile(json_path)
        nb = 3      # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 1      # Number of holonomic constraints
        nhc = 1*5 # Number of holonomic constraint equations

        nc = nhc + nb   # Number of constraint equations

        NTSDA = 0       # Number of TSDA force elements

        ux = [1; 0; 0]
        uy = [0; 1; 0]
        uz = [0; 0; 1]
        zer = zeros(3)



        # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
        # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
        # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
        SJDT = Array{Any}(undef, 22, nh)
        r_bodies1 = [0.1, 0.0, 0.1]
        r_bodies2 = [0.475, 0.1, 0.0]
        r_bodies3 = [0.75, 0.0, 0.0]
        si1pr1 = [0.,0.,0.]-r_bodies1
        sj2pr1 = [0.,0.,0.]

        
        SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz...]#z
        #SJDT[:, 2] = Any[10, 1, 0, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz...]#z
        #SJDT[:, 2] = Any[4, 1, 2, si1pr2..., sj2pr2..., 0, ux..., uz..., ux..., uz...]#z
        #SJDT[:, 3] = Any[4, 2, 3, si1pr3..., sj2pr3..., 0, ux..., uz..., ux..., uz...]#z

        #SMDT = hcat(json_data["SMDT"]["mass1"], json_data["SMDT"]["mass2"])
         
        SMDT = hcat(vcat(1.0, 10.E-03, 10.E-03, 10.E-05),
        vcat(1.0, 1.0, 1.0, 10.E-03),
        vcat(5.0, 10E-02, 10E-02, 10E-02)
        )
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        p10=[1., 0., 0., 0.]
        omeg1pr0 = [0, 0, 1]
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr1
        p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0

        # Initial positions for all bodies (r) and additional parameters (p), with their initial velocities (rd and pd)
        q0 = [r_bodies1..., 1., 0., 0., 0.,
        0.475, 0.1, 0.0, 1., 0., 0., 0.,
        0.75, 0.0, 0.0,1., 0., 0., 0.,
        ]

        qd0 = [r1d0..., p1d0...,
            0., 0., 0., 0.,0., 0., 0., 
            0., 0., 0., 0.,0., 0., 0.
        ]  

        ld_contact = []
        ld_damper = []#damper_g
        p_contact = Any[ld_damper,ld_contact,0]

        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0, p_contact

    end
end