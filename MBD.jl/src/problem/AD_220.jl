using LinearAlgebra
using JSON
using DelimitedFiles
using Dierckx

function AD220(app)
    if app >= 220 && app < 229
        script_dir = @__DIR__

        # Build the path to the JSON file
        json_path = joinpath(script_dir, "config.json")
        json_data = JSON.parsefile(json_path)
        nb = 2      # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 2      # Number of holonomic constraints
        nhc = 10 # Number of holonomic constraint equations

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
        r_bodies1 = [1., 1.1, 0.]
        si1pr1 = ([0.,0.,0.]-r_bodies1)
        sj2pr1 = [0.,0.,0.]
        if app == 220
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, -ux..., uy..., -ux..., uy...]#y
        end
        if app == 221
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz...]#z
        end
        if app == 222
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, uy..., uz..., uy..., uz...]#
        end
        if app == 223
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, uy..., ux..., uy..., ux...]#x
        end
        if app == 224
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, -uz..., ux..., -uz..., ux...]#x
        end
        if app == 225
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, uz..., uy..., uz..., uy...]#y
        end
        #SMDT = hcat(json_data["SMDT"]["mass1"], json_data["SMDT"]["mass2"])
        
        si1pr1 = ([1.5,1.1,0.]-r_bodies1)
        r_bodies2 = [2., 1.1, 0.]
        sj2pr1 = [0.,0.,0.]-r_bodies2
        SJDT[:, 2] = Any[4, 1, 2, si1pr1..., sj2pr1..., 0, -ux..., uy..., -ux..., uy...]#y
        SMDT = hcat(vcat(30, 90, 90, 90),vcat(30, 90, 90, 90))
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial positions for all bodies (r) and additional parameters (p), with their initial velocities (rd and pd)
        r_bodies = [r_bodies1, r_bodies2]
        p_initial = [1., 0., 0., 0.]  # Common initial 'p' for all bodies
        pd_initial = [0., 0., 0., 0.]  # Common initial 'pd' (velocity) for all 'p'

        # Initial velocities for 'r' bodies - first body has a unique velocity
        rd_initials = [0., 0., 0.]
        # Constructing q0 and qd0 using loops
        q0 = Float64[]
        qd0 = Float64[]

        for i in 1:2
            append!(q0, r_bodies[i]...)
            append!(q0, p_initial...)
            append!(qd0, rd_initials...)
            append!(qd0, pd_initial...)
        end

        contact_mg = Dict(
            "b" => 4,
            "pos" => -999.4630310344,
        )
        contact_lg = Dict(
            "b" => 5,
            "pos" => -11122.7978911861,
        )

        ld_contact =[contact_mg,contact_lg]# [contact_mg, contact_lg]
        damper_mg = Dict(
            "b" => 1,
            "damp" => [0,0,0],
        )
        damper_lg = Dict(
            "b" => 2,
            "damp" => [0,0,0],
        )

        ld_damper = [damper_lg,  damper_mg]#damper_g
        p_contact = Any[ld_damper, []]



        #println(qd0)

        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0, p_contact

    end
end

function AD230(app)
    if app >= 230 && app < 239
        script_dir = @__DIR__

        # Build the path to the JSON file
        json_path = joinpath(script_dir, "config.json")
        json_data = JSON.parsefile(json_path)
        nb = 1      # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 1      # Number of holonomic constraints
        nhc = 5 # Number of holonomic constraint equations

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
        r_bodies = [1., 1.1, 0.]
        si1pr1 = ([0.,0.,0.]-r_bodies)
        sj2pr1 = [0.,0.,0.]
        if app == 220
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, -ux..., uy..., -ux..., uy...]#y
        end
        if app == 221
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz...]#z
        end
        if app == 222
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, uy..., uz..., uy..., uz...]#
        end
        if app == 223
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, uy..., ux..., uy..., ux...]#x
        end
        if app == 224
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, -uz..., ux..., -uz..., ux...]#x
        end
        if app == 225
            SJDT[:, 1] = Any[4, 1, 0, si1pr1..., sj2pr1..., 0, uz..., uy..., uz..., uy...]#y
        end
        #SMDT = hcat(json_data["SMDT"]["mass1"], json_data["SMDT"]["mass2"])
         
        SMDT = hcat(vcat(30, 90, 90, 90))
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial positions for all bodies (r) and additional parameters (p), with their initial velocities (rd and pd)
        r_bodies = [1., 1.1, 0.]
        p_initial = [1., 0., 0., 0.]  # Common initial 'p' for all bodies
        pd_initial = [0., 0., 0., 0.]  # Common initial 'pd' (velocity) for all 'p'

        # Initial velocities for 'r' bodies - first body has a unique velocity
        rd_initials = [0., 0., 0.]
        # Constructing q0 and qd0 using loops
        q0 = Float64[]
        qd0 = Float64[]

        for i in 1:length(1)
            append!(q0, r_bodies...)
            append!(q0, p_initial...)
            append!(qd0, rd_initials...)
            append!(qd0, pd_initial...)
        end

        contact_mg = Dict(
            "b" => 4,
            "pos" => -999.4630310344,
        )
        contact_lg = Dict(
            "b" => 5,
            "pos" => -11122.7978911861,
        )

        ld_contact =[contact_mg,contact_lg]# [contact_mg, contact_lg]
        damper_mg = Dict(
            "b" => 4,
            "damp" => 10000,
        )
        damper_lg = Dict(
            "b" => 5,
            "damp" => 10000,
        )

        ld_damper = [damper_lg,  damper_mg]#damper_g
        p_contact = Any[[], []]



        #println(qd0)

        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0, p_contact

    end
end