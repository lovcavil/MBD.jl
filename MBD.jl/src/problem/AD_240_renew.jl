using LinearAlgebra
using JSON
using DelimitedFiles
using Dierckx


function AD240(app,contact_json)
    if app >= 240 && app < 245
        script_dir = @__DIR__

        # Build the path to the JSON file
        json_path = joinpath(script_dir, "config.json")
        json_data = JSON.parsefile(json_path)
        nb = 2      # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 3      # Number of holonomic constraints
        nhc = 11 # Number of holonomic constraint equations
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
        r_bodies1 = [0.5, 0.0, 0.5]
        r_bodies2 = [1.5, 0.0, 0.5]
        r_bodies3 = [1.5, 0.0, 0.0] 
        si1pr1 = [0.,0.,0.]-r_bodies1
        sj2pr1 = [0.,0.,0.]-zer
        si1pr2 = [1.,0.,1.]-r_bodies1
        sj2pr2 = [1.,0.,1.]-r_bodies2

        # SJDT(:,1)=[4;1;0;zer;0.1*ux+0.12*uy;0;ux;uz;ux;uz];   %Cyl-Body1 to Ground
        # SJDT(:,2)=[5;2;0;zer;zer;0;ux;uz;ux;uz];              %Tran-Body2 to Ground
        # SJDT(:,3)=[1;1;2;0.08*uy;0.02*uy;0.23;zer;zer;zer;zer];  %Dist-Body 1 to 2


        SJDT[:, 1] = Any[4, 1, 0, zer..., 0.1,0.12,0, 0, ux..., uz..., ux..., uz...]#z
        SJDT[:, 2] = Any[5, 2, 0, zer..., zer..., 0, ux..., uz..., ux..., uz...]#z
        SJDT[:, 3] = Any[1, 1, 2, 0,0.08,0,0,0.02,0, 0.23, zer..., zer..., zer..., zer...]#z

        #SMDT = hcat(json_data["SMDT"]["mass1"], json_data["SMDT"]["mass2"])
         
        SMDT = hcat(vcat(5.0, 0.2,0.2,0.2),
        vcat(5.0, 0.2,0.2,0.2),
        #vcat(5.0, 10E-02, 10E-02, 10E-02)
        )
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        p10=[1., 0., 0., 0.]
        omeg1pr0 = [0, 0, 0]
        r1d0 = mathfunction.ATran(p10) * mathfunction.atil(omeg1pr0) * si1pr1
        r1d0=[0, 0, 0]
        # p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        p1d0=0.5 * mathfunction.EEval(p10)' * 120 * uz
        r2d0=[0,0,9.378]

        # Initial positions for all bodies (r) and additional parameters (p), with their initial velocities (rd and pd)
        q0 = [0.1,0.12,0, 1., 0., 0., 0.,
        0,0,0.1027, 1., 0., 0., 0.,
        #r_bodies3...,1., 0., 0., 0.,
        ]

        qd0 = [r1d0..., p1d0...,
            r2d0..., 0.,0., 0., 0., 
        #    0., 0., 0., 0.,0., 0., 0.
        ]  

        ld_contact = []
        ld_damper = []#damper_g
        p_contact = Any[ld_damper,ld_contact,0,0]

        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0, p_contact

    end
    if app >=  245 && app < 246
        script_dir = @__DIR__

        # Build the path to the JSON file
        json_path = joinpath(script_dir, "config.json")
        json_data = JSON.parsefile(json_path)
        nb = 3      # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 5      # Number of holonomic constraints
        nhc = 17 # Number of holonomic constraint equations
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
        r_bodies1 = [0.1,0.12,0]
        r_bodies2 = [0,0.02,0.1027]
        r_bodies3 = [0,0,0.1027] 
        si1pr3 = [0,0,0.1027]-r_bodies2
        sj2pr3 = [0,0,0.1027]-r_bodies3
        #si1pr2 = [1.,0.,1.]-r_bodies1
        #sj2pr2 = [1.,0.,1.]-r_bodies2

        # SJDT(:,1)=[4;1;0;zer;0.1*ux+0.12*uy;0;ux;uz;ux;uz];   %Cyl-Body1 to Ground
        # SJDT(:,2)=[5;2;0;zer;zer;0;ux;uz;ux;uz];              %Tran-Body2 to Ground
        # SJDT(:,3)=[1;1;2;0.08*uy;0.02*uy;0.23;zer;zer;zer;zer];  %Dist-Body 1 to 2


        SJDT[:, 1] = Any[4, 1, 0, zer..., 0.1,0.12,0, 0, ux..., uz..., ux..., uz...]#z
        SJDT[:, 2] = Any[1, 1, 2, 0,0.08,0,0,0.00,0, 0.23, zer..., zer..., zer..., zer...]#z
        SJDT[:, 3] = Any[5, 2, 3, si1pr3..., sj2pr3..., 0, ux..., uy..., ux..., uy...]#z
        SJDT[:, 4] = Any[5, 3, 0, zer..., zer..., 0, ux..., uz..., ux..., uz...]#z
        SJDT[:, 5] = Any[1, 2, 3, 0,0.00,0,0,0.00,0, 0.08, zer..., zer..., zer..., zer...]#z
        #SMDT = hcat(json_data["SMDT"]["mass1"], json_data["SMDT"]["mass2"])
         
        SMDT = hcat(vcat(5.0, 0.2,0.2,0.2),
        vcat(5.0, 0.2,0.2,0.2),
        vcat(5.0, 0.2,0.2,0.2),
        )
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        p10=[1., 0., 0., 0.]
        omeg1pr0 = [0, 0, 0]

        r1d0=[0, 0, 0]
        # p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        p1d0=0.5 * mathfunction.EEval(p10)' * 120 * uz
        r2d0=[0,0,9.378]
        r3d0=[0,0,9.378]
        # Initial positions for all bodies (r) and additional parameters (p), with their initial velocities (rd and pd)
        q0 = [0.1,0.12,0, 1., 0., 0., 0.,
        0,0.02,0.1027, 1., 0., 0., 0.,
        0,0,0.1027, 1., 0., 0., 0.,
        #r_bodies3...,1., 0., 0., 0.,
        ]

        qd0 = [r1d0..., p1d0...,
            r2d0..., 0.,0., 0., 0., 
            r3d0..., 0.,0., 0., 0., 
        #    0., 0., 0., 0.,0., 0., 0.
        ]  

        ld_contact = []
        ld_damper = []#damper_g
        p_contact = Any[ld_damper,ld_contact,0,0]

        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0, p_contact

    end
    if app >=  246 && app < 249
        script_dir = @__DIR__

        # Build the path to the JSON file
        json_path = joinpath(script_dir, "config.json")
        json_data = JSON.parsefile(json_path)
        contact_json_path = joinpath(script_dir,"contact", "$contact_json.json")
        contact_json_data = JSON.parsefile(contact_json_path)
        nb = 3      # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 4      # Number of holonomic constraints
        nhc = 16 # Number of holonomic constraint equations
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
        r_bodies1 = [0.1,0.12,0]
        r_bodies2 = [0,0.02,0.1027]
        r_bodies3 = [0,0,0.1027] 
        si1pr3 = [0,0,0.1027]-r_bodies2
        sj2pr3 = [0,0,0.1027]-r_bodies3
        #si1pr2 = [1.,0.,1.]-r_bodies1
        #sj2pr2 = [1.,0.,1.]-r_bodies2

        # SJDT(:,1)=[4;1;0;zer;0.1*ux+0.12*uy;0;ux;uz;ux;uz];   %Cyl-Body1 to Ground
        # SJDT(:,2)=[5;2;0;zer;zer;0;ux;uz;ux;uz];              %Tran-Body2 to Ground
        # SJDT(:,3)=[1;1;2;0.08*uy;0.02*uy;0.23;zer;zer;zer;zer];  %Dist-Body 1 to 2


        SJDT[:, 1] = Any[4, 1, 0, zer..., 0.1,0.12,0, 0, ux..., uz..., ux..., uz...]#z
        SJDT[:, 2] = Any[1, 1, 2, 0,0.08,0,0,0.00,0, 0.23, zer..., zer..., zer..., zer...]#z
        SJDT[:, 3] = Any[5, 2, 3, si1pr3..., sj2pr3..., 0, ux..., uy..., ux..., uy...]#z
        SJDT[:, 4] = Any[5, 3, 0, zer..., zer..., 0, ux..., uz..., ux..., uz...]#z
        #SJDT[:, 5] = Any[1, 2, 3, 0,0.00,0,0,0.00,0, 0.08, zer..., zer..., zer..., zer...]#z
        #SMDT = hcat(json_data["SMDT"]["mass1"], json_data["SMDT"]["mass2"])
         
        SMDT = hcat(vcat(5.0, 0.2,0.2,0.2),
        vcat(5.0, 0.2,0.2,0.2),
        vcat(5.0, 0.2,0.2,0.2),
        )
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        p10=[1., 0., 0., 0.]
        omeg1pr0 = [0, 0, 0]

        r1d0=[0, 0, 0]
        # p1d0 = 0.5 * mathfunction.GEval(p10)' * omeg1pr0
        p1d0=0.5 * mathfunction.EEval(p10)' * 120 * uz
        r2d0=[0,0,9.378]
        r3d0=[0,0,9.378]
        # Initial positions for all bodies (r) and additional parameters (p), with their initial velocities (rd and pd)
        q0 = [0.1,0.12,0, 1., 0., 0., 0.,
        0,0.02,0.1027, 1., 0., 0., 0.,
        0,0,0.1027, 1., 0., 0., 0.,
        #r_bodies3...,1., 0., 0., 0.,
        ]

        qd0 = [r1d0..., p1d0...,
            r2d0..., 0.,0., 0., 0., 
            r3d0..., 0.,0., 0., 0., 
        #    0., 0., 0., 0.,0., 0., 0.
        ]  
        p_eff=contact_json_data["p_eff"]
        contact_onlyy = Dict(
            "type"=>"simple_y",
            "b" => 2,
            "pos" => 0,
            "p" => p_eff,
        )
        println(contact_onlyy)
        ld_contact = [contact_onlyy,]
        ld_damper = []#damper_g
        p_contact = Any[ld_damper,ld_contact,0,0]

        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0, p_contact

    end
end


