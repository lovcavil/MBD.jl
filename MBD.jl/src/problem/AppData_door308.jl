using LinearAlgebra
using JSON


function model_door_308_spline()
    script_dir = @__DIR__

    # Build the path to the JSON file
    json_path = joinpath(script_dir, "config.json")
    json_data = JSON.parsefile(json_path)
    nb = 3         # Number of bodies
    ngc = 7 * nb    # Number of generalized coordinates
    nh = 2   +7     # Number of holonomic constraints
    nhc = 10   +7    # Number of holonomic constraint equations
    nc = nhc + nb   # Number of constraint equations
    nv = ngc - nc
    nu = nc   
    NTSDA = 0       # Number of TSDA force elements

    ux = [1; 0; 0]
    uy = [0; 1; 0]
    uz = [0; 0; 1]
    zer = zeros(3)
    
    revrootMA = json_data["revrootMA"]
    revrootLA = json_data["revrootLA"]
    rDoor = json_data["rDoor"]
    rMA = json_data["rMA"]
    rLA = json_data["rLA"]

    # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
    # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
    # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
    # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
    SJDT = Array{Any}(undef, 28, nh)
    si1pr1 = (revrootMA-rDoor)
    sj2pr1 = (revrootMA-rMA)
    SJDT[:, 1] = Any[4, 1, 2, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz..., 0.1, 0.05, 0.0, 0.0, 1, 5]  

    si1pr2 = (revrootLA-rDoor)
    sj2pr2 = (revrootLA-rLA)
    SJDT[:, 2] = Any[4, 1, 3, si1pr2..., sj2pr2..., 0, ux..., uz..., ux..., uz..., 0.1, 0.05, 0.0, 0.0, 6, 5]  

    spl_dic = Dict(
    "MF" => fit_xycurve("GUIDE_M.csv",3),
    "MR" => fit_xycurve("GUIDE_M.csv",3), 
    "LF" => fit_xycurve("GUIDE_L.csv",3),
    "LR" => fit_xycurve("GUIDE_L.csv",3),
    "U" =>  fit_xycurve("GUIDE_U.csv",3)
    )
    
    fixrootrollerMG=json_data["fixrootrollerMG"]
    si1pr3= fixrootrollerMG-rMA
    SJDT[:, 3] = Any[1030, 2, 0, si1pr3..., zer..., 999.4630310344, zer..., zer..., zer..., zer..., 0.0, 0.0, 0.0, 0.0, 0, 0]

    fixrootrollerLG=json_data["fixrootrollerMG"]
    si1pr4= fixrootrollerLG-rLA
    SJDT[:, 4] = Any[1030, 3, 0, si1pr4..., zer..., -22.7978911861, zer..., zer..., zer..., zer..., 0.0, 0.0, 0.0, 0.0, 0, 0]

    fixrootrollerMF=json_data["fixrootrollerMF"]
    si1pr5= fixrootrollerMF-rMA
    SJDT[:, 5] = Any[1070, 2, 0, si1pr5..., zer..., spl_dic["MF"], zer..., zer..., zer..., zer..., 0.0, 0.0, 0.0, 0.0, 0, 0]

    fixrootrollerMR=json_data["fixrootrollerMR"]
    si1pr6= fixrootrollerMR-rMA
    SJDT[:, 6] = Any[1070, 2, 0, si1pr6..., zer..., spl_dic["MR"], zer..., zer..., zer..., zer..., 0.0, 0.0, 0.0, 0.0, 0, 0]

    fixrootrollerLF=json_data["fixrootrollerLF"]
    si1pr7= fixrootrollerLF-rLA
    SJDT[:, 7] = Any[1070, 3, 0, si1pr7..., zer..., spl_dic["LF"], zer..., zer..., zer..., zer..., 0.0, 0.0, 0.0, 0.0, 0, 0]

    fixrootrollerLR=json_data["fixrootrollerLR"]
    si1pr8= fixrootrollerLR-rLA
    SJDT[:, 8] = Any[1070, 3, 0, si1pr8..., zer..., spl_dic["LR"], zer..., zer..., zer..., zer..., 0.0, 0.0, 0.0, 0.0, 0, 0]

    fixrootrollerU=json_data["fixrootrollerU"]
    si1pr9= fixrootrollerU-rDoor
    SJDT[:, 9] = Any[1070, 1, 0, si1pr9..., zer..., spl_dic["U"], zer..., zer..., zer..., zer..., 0.0, 0.0, 0.0, 0.0, 0, 0]

    # # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
    # # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
    SMDT = hcat(json_data["SMDT"]["mass1"], json_data["SMDT"]["mass2"], json_data["SMDT"]["mass3"])
    #SMDT = hcat(vcat(30, 90, 90, 90))
    # STSDAT(12, 1): TSDA Data Table
    STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

    # Initial generalized coordinates

    p10 = [1, 0, 0, 0]
    p20 = [1, 0, 0, 0]
    p30 = [1, 0, 0, 0]

    q0 = [rDoor..., p10...,rMA..., p20...,rLA..., p30...]

    omeg1pr0 = [0, 0, 0]

    r1d0 = [-1000, 0, 0]
    r2d0 = [0, 0, 0]
    r3d0 = [0, 0, 0]
    p1d0 = [0, 0, 0, 0]
    p2d0 = [0, 0, 0, 0]
    p3d0 = [0, 0, 0, 0] 
    qd0 = [r1d0..., p1d0...,r2d0..., p2d0...,r3d0..., p3d0...]

    #println(qd0)
    return AppDataStruct("model_door_308",nb,ngc,nh,nc, nv,nu,   NTSDA,SJDT,SMDT,STSDAT,q0,qd0)
end