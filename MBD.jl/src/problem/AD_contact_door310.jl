using LinearAlgebra
using JSON
using DelimitedFiles
using Dierckx

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

function fit_xycurve(csvfile="", degree=3)
    # Define the function to get the directory of the file
    file_dir = @__DIR__ # Implement this function as needed

    # Load data from CSV
    full_path = joinpath(file_dir, csvfile)
    arr = readdlm(full_path, ',', skipstart=1)
    x = arr[:, 1]
    y = arr[:, 2]
    #z = arr[:, 3]
    sort_indices = sortperm(x)
    x_sorted = x[sort_indices]
    y_sorted = y[sort_indices]
    # Fit the spline
    spline = Spline1D(x_sorted, y_sorted, k=degree)
    return spline
end


function AD310(app)
    if app == 310
        script_dir = @__DIR__

        # Build the path to the JSON file
        json_path = joinpath(script_dir, "config.json")
        json_data = JSON.parsefile(json_path)
        nb = 3 + 7      # Number of bodies
        ngc = 7 * nb    # Number of generalized coordinates
        nh = 2 + 2 + 2 +2 +1 +3       # Number of holonomic constraints
        nhc = 12 + 12 + 12 +12 +6 +3  # Number of holonomic constraint equations
        nc = nhc + nb   # Number of constraint equations

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
        rMG = json_data["fixrootrollerMG"]
        rLG = json_data["fixrootrollerLG"]
        rMF = json_data["fixrootrollerMF"]
        rMR = json_data["fixrootrollerMR"]
        rLF = json_data["fixrootrollerLF"]
        rLR = json_data["fixrootrollerLR"]
        rU = json_data["fixrootrollerU"]

        # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
        # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
        # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
        # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
        SJDT = Array{Any}(undef, 22, nh)
        si1pr1 = (revrootMA - rDoor)
        sj2pr1 = (revrootMA - rMA)
        SJDT[:, 1] = Any[20, 1, 2, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz...]

        si1pr2 = (revrootLA - rDoor)
        sj2pr2 = (revrootLA - rLA)
        SJDT[:, 2] = Any[20, 1, 3, si1pr2..., sj2pr2..., 0, ux..., uz..., ux..., uz...]

        spl_dic = Dict(
            "MF" => fit_xycurve("GUIDE_M.csv", 3),
            "MR" => fit_xycurve("GUIDE_M.csv", 3),
            "LF" => fit_xycurve("GUIDE_L.csv", 3),
            "LR" => fit_xycurve("GUIDE_L.csv", 3),
            "U" => fit_xycurve("GUIDE_U.csv", 3)
        )

        revrootrollerMG = json_data["fixrootrollerMG"]
        si1pr3 = (revrootrollerMG - rMA)
        sj2pr3 = (revrootrollerMG - rMG)
        SJDT[:, 3] = Any[20, 2, 4, si1pr3..., sj2pr3..., 0, ux..., uz..., ux..., uz...]

        revrootrollerLG = json_data["fixrootrollerLG"]
        si1pr4 = (revrootrollerLG - rLA)
        sj2pr4 = (revrootrollerLG - rLG)
        SJDT[:, 4] = Any[20, 3, 5, si1pr4..., sj2pr4..., 0, ux..., uz..., ux..., uz...]

        revrootrollerMF = json_data["fixrootrollerMF"]
        si1pr5 = (revrootrollerMF - rMA)
        sj2pr5 = (revrootrollerMF - rMF)
        SJDT[:, 5] = Any[20, 2, 6, si1pr5..., sj2pr5..., 0, ux..., uz..., ux..., uz...]

        revrootrollerMR = json_data["fixrootrollerMR"]
        si1pr6 = (revrootrollerMR - rMA)
        sj2pr6 = (revrootrollerMR - rMR)
        SJDT[:, 6] = Any[20, 2, 7, si1pr6..., sj2pr6..., 0, ux..., uz..., ux..., uz...]

        revrootrollerLF = json_data["fixrootrollerLF"]
        si1pr7 = (revrootrollerLF - rLA)
        sj2pr7 = (revrootrollerLF - rLF)
        SJDT[:, 7] = Any[20, 3, 8, si1pr7..., sj2pr7..., 0, ux..., uz..., ux..., uz...]

        revrootrollerLR = json_data["fixrootrollerLR"]
        si1pr8 = (revrootrollerLR - rLA)
        sj2pr8 = (revrootrollerLR - rLR)
        SJDT[:, 8] = Any[20, 3, 9, si1pr8..., sj2pr8..., 0, ux..., uz..., ux..., uz...]

        revrootrollerU = json_data["fixrootrollerU"]
        si1pr9 = (revrootrollerU - rDoor)
        sj2pr9 = (revrootrollerU - rU)
        SJDT[:, 9] = Any[20, 1, 10, si1pr9..., sj2pr9..., 0, ux..., uz..., ux..., uz...]

        fixrootrollerMF = json_data["fixrootrollerMF"]
        si1pr5 = fixrootrollerMF - rMA
        SJDT[:, 10] = Any[1020, 6, 0, zer..., zer..., -792.8698947168, zer..., zer..., zer..., zer...]

        # fixrootrollerMR = [3641.0808112684, -790.5217978247, 1015.9232021959]
        # si1pr6 = fixrootrollerMR - rMA
        # SJDT[:, 11] = Any[1020, 7, 0, zer..., zer..., -790.5217978247, zer..., zer..., zer..., zer...]

        fixrootrollerLF = [2451.9295055222, -808.8362561184, 1.6216974129]
        si1pr7 = fixrootrollerLF - rLA
        SJDT[:, 11] = Any[1020, 8, 0, zer..., zer..., -808.8362561184, zer..., zer..., zer..., zer...]

        # fixrootrollerLR = [2503.9248514667; -808.8404593307; 2.3161700361]
        # si1pr8 = fixrootrollerLR - rLA
        # SJDT[:, 13] = Any[1020, 9, 0, zer..., zer..., -808.8404593307, zer..., zer..., zer..., zer...]

        fixrootrollerU = [2608.1721925959; -672.3724048378; 1385.4640400786]
        si1pr9 = fixrootrollerU - rDoor
        SJDT[:, 12] = Any[1020, 10, 0, zer..., zer..., -672.3724048378, zer..., zer..., zer..., zer...]

        # # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
        # # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
        SMDT = hcat(json_data["SMDT"]["mass1"], json_data["SMDT"]["mass2"], json_data["SMDT"]["mass3"],
         json_data["SMDT"]["mass4"], json_data["SMDT"]["mass5"],
         json_data["SMDT"]["mass6"], json_data["SMDT"]["mass7"],json_data["SMDT"]["mass8"], json_data["SMDT"]["mass9"],
         json_data["SMDT"]["mass10"]
         )
        #SMDT = hcat(vcat(30, 90, 90, 90))
        # STSDAT(12, 1): TSDA Data Table
        STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

        # Initial positions for all bodies (r) and additional parameters (p), with their initial velocities (rd and pd)
        r_bodies = Any[rDoor, rMA, rLA, rMG, rLG, rMF, rMR, rLF, rLR, rU]
        p_initial = [1., 0., 0., 0.]  # Common initial 'p' for all bodies
        pd_initial = [0., 0., 0., 0.]  # Common initial 'pd' (velocity) for all 'p'

        # Initial velocities for 'r' bodies - first body has a unique velocity
        rd_initials = [0., 0., 0.]
        # Constructing q0 and qd0 using loops
        q0 = Float64[]
        qd0 = Float64[]

        for i in 1:length(r_bodies)
            append!(q0, r_bodies[i]...)
            append!(q0, p_initial...)
            append!(qd0, rd_initials...)
            append!(qd0, pd_initial...)
        end

        contact_mg = Dict(
            "b" => 4,
            "pos" => 999.4630310344,
        )
        contact_lg = Dict(
            "b" => 5,
            "pos" => -22.7978911861,
        )
        contact_g = Dict(
            "b" => 1,
            "pos" => 634.2367244142,
        )
        contact_g1 = Dict(
            "b" => 1,
            "pos" => -634.2367244142,
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
        damper_g = Dict(
            "b" => 1,
            "damp" => 10000,
        )
        ld_damper = [damper_lg,  damper_mg]#damper_g
        p_contact = Any[ld_damper, ld_contact]



        #println(qd0)

        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0, p_contact

    end
end