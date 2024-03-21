using LinearAlgebra
include("../mathfunction_II7.jl")
function AD(app)
    if app ==205 ||app==209||app==210||app==211||app==212  # single Pendulum+plainer x=1  Spherical to Ground
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
        p_contact=0
        if app==205
            p_contact=Any[0,0]
        end
        if app==209
            p_contact=Any[1,1000]
            contact_mg=Dict(
                "b" => 1,
                "pos" => -0.5, 
                )
         
            ld_contact=[contact_mg,]
            damper=Dict(
                "b" => 1,
                "damp" => 1, 
                )
         
            ld_damper=[damper,]
            p_contact=Any[ld_damper,ld_contact]
        end
        if app==210
            contact_mg=Dict(
                "b" => 1,
                "pos" => -0.5, 
                )
            ld_contact=[contact_mg,]
            damper=Dict(
                "b" => 1,
                "damp" => 100, 
                )
         
            ld_damper=[damper,]
            p_contact=Any[ld_damper,ld_contact]
        end
        if app==211
            contact_mg=Dict(
                "b" => 1,
                "pos" => -0.5, 
                )
            ld_contact=[contact_mg,]
            damper=Dict(
                "b" => 1,
                "damp" => 500, 
                )
         
            ld_damper=[damper,]
            p_contact=Any[ld_damper,ld_contact]
        end
        if app==212
            contact_mg=Dict(
                "b" => 1,
                "pos" => -0.5, 
                )
            ld_contact=[contact_mg,]
            damper=Dict(
                "b" => 1,
                "damp" => 1000, 
                )
         
            ld_damper=[damper,]
            p_contact=Any[ld_damper,ld_contact]
        end
        return nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q0, qd0,p_contact
    end    
end

