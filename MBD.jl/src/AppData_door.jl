using LinearAlgebra

export AppData

function model_door_2()
    nb = 3         # Number of bodies
    ngc = 7 * nb    # Number of generalized coordinates
    nh = 2   +0     # Number of holonomic constraints
    nhc = 10   +0    # Number of holonomic constraint equations
    nc = nhc + nb   # Number of constraint equations
    NTSDA = 0       # Number of TSDA force elements

    ux = [1; 0; 0]
    uy = [0; 1; 0]
    uz = [0; 0; 1]
    zer = zeros(3)
    
    revrootMA=[2735.561844,	-754.2999088,	996.9649782]
    revrootLA=[1618.050638,	-687.9641442,	-23.00674629]

    rDoor = [2137.170608,    -849.1241534,	620.3160646]
    rMA =   [2794.903126,	-747.0129068,	983.1187542]
    rLA =    [1628.65274,	-671.9118677,	-17.63839392]
    r00_ = [0  , 0  , 0]
    # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
    # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
    # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
    # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
    SJDT = Array{Any}(undef, 22, nh)
    si1pr1 = (revrootMA-rDoor)
    sj2pr1 = (revrootMA-rMA)
    SJDT[:, 1] = Any[4, 1, 2, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz...]  

    si1pr2 = (revrootLA-rDoor)
    sj2pr2 = (revrootLA-rLA)
    SJDT[:, 2] = Any[4, 1, 3, si1pr2..., sj2pr2..., 0, ux..., uz..., ux..., uz...]  

    # fixrootrollerMG=[2853.523358,	-733.4159466,	986.1766265]
    # si1pr3= fixrootrollerMG-rMA
    # SJDT[:, 3] = Any[1030, 2, 0, si1pr3..., zer..., 986.1766265, zer..., zer..., zer..., zer...]

    # fixrootrollerLG=[1636.813753,	-651.1483516,	-23.00692676]
    # si1pr4= fixrootrollerLG-rLA
    # SJDT[:, 4] = Any[1030, 3, 0, si1pr4..., zer..., -23.00692676, zer..., zer..., zer..., zer...]

    # fixrootrollerMF=[2842.601738,	-713.5516563,	1001.72917    ]
    # si1pr5= fixrootrollerMF-rMA
    # SJDT[:, 5] = Any[1020, 2, 0, si1pr5..., zer..., -713.5516563, zer..., zer..., zer..., zer...]

    # # fixrootrollerMR=[2866.110593,	-751.9204626,	1002.139357    ]
    # # si1pr6= fixrootrollerMR-rMA
    # # SJDT[:, 6] = Any[1020, 2, 0, si1pr6..., zer..., -751.9204626, zer..., zer..., zer..., zer...]

    # fixrootrollerLF=[1613.653846,	-639.3795907,	1.772774656    ]
    # si1pr7= fixrootrollerLF-rLA
    # SJDT[:, 6] = Any[1020, 3, 0, si1pr7..., zer..., -639.3795907, zer..., zer..., zer..., zer...]

    # # fixrootrollerLR=[1659.986185,	-662.9870967,	1.772774656    ]
    # # si1pr8= fixrootrollerLR-rLA
    # # SJDT[:, 8] = Any[1020, 3, 0, si1pr8..., zer..., -662.9870967, zer..., zer..., zer..., zer...]

    # fixrootrollerU=[1762.934676,	-556.0009449,	1383.964647    ]
    # si1pr9= fixrootrollerU-rDoor
    # SJDT[:, 7] = Any[1020, 1, 0, si1pr9..., zer..., -556.0009449, zer..., zer..., zer..., zer...]

    # # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
    # # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
    SMDT = hcat(vcat(52.5, 4.37, 3.3, 1.09),
    vcat(0.5106279507, 1090.7366643991/1000000, 1018.4922811735/1000000, 273.6987895571/1000000),
    vcat(0.1767413708, 119.3068245091/1000000, 86.021251567/1000000, 54.5063908776/1000000))
    #SMDT = hcat(vcat(30, 90, 90, 90))
    # STSDAT(12, 1): TSDA Data Table
    STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

    # Initial generalized coordinates

    p10 = [1, 0, 0, 0]
    p20 = [1, 0, 0, 0]
    p30 = [1, 0, 0, 0]

    q0 = [rDoor..., p10...,rMA..., p20...,rLA..., p30...]

    omeg1pr0 = [0, 0, 0]

    r1d0 = [0, 0, 0]
    r2d0 = [0, 0, 0]
    r3d0 = [0, 0, 0]
    p1d0 = [0, 0, 0, 0]
    p2d0 = [0, 0, 0, 0]
    p3d0 = [0, 0, 0, 0] 
    qd0 = [r1d0..., p1d0...,r2d0..., p2d0...,r3d0..., p3d0...]

    #println(qd0)
    return AppDataStruct("model_door",nb,ngc,nh,nc,NTSDA,SJDT,SMDT,STSDAT,q0,qd0)
end


function model_door_3()
    nb = 3         # Number of bodies
    ngc = 7 * nb    # Number of generalized coordinates
    nh = 2   +7     # Number of holonomic constraints
    nhc = 10   +7    # Number of holonomic constraint equations
    nc = nhc + nb   # Number of constraint equations
    NTSDA = 0       # Number of TSDA force elements

    ux = [1; 0; 0]
    uy = [0; 1; 0]
    uz = [0; 0; 1]
    zer = zeros(3)
    
    revrootMA=[3580.8956927591, -906.558066787, 1012.9012888899]
    revrootLA=[2478.2110818085, -850.9747241734, -21.3449744949]

    rDoor = [2985.2252251857, -1000.9736541954, 634.2367244142]
    rMA =   [3603.1181651241, -851.633617828, 999.1776716127]
    rLA =    [2480.3079789101, -831.6969825016, -16.6079173541]

    # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
    # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
    # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
    # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
    SJDT = Array{Any}(undef, 22, nh)
    si1pr1 = (revrootMA-rDoor)
    sj2pr1 = (revrootMA-rMA)
    SJDT[:, 1] = Any[4, 1, 2, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz...]  

    si1pr2 = (revrootLA-rDoor)
    sj2pr2 = (revrootLA-rLA)
    SJDT[:, 2] = Any[4, 1, 3, si1pr2..., sj2pr2..., 0, ux..., uz..., ux..., uz...]  

    fixrootrollerMG=[3618.9570586645, -789.5653844227, 999.4630310344]
    si1pr3= fixrootrollerMG-rMA
    SJDT[:, 3] = Any[1030, 2, 0, si1pr3..., zer..., 999.4630310344, zer..., zer..., zer..., zer...]

    fixrootrollerLG=[2478.2407670739, -809.6981769731, -22.7978911861]
    si1pr4= fixrootrollerLG-rLA
    SJDT[:, 4] = Any[1030, 3, 0, si1pr4..., zer..., -22.7978911861, zer..., zer..., zer..., zer...]

    fixrootrollerMF=[3596.1590635603, -792.8698947168, 1014.6516952095]
    si1pr5= fixrootrollerMF-rMA
    SJDT[:, 5] = Any[1020, 2, 0, si1pr5..., zer..., -792.8698947168, zer..., zer..., zer..., zer...]

    fixrootrollerMR=[3641.0808112684, -790.5217978247, 1015.9232021959]
    si1pr6= fixrootrollerMR-rMA
    SJDT[:, 6] = Any[1020, 2, 0, si1pr6..., zer..., -790.5217978247, zer..., zer..., zer..., zer...]

    fixrootrollerLF=[2451.9295055222, -808.8362561184, 1.6216974129]
    si1pr7= fixrootrollerLF-rLA
    SJDT[:, 7] = Any[1020, 3, 0, si1pr7..., zer..., -808.8362561184, zer..., zer..., zer..., zer...]

    fixrootrollerLR=[2503.9248514667; -808.8404593307; 2.3161700361]
    si1pr8= fixrootrollerLR-rLA
    SJDT[:, 8] = Any[1020, 3, 0, si1pr8..., zer..., -808.8404593307, zer..., zer..., zer..., zer...]

    fixrootrollerU=[2608.1721925959; -672.3724048378; 1385.4640400786]
    si1pr9= fixrootrollerU-rDoor
    SJDT[:, 9] = Any[1020, 1, 0, si1pr9..., zer..., -672.3724048378, zer..., zer..., zer..., zer...]

    # # SMDT(4, nb): Mass Data Table (With diagonal inertia matrix)
    # # SMDT = [[m1, J11, J12, J13], ..., [mnb, Jnb1, Jnb2, Jnb3]]
    SMDT = hcat(vcat(52.5, 4.37, 3.3, 1.09),
    vcat(0.5106279507, 1090.7366643991/1000000, 1018.4922811735/1000000, 273.6987895571/1000000),
    vcat(0.1767413708, 119.3068245091/1000000, 86.021251567/1000000, 54.5063908776/1000000))
    #SMDT = hcat(vcat(30, 90, 90, 90))
    # STSDAT(12, 1): TSDA Data Table
    STSDAT = NTSDA == 0 ? zeros(12, NTSDA) : []  # Initialize if NTSDA == 0

    # Initial generalized coordinates

    p10 = [1, 0, 0, 0]
    p20 = [1, 0, 0, 0]
    p30 = [1, 0, 0, 0]

    q0 = [rDoor..., p10...,rMA..., p20...,rLA..., p30...]

    omeg1pr0 = [0, 0, 0]

    r1d0 = [10, 0, 0]
    r2d0 = [0, 0, 0]
    r3d0 = [0, 0, 0]
    p1d0 = [0, 0, 0, 0]
    p2d0 = [0, 0, 0, 0]
    p3d0 = [0, 0, 0, 0] 
    qd0 = [r1d0..., p1d0...,r2d0..., p2d0...,r3d0..., p3d0...]

    #println(qd0)
    return AppDataStruct("model_door",nb,ngc,nh,nc,NTSDA,SJDT,SMDT,STSDAT,q0,qd0)
end























function model_door()
    nb = 3         # Number of bodies
    ngc = 7 * nb    # Number of generalized coordinates
    nh = 2   +5     # Number of holonomic constraints
    nhc = 10   +5    # Number of holonomic constraint equations
    nc = nhc + nb   # Number of constraint equations
    NTSDA = 0       # Number of TSDA force elements

    ux = [1; 0; 0]
    uy = [0; 1; 0]
    uz = [0; 0; 1]
    zer = zeros(3)
    
    revrootMA=[2735.561844,	-754.2999088,	996.9649782]
    revrootLA=[1618.050638,	-687.9641442,	-23.00674629]

    rDoor = [2137.170608,    -849.1241534,	620.3160646]
    rMA =   [2794.903126,	-747.0129068,	983.1187542]
    rLA =    [1628.65274,	-671.9118677,	-17.63839392]
    r00_ = [0  , 0  , 0]
    # SJDT[:, k] = [t, i, j, sipr, sjpr, d, uxipr, uzipr, uxjpr, uzjpr]
    # k = joint number; t = joint type (1=Dist, 2=Sph, 3=Cyl, 4=Rev, 5=Tran, 
    # 6=Univ, 7=Strut, 8=Rev-Sph); i & j = bodies connected, i > 0; 1010=fxc
    # si & sjpr = vectors to Pi & Pj; d = dist.; uxipr, uzipr, uxjpr, uzjpr
    SJDT = Array{Any}(undef, 22, nh)
    si1pr1 = (revrootMA-rDoor)
    sj2pr1 = (revrootMA-rMA)
    SJDT[:, 1] = Any[4, 1, 2, si1pr1..., sj2pr1..., 0, ux..., uz..., ux..., uz...]  

    si1pr2 = (revrootLA-rDoor)
    sj2pr2 = (revrootLA-rLA)
    SJDT[:, 2] = Any[4, 1, 3, si1pr2..., sj2pr2..., 0, ux..., uz..., ux..., uz...]  

    fixrootrollerMG=[2853.523358,	-733.4159466,	986.1766265]
    si1pr3= fixrootrollerMG-rMA
    SJDT[:, 3] = Any[1030, 2, 0, si1pr3..., zer..., 986.1766265, zer..., zer..., zer..., zer...]

    fixrootrollerLG=[1636.813753,	-651.1483516,	-23.00692676]
    si1pr4= fixrootrollerLG-rLA
    SJDT[:, 4] = Any[1030, 3, 0, si1pr4..., zer..., -23.00692676, zer..., zer..., zer..., zer...]

    fixrootrollerMF=[2842.601738,	-713.5516563,	1001.72917    ]
    si1pr5= fixrootrollerMF-rMA
    SJDT[:, 5] = Any[1020, 2, 0, si1pr5..., zer..., -713.5516563, zer..., zer..., zer..., zer...]

    # fixrootrollerMR=[2866.110593,	-751.9204626,	1002.139357    ]
    # si1pr6= fixrootrollerMR-rMA
    # SJDT[:, 6] = Any[1020, 2, 0, si1pr6..., zer..., -751.9204626, zer..., zer..., zer..., zer...]

    fixrootrollerLF=[1613.653846,	-639.3795907,	1.772774656    ]
    si1pr7= fixrootrollerLF-rLA
    SJDT[:, 6] = Any[1020, 3, 0, si1pr7..., zer..., -639.3795907, zer..., zer..., zer..., zer...]

    # fixrootrollerLR=[1659.986185,	-662.9870967,	1.772774656    ]
    # si1pr8= fixrootrollerLR-rLA
    # SJDT[:, 8] = Any[1020, 3, 0, si1pr8..., zer..., -662.9870967, zer..., zer..., zer..., zer...]

    fixrootrollerU=[1762.934676,	-556.0009449,	1383.964647    ]
    si1pr9= fixrootrollerU-rDoor
    SJDT[:, 7] = Any[1020, 1, 0, si1pr9..., zer..., -556.0009449, zer..., zer..., zer..., zer...]

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

    q0 = [rDoor..., p10...,rMA..., p20...,rLA..., p30...]

    omeg1pr0 = [0, 0, 0]

    r1d0 = [0, 0, 0]
    r2d0 = [0, 0, 0]
    r3d0 = [0, 0, 0]
    p1d0 = [0, 0, 0, 0]
    p2d0 = [0, 0, 0, 0]
    p3d0 = [0, 0, 0, 0] 
    qd0 = [r1d0..., p1d0...,r2d0..., p2d0...,r3d0..., p3d0...]

    #println(qd0)
    return AppDataStruct("model_door",nb,ngc,nh,nc,NTSDA,SJDT,SMDT,STSDAT,q0,qd0)
end