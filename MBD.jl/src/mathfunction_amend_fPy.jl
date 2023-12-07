
export bbP2fy,bbP3fy,bbP4fy,bbPhify,bbPhiqfy
# cid=1020


function bbPhify(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    # if j == 0
    #     d12 = s2pr - r1 - A1 * s1pr
    #     Phi = ((d12'*E11* d12) - d^2) / 2
    # else
    #     r2, p2 = qPart(q, j)
    #     A2 = ATran(p2)
    #     d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
    #     Phi = ((d12'*E11* d12)  - d^2) / 2
    # end

    flag1=par[11]["p1"]
    d1 = ( r1 + A1 * s1pr)
    Phi=(d1[2]*d1[2]-d^2) / 2
    #println("Ahh Phi")
    return Phi
end


function bbPhiqfy(i, j, s1pr, s2pr, d,tn, q, par)

    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * E11*hcat(I3, BTran(p1, s1pr))
        Phiq2 = zeros(1, 7)
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * E11* hcat(I3, BTran(p1, s1pr))
        Phiq2 = d12' * E11* hcat(I3, BTran(p2, s2pr))
    end
    #Phiq1 = r1'*hcat(I3, zeros(3,4))
    flag1=par[11]["p1"]
    d1 = ( r1 + A1 * s1pr)
    flag2=par[11]["p2"]
    Phiq1 = hcat(zeros(1,1),d1[2],zeros(1,1), zeros(1,4))
    #println("bbPhiqfxc",Phiq1)
    #println("Ahh Phiq")
    return Phiq1, Phiq2
end

function bbP2fy(i, j, s1pr, s2pr, d,tn, q, qd, par)
    # Assuming par is a structure or a dictionary with the following keys.
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1x = BTran(xp1, s1pr)
    a1 = (xr1' + xp1' * BT1')
    a1E11 = a1 * E11
    P21 = zeros(1, 7)
    P22 = zeros(1, 7)

    # if j == 0
    #     d12 = s2pr - r1 - A1 * s1pr
    #     P21 = hcat(a1, a1 * BT1 - d12' * BT1x)
    # elseif j >= 1
    #     r2, p2 = qPart(q, j)
    #     xr2, xp2 = qPart(qd, j)
    #     A2 = ATran(p2)
    #     BT2 = BTran(p2, s2pr)
    #     BT2x = BTran(xp2, s2pr)
    #     a2 = (xr2' + xp2' * BT2')
    #     a2E11 = a2 * E11
    #     d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
    #     P21 = hcat((a1E11 - a2E11), (a1E11 - a2E11) * BT1 - d12' *E11* BT1x)
    #     P22 = hcat((a2E11 - a1E11), (a2E11 - a1E11) * BT2 + d12' *E11* BT2x)
    # end
    # # Special
    # #println("xr2=",xr2')
    # d1 = - r1 - A1 * s1pr
    # P21 = hcat(0,xr1[2],zeros(1, 5))

    flag1=par[11]["p1"]
    d1 = ( r1 + A1 * s1pr)
    flag3=par[11]["p3"]
    a1 = (xr1' + xp1' * BT1')
    a1E22 = a1 * E22
    flag4=par[11]["p4"]
    #P21 = hcat((a1E22), (a1E22) * BT1 - d12' *E22* BT1x)
    P21 = hcat((a1E22), (a1E22) * BT1 + d1' *E22* BT1x)
    P22 = zeros(1, 7)
    #println("P21=",P21)
    #println("-------------+++++++++++++",flag1,flag3,flag4)
    return P21, P22
end


function bbP3fy(i, j, s1pr, s2pr, d, tn, q, qd, par)
    I3 = Matrix{Float64}(I, 3, 3)  # Equivalent of eye(3) in MATLAB

    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    BT1 = BTran(p1, s1pr)
    BT1d = BTran(p1d, s1pr)
    a1bar = r1d + BT1 * p1d
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    a1barE11=a1bar'*E11
    # if j == 0
    #     P31 = hcat(p1d' * BT1d'*E11, 2 * (a1barE11) * BT1d + p1d'* BT1d' *E11 * BT1)
    #     P32 = zeros(1, 7)
    # end

    # if j >= 1
    #     r2, p2 = qPart(q, j)
    #     r2d, p2d = qPart(qd, j)
    #     BT2 = BTran(p2, s2pr)
    #     BT2d = BTran(p2d, s2pr)
    #     a2bar = r2d + BT2 * p2d
    #     a2barE11=a2bar*E11
    #     P31 = hcat(p1d' * BT1d'*E11 - p2d' * BT2d'*E11, 2 * (a1barE11 - a2barE11)' * BT1d +
    #            (p1d' * BT1d' - p2d' * BT2d') *E11* BT1)
    #     P32 = hcat(p2d' * BT2d'*E11 - p1d' * BT1d'*E11, 2 * (a2barE11 - a1barE11)' * BT2d +
    #            (p2d'* BT2d' - p1d' * BT1d') *E11 * BT2)
    # end
    # speical
    P31 = zeros(1, 7)
    #P31 = hcat(p1d' * BT1d'*E11, 2 * (a1barE11) * BT1d + p1d'* BT1d' *E11 * BT1)
    #P32 = zeros(1, 7)
    println("nnnnoooo3")
    return P31, P32
end


function bbP4fy(i, j, s1pr, s2pr, d, tn, q, etak, par)
    I3 = Matrix{Float64}(I, 3, 3)

    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3*E11,E11*BT1), hcat(BT1'*E11,BT1'*E11 * BT1 - KEval(s1pr, E11*d12)))
        P412 = zeros(7, 7)
        P422 = zeros(7, 7)
        d12 =  r1 + A1 * s1pr
        P411 = etak * vcat(hcat(I3*E22,E22*BT1), hcat(BT1'*E22,BT1'*E22 * BT1 + KEval(s1pr, E22*d12)))
        P412 = zeros(7, 7)
        P422 = zeros(7, 7)
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        BT2 = BTran(p2, s2pr)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3*E11, E11*BT1), hcat(BT1'*E11, BT1' *E11* BT1 - KEval(s1pr, E11*d12)))
        P412 = -etak * vcat(hcat(I3*E11, E11*BT2), hcat(BT1'*E11, BT1' *E11* BT2))
        P422 = etak * vcat(hcat(I3*E11, E11*BT2), hcat(BT2'*E11, BT2' *E11* BT2 + KEval(s2pr, E11*d12))    )
    end
    # P411 = etak * vcat(hcat(I3*E22,zeros(3,4)), hcat(zeros(4,3),- KEval(s1pr,  hcat(0,r1[2],zeros(1,1)) )))
    # P412 = zeros(7, 7)
    # P422 = zeros(7, 7)
    #println("nnnnoooo4")
    return P411, P412, P422
end

