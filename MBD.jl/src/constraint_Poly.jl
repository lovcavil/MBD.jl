
export bbP2_Poly,bbP3_Poly,bbP4_Poly,bbPhi_Poly,bbPhiq_Poly
# cid=1060


function bbPhi_Poly(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)

    d1 = ( r1 + A1 * s1pr)
    A=2
    B=-1
    C=0
    d=(A * d1[1] + B * d1[2] + C) / sqrt(A^2 + B^2)

    Phi=d

    return Phi
end


function bbPhiq_Poly(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    E1=[1 0 0]
    E2=[0 1 0]
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    d1 = ( r1 + A1 * s1pr)
    d11 = d1[1]
    d12 = d1[2]
    d1q=hcat(I3,BTran(p1,s1pr))
    # Phiq1 =2*d1'*E11*d1q +[0 -1 0 0 0 0 0];
    Phiq1 = 2*(d11)^1*E1*d1q-E2*d1q
    Phiq2 = zeros(1, 7)
    println("phiq",Phiq1)
    return Phiq1, Phiq2
end

function bbP2_Poly(i, j, s1pr, s2pr, d,tn, q, qd, par)
    # Assuming par is a structure or a dictionary with the following keys.
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    E11=[1 0 0;0 0 0;0 0 0]
    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1x = BTran(xp1, s1pr)

    d1 = ( r1 + A1 * s1pr)
    d1q=hcat(I,BTran(p1,s1pr))
    a1 = (xr1' + xp1' * BT1')

    # P21 = hcat((a1E11), (a1E11) * BT1 + d1' *E11* BT1x)
    k=2
    i=1
    eiT, eiq = ei(i,k,d1,d1q)
    P21 = a1*eiq+hcat([0 0 0],eiT*BT1x)
    P22 = zeros(1, 7)
    println("P2_fP")
    return P21, P22
end

function ei(i,k,d1,d1q)
    Ei=[0 0 0]
    if i==1
        Ei=[1 0 0]
    end
    if i==2
        Ei=[2 0 0]
    end
    if i==3
        Ei=[3 0 0]
    end
    eiT=k*(Ei*d1)[1, 1]^(k-1)*Ei
    eiq=Ei'*k*(k-1)*(Ei*d1)[1, 1]^(k-2)*Ei*d1q
    return eiT, eiq
end

function bbP3_Poly(i, j, s1pr, s2pr, d, tn, q, qd, par)
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

function bbP4_Poly(i, j, s1pr, s2pr, d, tn, q, etak, par)
    I3 = Matrix{Float64}(I, 3, 3)

    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)

    E22=[0 0 0;0 1 0;0 0 0]

    d12 =  r1 + A1 * s1pr;
    P411 = etak * vcat(hcat(I3*E22,E22*BT1), hcat(BT1'*E22,BT1'*E22 * BT1 + KEval(s1pr, E22*d12)))

    P412 = zeros(7, 7);
    P422 = zeros(7, 7);

    return P411, P412, P422
end

