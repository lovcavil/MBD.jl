
export bbP2fy,bbP3fy,bbP4fy,bbPhify,bbPhiqfy
# cid=1020


function bbPhify(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]

    d1 = ( r1 + A1 * s1pr)
    Phi=(d1[2]*d1[2]-d^2) / 2

    return Phi
end


function bbPhiqfy(i, j, s1pr, s2pr, d,tn, q, par)

    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)

    E22=[0 0 0;0 1 0;0 0 0]

    d1 = ( r1 + A1 * s1pr)
    Phiq1 =d1'*E22* hcat(I3,BTran(p1,s1pr));
    Phiq2 = zeros(1, 7)

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


    d1 = ( r1 + A1 * s1pr)

    a1 = (xr1' + xp1' * BT1')
    a1E22 = a1 * E22

    P21 = hcat((a1E22), (a1E22) * BT1 + d1' *E22* BT1x)
    P22 = zeros(1, 7)

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

    E22=[0 0 0;0 1 0;0 0 0]

    d12 =  r1 + A1 * s1pr;
    P411 = etak * vcat(hcat(I3*E22,E22*BT1), hcat(BT1'*E22,BT1'*E22 * BT1 + KEval(s1pr, E22*d12)))

    P412 = zeros(7, 7);
    P422 = zeros(7, 7);

    return P411, P412, P422
end

