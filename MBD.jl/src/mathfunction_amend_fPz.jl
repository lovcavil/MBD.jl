
export bbP2fz,bbP3fz,bbP4fz,bbPhifz,bbPhiqfz
# cid=1030


function bbPhifz(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    E33=[0 0 0;0 0 0;0 0 1]

    d1 = (r1 + A1 * s1pr)
    Phi=(d1[3]*d1[3]-d^2) / 2
    return Phi
end


function bbPhiqfz(i, j, s1pr, s2pr, d,tn, q, par)

    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    E33=[0 0 0;0 0 0;0 0 1]

    d1 = ( r1 + A1 * s1pr)

    Phiq1 = hcat(zeros(1,2),d1[3], zeros(1,4))
    Phiq2 = zeros(1, 7)

    return Phiq1, Phiq2
end

function bbP2fz(i, j, s1pr, s2pr, d,tn, q, qd, par)
    # Assuming par is a structure or a dictionary with the following keys.
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    E33=[0 0 0;0 0 0;0 0 1]
    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1x = BTran(xp1, s1pr)
    a1 = (xr1' + xp1' * BT1')

    P21 = zeros(1, 7)
    P22 = zeros(1, 7)

    d1 = ( r1 + A1 * s1pr)
    a1 = (xr1' + xp1' * BT1')
    a1E33 = a1 * E33

    P21 = hcat((a1E33), (a1E33) * BT1 + d1' *E33* BT1x)
    P22 = zeros(1, 7)

    return P21, P22
end


function bbP3fz(i, j, s1pr, s2pr, d, tn, q, qd, par)
    I3 = Matrix{Float64}(I, 3, 3)  # Equivalent of eye(3) in MATLAB

    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    BT1 = BTran(p1, s1pr)
    BT1d = BTran(p1d, s1pr)
    a1bar = r1d + BT1 * p1d
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    a1barE11=a1bar'*E11

    P31 = zeros(1, 7)
    #P31 = hcat(p1d' * BT1d'*E11, 2 * (a1barE11) * BT1d + p1d'* BT1d' *E11 * BT1)
    #P32 = zeros(1, 7)
    println("nnnnoooo3")
    return P31, P32
end


function bbP4fz(i, j, s1pr, s2pr, d, tn, q, etak, par)
    I3 = Matrix{Float64}(I, 3, 3)

    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    E22=[0 0 0;0 1 0;0 0 0]
    E33=[0 0 0;0 0 0;0 0 1]
    if j == 0
        d12 =  r1 + A1 * s1pr
        P411 = etak * vcat(hcat(I3*E33,E33*BT1), hcat(BT1'*E33,BT1'*E33 * BT1 + KEval(s1pr, E33*d12)))
        P412 = zeros(7, 7)
        P422 = zeros(7, 7)
    end
    return P411, P412, P422
end

