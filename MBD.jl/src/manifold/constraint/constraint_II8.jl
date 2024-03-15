

function bbP2dist(i, j, s1pr, s2pr, d,tn, q, qd, par)
    # Assuming par is a structure or a dictionary with the following keys.
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1x = BTran(xp1, s1pr)
    a1 = (xr1' + xp1' * BT1')

    P21 = zeros(1, 7)
    P22 = zeros(1, 7)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P21 = [a1, a1 * BT1 - d12' * BT1x]
    elseif j >= 1
        r2, p2 = qPart(q, j)
        xr2, xp2 = qPart(qd, j)
        A2 = ATran(p2)
        BT2 = BTran(p2, s2pr)
        BT2x = BTran(xp2, s2pr)
        a2 = (xr2' + xp2' * BT2')
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P21 = hcat((a1 - a2), (a1 - a2) * BT1 - d12' * BT1x)
        P22 = hcat((a2 - a1), (a2 - a1) * BT2 + d12' * BT2x)
    end

    return P21, P22
end

function bbP2dot1(i, j, a1pr, a2pr,tn, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    
    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, a1pr)
    BTx1 = BTran(xp1, a1pr)
    
    P21 = zeros(1, 7)
    P22 = zeros(1, 7)
    
    if j == 0
        P21 = hcat(zeros(1, 3), a2pr' * BTx1)
    elseif j >= 1
        r2, p2 = qPart(q, j)
        xr2, xp2 = qPart(qd, j)
        A2 = ATran(p2)
        BT2 = BTran(p2, a2pr)
        BTx2 = BTran(xp2, a2pr)
        P21 = hcat(zeros(1, 3), a2pr' * A2' * BTx1 + xp2' * BT2' * BT1)
        P22 = hcat(zeros(1, 3), a1pr' * A1' * BTx2 + xp1' * BT1' * BT2)
    end
    
    return P21, P22
end


function bbP2dot2(i, j, a2pr, s1pr, s2pr,tn, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    BT1s1 = BTran(p1, s1pr)
    xr1, xp1 = qPart(qd, i)
    BTx1s1 = BTran(xp1, s1pr)
    
    P21 = zeros(1, 7)
    P22 = zeros(1, 7)
    
    if j == 0
        P21 = -hcat(zeros(1, 3), a2pr' * BTx1s1)
    elseif j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        xr2, xp2 = qPart(qd, j)
        BT2s2 = BTran(p2, s2pr)
        BT2a2 = BTran(p2, a2pr)
        BTx2s2 = BTran(xp2, s2pr)
        BTx2a2 = BTran(xp2, a2pr)
        c = xp2' * BT2a2' * BT2s2 + d12' * BTx2a2 + a2pr' * A2' * BTx2s2 + (xr2' + xp2' * BT2s2') * BT2a2
        P21 = -hcat(xp2' * BT2a2', a2pr' * A2' * BTx1s1 + xp2' * BT2a2' * BT1s1)
        P22 = hcat(xp2' * BT2a2', c - (xr1' + xp1' * BT1s1') * BT2a2)
    end
    
    return P21, P22
end

function bbP2fxc(i, j, s1pr, s2pr, d,tn, q, qd, par)
    # Assuming par is a structure or a dictionary with the following keys.
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    E11=[1 0 0;0 0 0;0 0 0]
    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1x = BTran(xp1, s1pr)
    a1 = (xr1' + xp1' * BT1')
    a1E11 = a1 * E11
    P21 = zeros(1, 7)
    P22 = zeros(1, 7)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P21 = hcat(a1, a1 * BT1 - d12' * BT1x)
    elseif j >= 1
        r2, p2 = qPart(q, j)
        xr2, xp2 = qPart(qd, j)
        A2 = ATran(p2)
        BT2 = BTran(p2, s2pr)
        BT2x = BTran(xp2, s2pr)
        a2 = (xr2' + xp2' * BT2')
        a2E11 = a2 * E11
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P21 = hcat((a1E11 - a2E11), (a1E11 - a2E11) * BT1 - d12' *E11* BT1x)
        P22 = hcat((a2E11 - a1E11), (a2E11 - a1E11) * BT2 + d12' *E11* BT2x)
    end
    # Special
    #println("xr1=",xr1')
    P21 = hcat(xr1[1],zeros(1, 6))
    return P21, P22
end



function bbP2RotDr(i, j, vx1pr, vy1pr, vx2pr, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    s = bbPhidot1(i, j, vy1pr, vx2pr,tn, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr,tn, q, par)

    if abs(c) >= abs(s)
        P21, P22 = bbP2dot1(i, j, vy1pr, vx2pr,tn, q, qd, par)
    else
        P21, P22 = bbP2dot1(i, j, vx1pr, vx2pr,tn, q, qd, par)
    end

    return P21, P22
end

# Example usage:
# You would need to define the variables `i`, `j`, `vx1pr`, `vy1pr`, `vx2pr`, `q`, `qd`, and `par`
# appropriately before calling this function.
# P21, P22 = bbP2RotDr(i, j, vx1pr, vy1pr, vx2pr, q, qd, par)
function bbP2sph(i, j, s1pr, s2pr,tn, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    xr1, xp1 = qPart(qd, i)
    BT1x = BTran(xp1, s1pr)

    P21 = zeros(3, 7)
    P22 = zeros(3, 7)
    
    if j == 0
        P21[:, 4:7] = -BT1x
    elseif j >= 1
        xr2, xp2 = qPart(qd, j)
        BT2x = BTran(xp2, s2pr)
        P21[:, 4:7] = -BT1x
        P22[:, 4:7] = BT2x
    end

    return P21, P22
end

function bbP3dist(i, j, s1pr, s2pr, d, tn, q, qd, par)
    I3 = Matrix{Float64}(I, 3, 3)  # Equivalent of eye(3) in MATLAB

    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    BT1 = BTran(p1, s1pr)
    BT1d = BTran(p1d, s1pr)
    a1bar = r1d + BT1 * p1d
    if j == 0
        P31 = hcat(p1d' * BT1d', 2 * (a1bar)' * BT1d + p1d' * BT1d' * BT1)
        P32 = zeros(1, 7)
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        r2d, p2d = qPart(qd, j)
        BT2 = BTran(p2, s2pr)
        BT2d = BTran(p2d, s2pr)
        a2bar = r2d + BT2 * p2d
        P31 = hcat(p1d' * BT1d' - p2d' * BT2d', 2 * (a1bar - a2bar)' * BT1d +
               (p1d' * BT1d' - p2d' * BT2d') * BT1)
        P32 = hcat(p2d' * BT2d' - p1d' * BT1d', 2 * (a2bar - a1bar)' * BT2d +
               (p2d' * BT2d' - p1d' * BT1d') * BT2)
    end

    return P31, P32
end


function bbP3dot1(i, j, a1pr, a2pr, tn, q, qd, par)
    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    BT1 = BTran(p1, a1pr)
    BT1d = BTran(p1d, a1pr)

    P31 = zeros(1, 7)
    P32 = zeros(1, 7)

    if j >= 1
        r2, p2 = qPart(q, j)
        r2d, p2d = qPart(qd, j)
        BT2 = BTran(p2, a2pr)
        BT2d = BTran(p2d, a2pr)
        P31 = hcat(0, 0, 0, 2 * (p2d'* BT2' * BT1d) + (p2d'* BT2d' * BT1))# from source code error
        P32 = hcat(0, 0, 0, 2 * (p1d'* BT1' * BT2d) + (p1d'* BT1d' * BT2))
    end 

    return P31, P32
end
function bbP3dot2(i, j, a2pr, s1pr, s2pr, tn, q, qd, par)
    # nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)

    I3 = Matrix{Float64}(I, 3, 3)
    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1d = BTran(p1d, s1pr)

    P31 = zeros(1, 7)
    P32 = zeros(1, 7)

    if j >= 1
        r2, p2 = qPart(q, j)
        r2d, p2d = qPart(qd, j)
        A2 = ATran(p2)
        BT2s2 = BTran(p2, s2pr)
        BT2a2 = BTran(p2, a2pr)
        BT2s2d = BTran(p2d, s2pr)
        BT2a2d = BTran(p2d, a2pr)
        d = p2d' * BT2a2d' * BT2s2 + p2d' * BT2s2d' * BT2a2 + 2 * r2d' * BT2a2d +
            2 * p2d' * BT2s2' * BT2a2d + 2 * p2d' * BT2a2' * BT2s2d
        P31 = -hcat(p2d' * BT2a2', 2 * p2d' * BT2a2' * BT1d + p2d' * BT2a2d' * BT1)
        P32 = hcat(p2d' * BT2a2', -p1d' * BT1d' * BT2a2 - 2 * (r1d' + p1d' * BT1') * BT2a2d + d)
    end

    return P31, P32
end

function bbP3fxc(i, j, s1pr, s2pr, d, tn, q, qd, par)
    I3 = Matrix{Float64}(I, 3, 3)  # Equivalent of eye(3) in MATLAB

    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    BT1 = BTran(p1, s1pr)
    BT1d = BTran(p1d, s1pr)
    a1bar = r1d + BT1 * p1d
    E11=[1 0 0;0 0 0;0 0 0]
    a1barE11=a1bar'*E11
    if j == 0
        P31 = hcat(p1d' * BT1d'*E11, 2 * (a1barE11) * BT1d + p1d'* BT1d' *E11 * BT1)
        P32 = zeros(1, 7)
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        r2d, p2d = qPart(qd, j)
        BT2 = BTran(p2, s2pr)
        BT2d = BTran(p2d, s2pr)
        a2bar = r2d + BT2 * p2d
        a2barE11=a2bar*E11
        P31 = hcat(p1d' * BT1d'*E11 - p2d' * BT2d'*E11, 2 * (a1barE11 - a2barE11)' * BT1d +
               (p1d' * BT1d' - p2d' * BT2d') *E11* BT1)
        P32 = hcat(p2d' * BT2d'*E11 - p1d' * BT1d'*E11, 2 * (a2barE11 - a1barE11)' * BT2d +
               (p2d'* BT2d' - p1d' * BT1d') *E11 * BT2)
    end
    # speical
    P31 = zeros(1, 7)
    println("nnnnoooo")
    return P31, P32
end

function bbP4dist(i, j, s1pr, s2pr, d, tn, q, etak, par)
    I3 = Matrix{Float64}(I, 3, 3)

    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3,BT1), hcat(BT1',BT1' * BT1 - KEval(s1pr, d12)))
        P412 = zeros(7, 7)
        P422 = zeros(7, 7)
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        BT2 = BTran(p2, s2pr)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3, BT1), hcat(BT1', BT1' * BT1 - KEval(s1pr, d12)))
        P412 = -etak * vcat(hcat(I3, BT2), hcat(BT1', BT1' * BT2))
        P422 = etak * vcat(hcat(I3, BT2), hcat(BT2', BT2' * BT2 + KEval(s2pr, d12))    )
    end

    return P411, P412, P422
end

function bbP4dot1(i, j, a1pr, a2pr, tn, q, etak, par)
    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, a1pr)
    A1 = ATran(p1)

    P411 = zeros(7, 7)
    P412 = zeros(7, 7)
    P422 = zeros(7, 7)

    if j == 0
        P411 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3), KEval(a1pr, a2pr)))
    elseif j >= 1
        r2, p2 = qPart(q, j)
        BT2 = BTran(p2, a2pr)
        A2 = ATran(p2)

        P411 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3), KEval(a1pr, A2 * a2pr)))
        P412 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3), BT1' * BT2))
        P422 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3), KEval(a2pr, A1 * a1pr)))
    end

    return P411, P412, P422
end
function bbP4dot2(i, j, a2pr, s1pr, s2pr, tn, q, etak, par)
    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)

    P411 = zeros(7, 7)
    P412 = zeros(7, 7)
    P422 = zeros(7, 7)

    if j == 0 #to ch8
        P411 = -etak * vcat(zeros(3, 7), hcat(zeros(4, 3), KEval(s1pr, a2pr)))
    elseif j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        BT2a2 = BTran(p2, a2pr)
        BT2s2 = BTran(p2, s2pr)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        e = KEval(s2pr, A2 * a2pr) + KEval(a2pr, d12) + BT2s2' * BT2a2 + BT2a2' * BT2s2
        P411 = -etak * vcat(zeros(3, 7), hcat(zeros(4, 3),KEval(s1pr, A2 * a2pr)))#to ch8
        P412 = -etak * vcat(hcat(zeros(3, 3), BT2a2), hcat(zeros(4, 3), BT1' * BT2a2))
        P422 = etak * vcat(hcat(zeros(3, 3), BT2a2), hcat(BT2a2', e))
    end

    return P411, P412, P422
end
function bbP4fxc(i, j, s1pr, s2pr, d, tn, q, etak, par)
    I3 = Matrix{Float64}(I, 3, 3)

    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3*E11,E11*BT1), hcat(BT1'*E11,BT1'*E11 * BT1 - KEval(s1pr, E11*d12)))
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
    P411 = etak * vcat(hcat(I3*E11,zeros(3,4)), hcat(zeros(4,3),- KEval(s1pr,  vcat(r1[1],zeros(2,1)) )))
    println("nnnnoooo4")
    return P411, P412, P422
end

function bbP4sph(i, j, s1pr, s2pr, tn, q, etak, par)
    # P412 is identically zero
    P412 = zeros(7, 7)

    if j == 0
        P411 = vcat(zeros(3, 7), hcat(zeros(4, 3),-KEval(s1pr, etak)))
        P422 = zeros(7, 7)
    end

    if j >= 1
        P411 = vcat(zeros(3, 7), hcat(zeros(4, 3),-KEval(s1pr, etak)))
        P422 = vcat(zeros(3, 7), hcat(zeros(4, 3),KEval(s2pr, etak)))
    end

    return P411, P412, P422
end

# Example usage:
# You would need to define the variables `i`, `j`, `s1pr`, `s2pr`, `q`, `qd`, and `par`
# appropriately before calling this function.
# P21, P22 = bbP2sph(i, j, s1pr, s2pr, q, qd, par)
function bbPhidist(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phi = (dot(d12, d12) - d^2) / 2
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phi = (dot(d12, d12) - d^2) / 2
    end
    return Phi
end

function bbPhiqdist(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * hcat(I3, BTran(p1, s1pr))
        Phiq2 = zeros(1, 7)
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * hcat(I3, BTran(p1, s1pr))
        Phiq2 = d12' * hcat(I3, BTran(p2, s2pr))
    end
    return Phiq1, Phiq2
end

function bbPhifxc(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phi = ((d12'*E11* d12) - d^2) / 2
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phi = ((d12'*E11* d12)  - d^2) / 2
    end
    Phi=(r1[1]*r1[1]-1) / 2
    return Phi
end

function bbPhiqfxc(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
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
    #println("bbPhiqfxc-r1",r1)
    Phiq1 = r1'*hcat(I3, zeros(3,4))
    Phiq1 = hcat(r1[1],zeros(1,2), zeros(1,4))
    #println("bbPhiqfxc",Phiq1)
    return Phiq1, Phiq2
end

function bbPhidot1(i, j, a1pr, a2pr,tn, q, par)
    # Unpack parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Get position and orientation of body i
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)

    Phi = 0.0
    if j == 0
        Phi = a1pr' * A1' * a2pr
    elseif j >= 1
        # Get position and orientation of body j
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        Phi = a1pr' * A1' * A2 * a2pr
    end

    return Phi
end

function bbPhidot2(i, j, a2pr, s1pr, s2pr,tn, q, par)
    # Unpack parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Get position and orientation of body i
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)

    Phi = 0.0
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phi = a2pr' * d12
    elseif j >= 1
        # Get position and orientation of body j
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phi = a2pr' * A2' * d12
    end

    return Phi
end
function bbPhiqdot1(i, j, a1pr, a2pr,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, a1pr)

    Phiq1 = zeros(1, 7)
    Phiq2 = zeros(1, 7)

    if j == 0
        Phiq1[4:7] = (a2pr' * BT1)'
    elseif j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        BT2 = BTran(p2, a2pr)
        Phiq1[4:7] = (a2pr' * A2' * BT1)'
        Phiq2[4:7] = (a1pr' * A1' * BT2)'
    end

    return Phiq1, Phiq2
end
function bbPhiqdot2(i, j, a2pr, s1pr, s2pr,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)

    Phiq1 = zeros(1, 7)
    Phiq2 = zeros(1, 7)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phiq1[1:7] = (-a2pr' * hcat(I3, BT1))'
    elseif j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        BT2s2 = BTran(p2, s2pr)
        BT2a2 = BTran(p2, a2pr)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phiq1[1:7] = (-a2pr' * A2' * hcat(I3, BT1))'
        Phiq2[1:7] = (a2pr' * A2' * hcat(I3, BT2s2) + d12' * hcat(zeros(3, 3), BT2a2))'
    end

    return Phiq1, Phiq2
end
function bbPhiRotDr(i, j, vx1pr, vy1pr, vx2pr, q, par)
    # Unpack parameters; parPart is assumed to be a function that extracts these values
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Compute s and c using the bbPhidot1 function which needs to be defined in Julia
    s = bbPhidot1(i, j, vy1pr, vx2pr,tn, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr,tn, q, par)

    # Use conditional statement to decide which value to return
    Phi = abs(c) >= abs(s) ? s : c

    return Phi
end

function bbPhiqRotDr(i, j, vx1pr, vy1pr, vx2pr,tn, q, par)
    # Unpack parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Compute scalar values of the constraint derivatives
    s = bbPhidot1(i, j, vy1pr, vx2pr,tn, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr,tn, q, par)

    # Decide which derivative to use based on the magnitude comparison
    if abs(c) >= abs(s)
        return bbPhiqdot1(i, j, vy1pr, vx2pr,tn, q, par)
    else
        return bbPhiqdot1(i, j, vx1pr, vx2pr,tn, q, par)
    end
end
function bbPhisph(i, j, s1pr, s2pr, tn, q, par)
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phi = d12
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phi = d12
    end

    return Phi
end


function bbPhiqsph(i, j, s1pr, s2pr,tn, q, par)
    # Unpack parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Identity matrix in Julia
    I3 = Matrix{Float64}(I, 3, 3)

    # Constraint Jacobians
    Phiq1 = zeros(3, 7)
    Phiq2 = zeros(3, 7)

    # Get the first part of the generalized coordinate for body i
    r1, p1 = qPart(q, i)

    if j == 0
        Phiq1[:, 1:7] = -hcat(I3, BTran(p1, s1pr))
        # Phiq2 remains zero as initialized
    elseif j >= 1
        r2, p2 = qPart(q, j)
        Phiq1[:, 1:7] = -hcat(I3, BTran(p1, s1pr))
        Phiq2[:, 1:7] = hcat(I3, BTran(p2, s2pr))
    end

    return Phiq1, Phiq2
end
