# mathfunction.jl
module mathfunction
using Test
using LinearAlgebra

export add_constraint!
export atilde
export ATran,BTran,qPart,bbP2dist,bbP2dot1,bbP2dot2,bbP2RotDr,bbP2sph
export bbPhidist,bbPhiqdist,bbPhidot1,bbPhidot2,bbPhiqdot1,bbPhiqdot2
export bbPhiqRotDr,bbPhiqsph,InitConfig,PhiqEval,PhiEval

# Define the function from the previous translation
function add_constraint!(A, B, m, n)
    r, s = size(B, 1), size(B, 2)
    for i in 1:r
        for j in 1:s
            A[m+i, n+j] += B[i, j]
        end
    end
    return A
end


function atilde(a)
    return [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
end
#= a = [1, 2, 3]
atilde_matrix = atilde(a)
println(atilde_matrix) =#


function ATran(p)
    e0 = p[1]
    e = p[2:4]
    I3 = I  # The identity matrix
    etil = atilde(e)
    AT = (e0^2 - dot(e, e)) * I3 + 2 * e * e' + 2 * e0 * etil
    return AT
end

#= p = [1, 2, 3,0]
at = ATran(p)
println(at) =#

function BTran(p, apr)
    e0 = p[1]
    e = p[2:4]
    I3 = I  # Identity matrix in Julia is represented by I
    etil = atilde(e)
    BT = 2 * hcat((e0 * I3 + etil) * apr, e * apr' - (e0 * I3 + etil) * atilde(apr))
    return BT
end

#= # Example usage:
p = [1, 2, 3, 4]  # Example quaternion components
apr = [5, 6, 7]   # Example vector
BT_matrix = BTran(p, apr)
println("The BTran matrix is:")
println(BT_matrix) =#

function qPart(q, i)
    r = q[7*(i-1)+1:7*(i-1)+3]
    p = q[7*(i-1)+4:7*(i-1)+7]
    return r, p
end

# Example usage:
q = rand(14) # Example vector of length 14, which can hold data for two 'i' indices
i = 2        # We want to extract the second set of (r, p)

#= r, p = qPart(q, i)
println("Vector r: ", r)
println("Quaternion p: ", p) =#

function parPart(par)
    nb = par[1]
    ngc = par[2]
    nh = par[3]
    nhc = par[4]
    nd = par[5]
    qtol = par[6]
    app = par[7]

    return nb, ngc, nh, nhc, nd, qtol, app
end

# Example usage:
par = [10, 20, 30, 40, 50, 0.001, 100]  # Example parameter array
nb, ngc, nh, nhc, nd, qtol, app = parPart(par)
println("nb: ", nb)
println("ngc: ", ngc)
println("nh: ", nh)
println("nhc: ", nhc)
println("nd: ", nd)
println("qtol: ", qtol)
println("app: ", app)

function bbP2dist(i, j, s1pr, s2pr, d, q, qd, par)
    # Assuming par is a structure or a dictionary with the following keys.
    nb, ngc, nh, nhc, nd, qtol, app = par[:nb], par[:ngc], par[:nh], par[:nhc], par[:nd], par[:qtol], par[:app]

    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1x = BTran(xp1, s1pr)
    a1 = (xr1' + xp1' * BT1)

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
        a2 = (xr2' + xp2' * BT2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P21 = [(a1 - a2), (a1 - a2) * BT1 - d12' * BT1x]
        P22 = [(a2 - a1), (a2 - a1) * BT2 + d12' * BT2x]
    end

    return P21, P22
end

function bbP2dot1(i, j, a1pr, a2pr, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)
    
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


function bbP2dot2(i, j, a2pr, s1pr, s2pr, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)
    
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


function bbP2RotDr(i, j, vx1pr, vy1pr, vx2pr, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

    s = bbPhidot1(i, j, vy1pr, vx2pr, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr, q, par)

    if abs(c) >= abs(s)
        P21, P22 = bbP2dot1(i, j, vy1pr, vx2pr, q, qd, par)
    else
        P21, P22 = bbP2dot1(i, j, vx1pr, vx2pr, q, qd, par)
    end

    return P21, P22
end

# Example usage:
# You would need to define the variables `i`, `j`, `vx1pr`, `vy1pr`, `vx2pr`, `q`, `qd`, and `par`
# appropriately before calling this function.
# P21, P22 = bbP2RotDr(i, j, vx1pr, vy1pr, vx2pr, q, qd, par)
function bbP2sph(i, j, s1pr, s2pr, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

    xr1, xp1 = qPart(qd, i)
    BT1x = BTran(xp1, s1pr)

    P21 = zeros(3, 7)
    P22 = zeros(3, 7)
    
    if j == 0
        P21[:, 4:6] = -BT1x
    elseif j >= 1
        xr2, xp2 = qPart(qd, j)
        BT2x = BTran(xp2, s2pr)
        P21[:, 4:6] = -BT1x
        P22[:, 4:6] = BT2x
    end

    return P21, P22
end

# Example usage:
# You would need to define the variables `i`, `j`, `s1pr`, `s2pr`, `q`, `qd`, and `par`
# appropriately before calling this function.
# P21, P22 = bbP2sph(i, j, s1pr, s2pr, q, qd, par)
function bbPhidist(i, j, s1pr, s2pr, d, q, par)
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)
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

function bbPhiqdist(i, j, s1pr, s2pr, d, q, par)
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * [I3, BTran(p1, s1pr)]
        Phiq2 = zeros(1, 7)
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * [I3, BTran(p1, s1pr)]
        Phiq2 = d12' * [I3, BTran(p2, s2pr)]
    end
    return Phiq1, Phiq2
end
function bbPhidot1(i, j, a1pr, a2pr, q, par)
    # Unpack parameters
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

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

function bbPhidot2(i, j, a2pr, s1pr, s2pr, q, par)
    # Unpack parameters
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

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
function bbPhiqdot1(i, j, a1pr, a2pr, q, par)
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

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
function bbPhiqdot2(i, j, a2pr, s1pr, s2pr, q, par)
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

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
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

    # Compute s and c using the bbPhidot1 function which needs to be defined in Julia
    s = bbPhidot1(i, j, vy1pr, vx2pr, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr, q, par)

    # Use conditional statement to decide which value to return
    Phi = abs(c) >= abs(s) ? s : c

    return Phi
end

function bbPhiqRotDr(i, j, vx1pr, vy1pr, vx2pr, q, par)
    # Unpack parameters
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

    # Compute scalar values of the constraint derivatives
    s = bbPhidot1(i, j, vy1pr, vx2pr, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr, q, par)

    # Decide which derivative to use based on the magnitude comparison
    if abs(c) >= abs(s)
        return bbPhiqdot1(i, j, vy1pr, vx2pr, q, par)
    else
        return bbPhiqdot1(i, j, vx1pr, vx2pr, q, par)
    end
end
function bbPhiqsph(i, j, s1pr, s2pr, q, par)
    # Unpack parameters
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

    # Identity matrix in Julia
    I3 = Matrix{Float64}(I, 3, 3)

    # Constraint Jacobians
    Phiq1 = zeros(3, 7)
    Phiq2 = zeros(3, 7)

    # Get the first part of the generalized coordinate for body i
    r1, p1 = qPart(q, i)

    if j == 0
        Phiq1[:, 1:7] = -[I3, BTran(p1, s1pr)]
        # Phiq2 remains zero as initialized
    elseif j >= 1
        r2, p2 = qPart(q, j)
        Phiq1[:, 1:7] = -[I3, BTran(p1, s1pr)]
        Phiq2[:, 1:7] = [I3, BTran(p2, s2pr)]
    end

    return Phiq1, Phiq2
end

function CylPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    ux1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    uz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    ux2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]
    uz2pr = [SJDT[20, k], SJDT[21, k], SJDT[22, k]]

    return i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr
end
function DistDrPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    return i, j, s1pr, s2pr
end
function DistPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    d = SJDT[10, k]
    return i, j, s1pr, s2pr, d
end
function EEval(p)
    e0 = p[1]
    e = p[2:4]
    Ebar = [-e atil(e) + e0 * I]  # I in Julia is the identity matrix
    return Ebar
end

# Helper function 'atil' that was used in EEval
function atil(a)
    return [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
end
function GamEval(tn, q, qd, SJDT, par)
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)
    P, Pst, Pstt = P5Eval(tn, q, SJDT, par)
    P2 = P2Eval(q, qd, SJDT, par)
    Gam = P2 * qd + Pstt
    return Gam
end
function GEval(x)
    e0 = x[1]
    e = x[2:4]
    etil = atil(e)
    G = [-e -etil + e0 * I]
    return G
end
function GamsqqdEval(tn, q, qd, SJDT, par)
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)
    P2 = P2Eval(tn, q, qd, SJDT, par)
    P3 = P3Eval(tn, q, qd, SJDT, par)
    Gamsq = P3
    Gamsqd = 2 * P2
    return Gamsq, Gamsqd
end

function InitConfig(rO, rP, rQ)
    # Define A using point definition of Sect 2.5.5.
    f = (1 / norm(rP - rO)) * (rP - rO)
    h = (1 / norm(atil(f) * (rQ - rO))) * atil(f) * (rQ - rO)
    g = -atil(f) * h
    A = [f g h]

    # Compute Euler parameter vector p in Section 2.5.4
    if norm(A' - A) > 1e-10
        trA = A[1, 1] + A[2, 2] + A[3, 3]
        e0 = 0.5 * sqrt(trA + 1)
        e = (1 / (2 * sqrt(trA + 1))) * [A[3, 2] - A[2, 3], A[1, 3] - A[3, 1], A[2, 1] - A[1, 2]]
    end
    if norm(A' - A) ≤ 1e-10
        e0 = 1
        e = zeros(3)
    end        
    p = [e0; e]

    # Evaluate q
    q = [rO; p]
    return q, p, A
end
function KEval(apr, b)
    K = 2 * [apr' * b apr' * atil(b);
              atil(apr) * b apr * b' + b * apr' - (apr' * b) * I]
    return K
end
function P2Eval(q, qd, SJDT, par)
    # Retrieve parameters from 'par'
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

    # Initialize P2 as a zero matrix of size ngc x ngc
    P2 = zeros(ngc, ngc)

    # Initialize joint number and constraint counter
    k = 1
    m = 0

    # Loop over each joint/constraint
    while k <= nh + nd
        constraintType = SJDT[1, k]
        # Process according to constraint type
        if constraintType == 1
            # Check if the constraint type is a Distance Constraint
            # Extract parameters for the Distance Constraint
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

            # Compute P21 and P22 for the Distance Constraint
            P21, P22 = bbP2dist(i, j, s1pr, s2pr, d, q, qd, par)

            # Add P21 to P2 at the appropriate location
            P2 = add_constraint!(P2, P21, m, 7 * (i - 1))

            # If j is not zero, add P22 as well
            if j >= 1
                P2 = add_constraint!(P2, P22, m, 7 * (j - 1))
            end

            # Increment the constraint counter
            m += 1
        # Check if the constraint type is a Spherical Constraint
        elseif constraintType == 2
            # Extract parameters for the Spherical Constraint
            i, j, s1pr, s2pr = SphPart(k, SJDT)

            # Compute P21 and P22 for the Spherical Constraint
            P21, P22 = bbP2sph(i, j, s1pr, s2pr, q, qd, par)

            # Add P21 to P2 at the appropriate location
            P2 = add_constraint!(P2, P21, m, 7 * (i - 1))

            # If j is not zero, add P22 as well
            if j >= 1
                P2 = add_constraint!(P2, P22, m, 7 * (j - 1))
            end

            # Increment the constraint counter by 3 for a Spherical Constraint
            m += 3
        # Cylindrical Constraint
        elseif constraintType == 3
            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = CylPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            P2cyl11, P2cyl12 = bbP2dot2(i, j, vx2pr, s1pr, s2pr, q, qd, par)
            P2cyl21, P2cyl22 = bbP2dot2(i, j, vy2pr, s1pr, s2pr, q, qd, par)
            P2cyl31, P2cyl32 = bbP2dot1(i, j, vz1pr, vx2pr, q, qd, par)
            P2cyl41, P2cyl42 = bbP2dot1(i, j, vz1pr, vy2pr, q, qd, par)

            P21k = vcat(P2cyl11, P2cyl21, P2cyl31, P2cyl41)
            P22k = vcat(P2cyl12, P2cyl22, P2cyl32, P2cyl42)

            P2 = add_constraint!(P2, P21k, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22k, m, 7 * (j - 1))
            end
            m += 4
        # Revolute Constraint
        elseif constraintType == 4
            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = RevPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            P2cyl11, P2cyl12 = bbP2dot2(i, j, vx2pr, s1pr, s2pr, q, qd, par)
            P2cyl21, P2cyl22 = bbP2dot2(i, j, vy2pr, s1pr, s2pr, q, qd, par)
            P2cyl31, P2cyl32 = bbP2dot1(i, j, vz1pr, vx2pr, q, qd, par)
            P2cyl41, P2cyl42 = bbP2dot1(i, j, vz1pr, vy2pr, q, qd, par)
            P2rev51, P2rev52 = bbP2dot2(i, j, vz2pr, s1pr, s2pr, q, qd, par)

            P21k = vcat(P2cyl11, P2cyl21, P2cyl31, P2cyl41, P2rev51)
            P22k = vcat(P2cyl12, P2cyl22, P2cyl32, P2cyl42, P2rev52)

            P2 = add_constraint!(P2, P21k, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22k, m, 7 * (j - 1))
            end
            m += 5

        # Translational Constraint
        elseif constraintType == 5
            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = TranPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            P2cyl11, P2cyl12 = bbP2dot2(i, j, vx2pr, s1pr, s2pr, q, qd, par)
            P2cyl21, P2cyl22 = bbP2dot2(i, j, vy2pr, s1pr, s2pr, q, qd, par)
            P2cyl31, P2cyl32 = bbP2dot1(i, j, vz1pr, vx2pr, q, qd, par)
            P2cyl41, P2cyl42 = bbP2dot1(i, j, vz1pr, vy2pr, q, qd, par)
            P2tran51, P2tran52 = bbP2dot1(i, j, vy1pr, vx2pr, q, qd, par)

            P21k = vcat(P2cyl11, P2cyl21, P2cyl31, P2cyl41, P2tran51)
            P22k = vcat(P2cyl12, P2cyl22, P2cyl32, P2cyl42, P2tran52)

            P2 = add_constraint!(P2, P21k, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22k, m, 7 * (j - 1))
            end
            m += 5

        elseif constraintType == 10
            # Rotation Driver
            i, j, vx1pr, vy1pr, vx2pr = RotDrPart(k, SJDT)
            P21, P22 = bbP2RotDr(i, j, vx1pr, vy1pr, vx2pr, q, qd, par)
            P2 = add_constraint!(P2, P21, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22, m, 7 * (j - 1))
            end
            m += 1
        end



        k += 1
    end

    # Euler Parameter Normalization Constraints
    i = 1
    while i <= nb
        xr1, xp1 = qPart(qd, i)
        P21 = hcat(zeros(1, 3),xp1')
        P2 = add_constraint!(P2, P21, m, 7 * (i - 1))
        i += 1
        m += 1
    end

    return P2
end

function P5Eval(tn, q, SJDT, par)
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)  # Unpack parameters

    P = zeros(nhc + 2 + nb, 1)
    Pst = zeros(nhc + 2 + nb, 1)
    Pstt = zeros(nhc + 2 + nb, 1)

    if app == 1  # 2 Bar with 2 rotational drivers
        # Rot driver body1 to ground
        k = 3
        i, j, u1pr, v1pr, u2pr = RotDrPart(k, SJDT)
        omega = 10
        theta10 = -omega * tn  # Angle from bod 1, to ground
        theta10d = -omega
        theta10dd = 0
        r1, p1 = qPart(q, i)
        A1 = ATran(p1)
        PD, PDd, PDdd = computeDriverInputs(u1pr, v1pr, u2pr, A1, theta10, theta10d, theta10dd, j, q)

        PD1 = PD
        PDd1 = PDd
        PDdd1 = PDdd

        # Rot driver body1 to body 2
        k = 4
        i, j, u1pr, v1pr, u2pr = RotDrPart(k, SJDT)
        theta12 = omega * tn  # Angle from bod 1, to bod 2
        theta12d = omega
        theta12dd = 0
        r1, p1 = qPart(q, i)
        A1 = ATran(p1)
        PD, PDd, PDdd = computeDriverInputs(u1pr, v1pr, u2pr, A1, theta12, theta12d, theta12dd, j, q)

        PD2 = PD
        PDd2 = PDd
        PDdd2 = PDdd

        P = vcat(zeros(nhc, 1), PD1, PD2, zeros(nb, 1))
        Pst = vcat(zeros(nhc, 1), PDd1, PDd2, zeros(nb, 1))
        Pstt = vcat(zeros(nhc, 1), PDdd1, PDdd2, zeros(nb, 1))
    end

    return P, Pst, Pstt
end

function computeDriverInputs(u1pr, v1pr, u2pr, A1, theta, thetad, thetadd, j, q)
    c = s = PD = PDd = PDdd = 0
    if j == 0
        c = dot(u1pr, A1' * u2pr)
        s = dot(v1pr, A1' * u2pr)
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        c = dot(u1pr, A1' * A2 * u2pr)
        s = dot(v1pr, A1' * A2 * u2pr)
    end

    if abs(c) >= abs(s)
        PD = -sin(theta)
        PDd = -thetad * cos(theta)
        PDdd = -thetadd * cos(theta) + thetad^2 * sin(theta)
    else
        PD = -cos(theta)
        PDd = thetad * sin(theta)
        PDdd = thetadd * sin(theta) + thetad^2 * cos(theta)
    end

    return PD, PDd, PDdd
end




function PhiEval(tn, q, SJDT, par)
    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)
    Phi = zeros(ngc)
    I3 = Matrix{Float64}(I, 3, 3)
    k = 1        # Joint No.
    m = 0        # Constraint Counter - 1
    while k <= nh + nd
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            Phik = bbPhidist(i, j, s1pr, s2pr, d, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end

        # Spherical Constraint
        if SJDT[1, k] == 2
            i, j, s1pr, s2pr = SphPart(k, SJDT)
            Phik = bbPhisph(i, j, s1pr, s2pr, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 3
        end

        # Cylindrical Constraint
        if SJDT[1, k] == 3
            i, j, s1pr, s2pr, u1pr, vz1pr, u2pr, vz2pr = CylPart(k, SJDT)
            v1pr = atil(vz1pr) * u1pr
            vy2pr = atil(vz2pr) * u2pr
            Phicyl1 = bbPhidot2(i, j, u2pr, s1pr, s2pr, q, par)
            Phicyl2 = bbPhidot2(i, j, vy2pr, s1pr, s2pr, q, par)
            Phicyl3 = bbPhidot1(i, j, vz1pr, u2pr, q, par)
            Phicyl4 = bbPhidot1(i, j, vz1pr, vy2pr, q, par)
            Phik = vcat(Phicyl1, Phicyl2, Phicyl3, Phicyl4)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 4
        end
        # Revolute Constraint
        if SJDT[1, k] == 4
            i, j, s1pr, s2pr, u1pr, vz1pr, u2pr, vz2pr = RevPart(k, SJDT)
            v1pr = atil(vz1pr) * u1pr
            vy2pr = atil(vz2pr) * u2pr

            Phicyl1 = bbPhidot2(i, j, u2pr, s1pr, s2pr, q, par)
            Phicyl2 = bbPhidot2(i, j, vy2pr, s1pr, s2pr, q, par)
            Phicyl3 = bbPhidot1(i, j, vz1pr, u2pr, q, par)
            Phicyl4 = bbPhidot1(i, j, vz1pr, vy2pr, q, par)
            Phirev5 = bbPhidot2(i, j, vz2pr, s1pr, s2pr, q, par)

            Phik = vcat(Phicyl1, Phicyl2, Phicyl3, Phicyl4, Phirev5)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 5
        end
        # Rotation Driver
        if SJDT[1, k] == 10
            i, j, u1pr, v1pr, u2pr = RotDrPart(k, SJDT)
            Phik = bbPhiRotDr(i, j, u1pr, v1pr, u2pr, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end

        # ... (Continue with the rest of the constraints as per the given code)

        k += 1
    end

    # Euler Parameter Normalization Constraints
    j = 1
    while j <= nb
        r, p = qPart(q, j)
        Phik = (dot(p, p) - 1) / 2
        Phi = add_constraint!(Phi, Phik, m, 0)
        j += 1
        m += 1
    end

    P, Pst, Pstt = P5Eval(tn, q, SJDT, par)
    Phi += P

    return Phi
end




function PhiqEval(q, SJDT, par)

    nb, ngc, nh, nhc, nd, qtol, app = parPart(par)

    Phiq = zeros(Float64,ngc, ngc)

    I3 = Matrix{Float64}(I, 3, 3)
    k = 1  # Joint number
    m = 0  # Constraint equation counter - 1

    while k <= nh + nd
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqdist(i, j, s1pr, s2pr, d, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        # Revolute Constraint
        if SJDT[1, k] == 4
            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = RevPart(k, SJDT)
            vy1pr = cross(vz1pr, vx1pr)
            vy2pr = cross(vz2pr, vx2pr)
        
            Phiqcyl11, Phiqcyl12 = bbPhiqdot2(i, j, vx2pr, s1pr, s2pr, q, par)
            Phiqcyl21, Phiqcyl22 = bbPhiqdot2(i, j, vy2pr, s1pr, s2pr, q, par)
            Phiqcyl31, Phiqcyl32 = bbPhiqdot1(i, j, vz1pr, vx2pr, q, par)
            Phiqcyl41, Phiqcyl42 = bbPhiqdot1(i, j, vz1pr, vy2pr, q, par)
            Phiqrev51, Phiqrev52 = bbPhiqdot2(i, j, vz2pr, s1pr, s2pr, q, par)
        
            Phiq1k = vcat(Phiqcyl11, Phiqcyl21, Phiqcyl31, Phiqcyl41, Phiqrev51)
            Phiq2k = vcat(Phiqcyl12, Phiqcyl22, Phiqcyl32, Phiqcyl42, Phiqrev52)
        
            Phiq = add_constraint!(Phiq, Phiq1k, m, 7 * (i - 1))
        
            if j ≥ 1
                Phiq = add_constraint!(Phiq, Phiq2k, m, 7 * (j - 1))
            end
        
            m += 5
        end
        # Rotation Driver
        if SJDT[1, k] == 10
            i, j, vx1pr, vy1pr, vx2pr = RotDrPart(k, SJDT)
            Phiq1, Phiq2 = bbPhiqRotDr(i, j, vx1pr, vy1pr, vx2pr, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7 * (i - 1))
        
            if j ≥ 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7 * (j - 1))
            end
            m += 1
        end

        # ... (Continue with the rest of the constraints as per the given code)

        k += 1
    end

    # Euler Parameter Normalization Constraints
    i = 1
    while i <= nb
        r1, p1 = qPart(q, i)
        Phiq1 = hcat(zeros(1, 3),p1')
        Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))
        m += 1
        i += 1
    end

    return Phiq
end

function PhiqEval(tn, q, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    Phiq = zeros(nc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)
    k = 1  # Joint No.
    m = 0  # Constraint equation Counter - 1

    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            Phiq1, Phiq2 = bbPhiqdist(i, j, s1pr, s2pr, d, tn, q, par)
            Phiq = Add(Phiq, Phiq1, m, 7 * (i - 1))

            if j >= 1
                Phiq = Add(Phiq, Phiq2, m, 7 * (j - 1))
            end

            m += 1
        end

        # [Continue with other constraints, similar to MATLAB but adapted to Julia syntax]

        k += 1
    end

    # Euler Parameter Normalization Constraints
    i = 1
    while i <= nb
        r1, p1 = qPart(q, i)
        Phiq1 = [zeros(1, 3) p1']
        Phiq = Add(Phiq, Phiq1, m, 7 * (i - 1))
        m += 1
        i += 1
    end

    return Phiq
end


function RevPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    vx1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    vz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    vx2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]
    vz2pr = [SJDT[20, k], SJDT[21, k], SJDT[22, k]]

    return i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr
end

function RotDrPart(k, SJDT)
    # Data derived from the host revolute or cylindrical joint data table
    i = SJDT[2, k]
    j = SJDT[3, k]
    vx1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    vz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    vy1pr = cross(vz1pr, vx1pr)  # using cross product to get the perpendicular vector
    vx2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]

    return i, j, vx1pr, vy1pr, vx2pr
end

# You would need to implement or provide the equivalent methods such as parPart, DistPart, Add, bbPhiqdist, qPart, etc.


end # module mathfunction