# mathfunction.jl

using Test
using LinearAlgebra
using BlockDiagonals
using CSV, DataFrames
include("./partition/partition_II8.jl")
include("./eval/eval_II8.jl")
include("./eval/P2.jl")
include("./eval/PhiEval.jl")
include("./eval/PhiqEval.jl")
include("./eval/QA.jl")
include("./eval/QAsLamEval.jl")
include("./eval/QAsqqd.jl")
include("./eval/QAwEval.jl")
include("./eval/QAwsLam.jl")
include("./eval/SfrSfrpr.jl")
include("./constraint/constraint_II8.jl")
include("./constraint/constraint_fPy.jl")
include("./constraint/constraint_fPz.jl")
include("./constraint/constraint_fP.jl")
include("./constraint/constraint_Poly.jl")
include("./constraint/constraint_Spline.jl")
export add_constraint!
export atilde
export ATran,BTran,qPart,bbP2dist,bbP2dot1,bbP2dot2,bbP2RotDr,bbP2sph
export bbPhidist,bbPhiqdist,bbPhidot1,bbPhidot2,bbPhiqdot1,bbPhiqdot2
export bbPhiqRotDr,bbPhiqsph,InitConfig

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

function CTran(p, apr)
    e0 = p[1]
    e = p[2:4]
    etil = atil(e)  # Assuming 'atil' is previously defined
    CT = 2 * hcat((e0 * I - etil) * apr, e * apr' + (e0 * I - etil) * atil(apr))
    return CT
end

#= # Example usage:
p = [1, 2, 3, 4]  # Example quaternion components
apr = [5, 6, 7]   # Example vector
BT_matrix = BTran(p, apr)
println("The BTran matrix is:")
println(BT_matrix) =#



# Example usage:
#q = rand(14) # Example vector of length 14, which can hold data for two 'i' indices
#i = 2        # We want to extract the second set of (r, p)

#= r, p = qPart(q, i)
println("Vector r: ", r)
println("Quaternion p: ", p) =#





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


function GEval(x)
    e0 = x[1]
    e = x[2:4]
    etil = atil(e)
    G = hcat(-e,-etil + e0 * I)
    return G
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


function TEval(a)
    T =vcat( hcat(0,-a'),hcat( a,-atil(a)))
    return T
end
function ODEfunct(tn, q, qd, SMDT, STSDAT, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    Gam = GamEval(tn, q, qd, SJDT, par)
    QA = QAEval(tn, q, qd, SMDT, STSDAT, par)
    S = SEval(q, qd, SMDT, par)
    RHS = vcat(QA + S, -Gam)
    M = MEval(q, SMDT, par)
    Phiq = PhiqEval(tn, q, SJDT, par)
    E = vcat(hcat(M,Phiq'), hcat(Phiq,zeros(nc, nc)))
    ECond = cond(E)
    println("cond",ECond)
    println("r",rank(E))
    # x = E \ RHS
    x = pinv(E) * RHS

    qdd = zeros(ngc)
    for i in 1:ngc
        qdd[i] = x[i]
    end
    Lam = zeros(nc)
    for i in 1:nc
        Lam[i] = x[ngc + i]
    end

    return qdd, Lam, ECond
end


function P3Eval(tn, q, qd, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    P3 = zeros(nc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix of size 3x3
    k = 1   # Joint number
    m = 0   # Constraint counter - 1

    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            P31, P32 = bbP3dist(i, j, s1pr, s2pr, d, tn, q, qd, par)
            P3 = add_constraint!(P3, P31, m, 7 * (i - 1))

            if j >= 1
                P3 = add_constraint!(P3, P32, m, 7 * (j - 1))
            end

            m = m + 1
        end

        # Spherical Constraint
        if SJDT[1, k] == 2
            i, j, s1pr, s2pr = SphPart(k, SJDT)
            # P3 spherical contribution is zero, no nonzero terms to add

            m = m + 3
        end

        # Rev Constraint
        if SJDT[1, k] == 4
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = RevPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            P3cyl11, P3cyl12 = bbP3dot2(i, j, ux2pr, s1pr, s2pr, tn, q, qd, par)
            P3cyl21, P3cyl22 = bbP3dot2(i, j, uy2pr, s1pr, s2pr, tn, q, qd, par)
            P3cyl31, P3cyl32 = bbP3dot1(i, j, uz1pr, ux2pr, tn, q, qd, par)
            P3cyl41, P3cyl42 = bbP3dot1(i, j, uz1pr, uy2pr, tn, q, qd, par)
            P3rev51, P3rev52 = bbP3dot2(i, j, uz2pr, s1pr, s2pr, tn, q, qd, par)
        
            P31k = [P3cyl11; P3cyl21; P3cyl31; P3cyl41; P3rev51]
            P32k = [P3cyl12; P3cyl22; P3cyl32; P3cyl42; P3rev52]
        
            P3 = add_constraint!(P3, P31k, m, 7 * (i - 1))
        
            if j >= 1
                P3 = add_constraint!(P3, P32k, m, 7 * (j - 1))
            end
        
            m += 5
        end
        # Translational Constraint
        if SJDT[1, k] == 5
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = TranPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            P3cyl11, P3cyl12 = bbP3dot2(i, j, ux2pr, s1pr, s2pr, tn, q, qd, par)
            P3cyl21, P3cyl22 = bbP3dot2(i, j, uy2pr, s1pr, s2pr, tn, q, qd, par)
            P3cyl31, P3cyl32 = bbP3dot1(i, j, uz1pr, ux2pr, tn, q, qd, par)
            P3cyl41, P3cyl42 = bbP3dot1(i, j, uz1pr, uy2pr, tn, q, qd, par)
            P3tran51, P3tran52 = bbP3dot1(i, j, uy1pr, ux2pr, tn, q, qd, par)
        
            P31k = [P3cyl11; P3cyl21; P3cyl31; P3cyl41; P3tran51]
            P32k = [P3cyl12; P3cyl22; P3cyl32; P3cyl42; P3tran52]
        
            P3 = add_constraint!(P3, P31k, m, 7 * (i - 1))
        
            if j >= 1
                P3 = add_constraint!(P3, P32k, m, 7 * (j - 1))
            end
        
            m += 5
        end
        # fxc Constraint
        if SJDT[1, k] == 1010
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            P31, P32 = bbP3fxc(i, j, s1pr, s2pr, d, tn, q, qd, par)
            P3 = add_constraint!(P3, P31, m, 7 * (i - 1))

            if j >= 1
                P3 = add_constraint!(P3, P32, m, 7 * (j - 1))
            end

            m = m + 1
        end                
        if SJDT[1, k] == 1020
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            P31, P32 = bbP3fy(i, j, s1pr, s2pr, d, tn, q, qd, par)
            P3 = add_constraint!(P3, P31, m, 7 * (i - 1))

            if j >= 1
                P3 = add_constraint!(P3, P32, m, 7 * (j - 1))
            end

            m = m + 1
        end           
        k = k + 1
    end

    return P3
end



function P4Eval(tn, q, eta, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    P4 = zeros(ngc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)
    Z1 = zeros(3, 1)
    Z3 = zeros(3, 3)
    k = 1    # Joint number
    m = 0    # Address in vector eta

    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            etak = eta[m + 1]

            P411, P412, P422 = bbP4dist(i, j, s1pr, s2pr, d, tn, q, etak, par)
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))

            if j >= 1
                P421 = transpose(P412)
                P4 = add_constraint!(P4, P412, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 1
        end

        # Spherical Constraint
        if SJDT[1, k] == 2
            i, j, s1pr, s2pr = SphPart(k, SJDT)
            etak = eta[m + 1 : m + 3]

            P411, P412, P422 = bbP4sph(i, j, s1pr, s2pr, tn, q, etak, par)    
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))    

            if j >= 1
                # P412 and P421 are zero
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 3
        end
        # Rev
        if SJDT[1, k] == 4
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = RevPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            P4cyl111, P4cyl112, P4cyl122 = bbP4dot2(i, j, ux2pr, s1pr, s2pr, tn, q, eta[m + 1], par)
            P4cyl211, P4cyl212, P4cyl222 = bbP4dot2(i, j, uy2pr, s1pr, s2pr, tn, q, eta[m + 2], par)
            P4cyl311, P4cyl312, P4cyl322 = bbP4dot1(i, j, uz1pr, ux2pr, tn, q, eta[m + 3], par)
            P4cyl411, P4cyl412, P4cyl422 = bbP4dot1(i, j, uz1pr, uy2pr, tn, q, eta[m + 4], par)
            P4rev511, P4rev512, P4rev522 = bbP4dot2(i, j, uz2pr, s1pr, s2pr, tn, q, eta[m + 5], par)
        
            P411k = P4cyl111 + P4cyl211 + P4cyl311 + P4cyl411 + P4rev511
            P412k = P4cyl112 + P4cyl212 + P4cyl312 + P4cyl412 + P4rev512
            P421k = transpose(P412k)
            P422k = P4cyl122 + P4cyl222 + P4cyl322 + P4cyl422 + P4rev522
        
            P4 = add_constraint!(P4, P411k, 7 * (i - 1), 7 * (i - 1))
        
            if j >= 1
                P4 = add_constraint!(P4, P412k, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421k, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422k, 7 * (j - 1), 7 * (j - 1))
            end
        
            m += 5
        end
        # Tra
        if SJDT[1, k] == 5
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = TranPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            P4cyl111, P4cyl112, P4cyl122 = bbP4dot2(i, j, ux2pr, s1pr, s2pr, tn, q, eta[m + 1], par)
            P4cyl211, P4cyl212, P4cyl222 = bbP4dot2(i, j, uy2pr, s1pr, s2pr, tn, q, eta[m + 2], par)
            P4cyl311, P4cyl312, P4cyl322 = bbP4dot1(i, j, uz1pr, ux2pr, tn, q, eta[m + 3], par)
            P4cyl411, P4cyl412, P4cyl422 = bbP4dot1(i, j, uz1pr, uy2pr, tn, q, eta[m + 4], par)
            P4tran511, P4tran512, P4tran522 = bbP4dot1(i, j, uy1pr, ux2pr, tn, q, eta[m + 3], par)
        
            P411k = P4cyl111 + P4cyl211 + P4cyl311 + P4cyl411 + P4tran511
            P412k = P4cyl112 + P4cyl212 + P4cyl312 + P4cyl412 + P4tran512
            P421k = transpose(P412k)
            P422k = P4cyl122 + P4cyl222 + P4cyl322 + P4cyl422 + P4tran522
        
            P4 = add_constraint!(P4, P411k, 7 * (i - 1), 7 * (i - 1))
        
            if j >= 1
                P4 = add_constraint!(P4, P412k, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421k, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422k, 7 * (j - 1), 7 * (j - 1))
            end
        
            m += 5
        end
        # fxc Constraint
        if SJDT[1, k] == 1010
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            etak = eta[m + 1]

            P411, P412, P422 = bbP4fxc(i, j, s1pr, s2pr, d, tn, q, etak, par)
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))

            if j >= 1
                P421 = transpose(P412)
                P4 = add_constraint!(P4, P412, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 1
        end        
        if SJDT[1, k] == 1020
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            etak = eta[m + 1]

            P411, P412, P422 = bbP4fy(i, j, s1pr, s2pr, d, tn, q, etak, par)
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))

            if j >= 1
                P421 = transpose(P412)
                P4 = add_constraint!(P4, P412, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 1
        end 
        if SJDT[1, k] == 1030
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            etak = eta[m + 1]

            P411, P412, P422 = bbP4fz(i, j, s1pr, s2pr, d, tn, q, etak, par)
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))

            if j >= 1
                P421 = transpose(P412)
                P4 = add_constraint!(P4, P412, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 1
        end       
        k += 1
    end

    # Euler Parameter Normalization Constraint
    i = 1
    while i <= nb
        etak = eta[m + 1]
        P411 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3),I))
        P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))
        i += 1
        m += 1
    end

    return P4
end


function P5Eval(tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Enter Constraint t derivatives of P(t, q, qd, par); Default is Zeros
    Pst = zeros(nc)
    Pstt = zeros(nc)
    Pstq = zeros(nc, ngc)
    Psttq = zeros(nc, ngc)

    return Pst, Pstt, Pstq, Psttq
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

function csign(a, par)
    nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)
    csa=0
    dcsa=0
    # Continuous sign function
    if abs(a) > vt
        csa = sign(a)
    elseif abs(a) <= vt
        csa = sin(pi * a / (2 * vt))
    end

    # Derivative of continuous sign function
    if abs(a) > vt
        dcsa = 0
    elseif abs(a) <= vt
        dcsa = (pi / (2 * vt)) * cos(pi * a / (2 * vt))
    end

    return csa, dcsa
end

function Param(n, tnm, Q, Qd, Qdd, SJDT, par, jRepar)
    nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)

    jRepar += 1
    q0 = Q[:, n - 1]
    Phiq = PhiqEval(tnm, q0, SJDT, par)
    U = Phiq'
    B = inv(U' * U)
    V = nullspace(Phiq)
    v = zeros(nv)
    qd = Qd[:, n - 1]
    vd = V' * qd
    qdd = Qdd[:, n - 1]
    vdd = V' * qdd

    return v, vd, vdd, q0, U, V, B, jRepar
end


function CorrectB(tn, q, B, U, SJDT, par)
    nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)

    i = 1
    Cnorm = Btol + 1
    #I = Matrix(I, nu, nu)  # Identity matrix
    Phiq = PhiqEval(tn, q, SJDT, par)
    
    while Cnorm > Btol
        delB = -B * Phiq * U * B + B
        B += delB
        C = B * Phiq * U - I
        Cnorm = norm(C)
        i += 1
    end

    Biter = i - 1

    return B, Biter
end
