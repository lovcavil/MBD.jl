
function GamEval(tn, q, qd, SJDT, par)
    Pst, Pstt, Pstq, Psttq = P5Eval(tn, q, par)
    P2 = P2Eval(tn, q, qd, SJDT, par)
    Gam = P2 * qd + Pstt

    return Gam
end

function GamsqqdEval(tn, q, qd, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    P2 = P2Eval(tn, q, qd, SJDT, par)
    P3 = P3Eval(tn, q, qd, SJDT, par)
    Gamsq = P3
    Gamsqd = 2 * P2
    return Gamsq, Gamsqd
end

function KEval(apr, b)
    # Assuming `atil` function exists or is replaced by an equivalent Julia function
    K = 2 * vcat(hcat(dot(apr, b),     apr' * atil(b)),
             hcat(atil(apr) * b,   apr * b' + b * apr' - dot(apr, b) * I))
    return K
end


function MEval(q, SMDT, par)
    nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt=parPart(par)
    M = zeros(ngc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix

    for i in 1:nb
        m = SMDT[1, i]
        J = Diagonal([SMDT[2, i], SMDT[3, i], SMDT[4, i]])
        r, p = qPart(q, i)
        G = GEval(p)
        Mi = BlockDiagonal([m * I3, 4 * G' * J * G])  # blockdiag function to create block diagonal matrix
        M = add_constraint!(M, Mi, 7 * (i - 1), 7 * (i - 1))
    end

    return M
end

function M2Eval(q, mu, SMDT, par)
    nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt=parPart(par)

    M2 = zeros(ngc, ngc)

    for i in 1:nb
        r, p = qPart(q, i)
        mur, mup = qPart(mu, i)
        m = SMDT[1, i]
        J = Diagonal([SMDT[2, i], SMDT[3, i], SMDT[4, i]])
        G = GEval(p)
        Gmu = GEval(mup)
        M2i = vcat(zeros(3, 7), hcat(zeros(4, 3), TEval(4 * J * G * mup) - 4 * G' * J * Gmu))
        M2 = add_constraint!(M2, M2i, 7 * (i - 1), 7 * (i - 1))
    end

    return M2
end


function SEval(q, qd, SMDT, par)
    nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt = parPart(par)
    S = zeros(ngc)
    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix in Julia
    for i in 1:nb
        m = SMDT[1, i]
        J = diagm([SMDT[2, i], SMDT[3, i], SMDT[4, i]])
        r, p = qPart(q, i)
        rd, pd = qPart(qd, i)
        G = GEval(p)
        Gd = GEval(pd)
        Si = vcat([0; 0; 0], 8 * Gd' * J * Gd * p)  # Adjust if necessary for correct matrix/vector sizes
        # Add Si to S (Assuming 'Add' is defined elsewhere)
        S = add_constraint!(S, Si, 7 * (i - 1), 0)  # Adjust this line if the 'Add' function works differently in Julia
    end
    return S
end


function SqqdEval(q, qd, SMDT, par)
    nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt = parPart(par)

    Ssq = zeros(ngc, ngc)
    Ssqd = zeros(ngc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix of size 3x3

    for i = 1:nb
        m = SMDT[1, i]
        J = Diagonal([SMDT[2, i], SMDT[3, i], SMDT[4, i]])
        r, p = qPart(q, i)
        rd, pd = qPart(qd, i)
        G = GEval(p)
        Gd = GEval(pd)
        Ssqi = vcat(zeros(3, 7), hcat(zeros(4, 3),8 * Gd' * J * Gd))
        Ssq = add_constraint!(Ssq, Ssqi, 7 * (i - 1), 7 * (i - 1))
        Ssqdi = vcat(zeros(3, 7), hcat(zeros(4, 3),8 * TEval(J * Gd * p) - 8 * Gd' * J * G))
        Ssqd = add_constraint!(Ssqd, Ssqdi, 7 * (i - 1), 7 * (i - 1))
    end

    return Ssq, Ssqd
end
