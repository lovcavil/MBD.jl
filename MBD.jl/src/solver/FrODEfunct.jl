function FrODEfunct(n, tn, v, vd, vdd, Lam, u, q0, V, U, B, SJDT, SMDT, STSDAT, par)
    nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)

    # Evaluate terms that do not depend on vdd and Lam
    u, Iteru = usolv(tn, u, v, q0, SJDT, V, U, B, par)
    q = q0 .+ V * v .- U * u
    Phiq = PhiqEval(tn, q, SJDT, par)
    D = (Matrix(I, ngc, ngc) - U * B * Phiq) * V
    qd = D * vd
    M = MEval(q, SMDT, par)
    Gam = GamEval(tn, q, qd, SJDT, par)
    B, Biter = CorrectB(tn, q, B, U, SJDT, par)
    Pvdd = [Matrix(I, nv, nv) zeros(nv, nc)]
    PLam = [zeros(nc, nv) Matrix(I, nc, nc)]

    # Solve for vdd and Lam
    i = 1
    err = intol + 1
    E=[]
    while err > intol
        # Jacobian Evaluation
        QAsLam = QAsLamEval(tn, q, qd, Lam, SJDT, STSDAT, par)
        E = [M * D Phiq' - QAsLam]

        # Residual Calculation
        QA = QAEval(tn,q,qd,Lam,SMDT,SJDT,STSDAT,par)
        S = SEval(q, qd, SMDT, par)
        R = M * D * vdd + Phiq' * Lam - M * U * B * Gam - QA - S

        # Newton Correction
        x = -E \ R
        vdd += Pvdd * x
        Lam += PLam * x

        err = norm(R) + norm(x)
        i += 1
    end

    jodeiter = i - 1
    ECond = cond(E)

    return vdd, Lam, jodeiter, ECond
end
