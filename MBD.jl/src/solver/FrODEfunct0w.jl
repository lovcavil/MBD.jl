function FrODEfunct0w(q0, qd0, vdd, Lam, V, U, B, SJDT, SMDT, STSDAT, par, w, N)

    nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)

    # Evaluate terms that do not depend on vdd and Lam
    Phiq = PhiqEval(0, q0, SJDT, par)
    D = (I - U * B * Phiq) * V
    M = MEval(q0, SMDT, par)
    Gam = GamEval(0, q0, qd0, SJDT, par)
    Pvdd = hcat(I(nv), zeros(nv, nc))
    PLam = hcat(zeros(nc, nv), I(nc))
    S = SEval(q0, qd0, SMDT, par)

    # Solve for vdd and Lam
    i = 1  # Set solution iteration counter
    err = intol + 1
    E=[]
    while err > intol
        # Jacobian Evaluation
        QAwsLam = QAwsLamEval(0, q0, qd0, Lam, SJDT, STSDAT, par, w, N)
        E = hcat(M * D, Phiq' - QAwsLam)

        # Residual Calculation
        QAw = QAwEval(0, q0, qd0, Lam, SMDT, SJDT, STSDAT, par, w, N)
        R = M * D * vdd + Phiq' * Lam - M * U * B * Gam - QAw - S

        # Newton Correction
        x = -E \ R
        vdd += Pvdd * x
        Lam += PLam * x

        err = norm(R)
        i += 1
    end

    jodeiter = i - 1
    ECond = cond(E)

    return vdd, Lam, jodeiter, ECond
end
