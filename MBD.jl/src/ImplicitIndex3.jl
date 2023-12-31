include("mathfunction_II7.jl")
using CSV, DataFrames
#using mathfunction: einv
function ImplicitIndex3(n, tn, Q, Qd, Qdd, LLam, h, hmax,
    SMDT, STSDAT, SJDT, par, alpha, nch)

    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = mathfunction.parPart(par)

    beta = ((1 - alpha)^2) / 4
    delta = (1 - 2 * alpha) / 2

    tnm = tn - h
    qnm = Q[:, n-1]
    qdnm = Qd[:, n-1]
    qddnm = Qdd[:, n-1]
    Lamnm = LLam[:, n-1]
    q = Q[:, n-1]
    qd = Qd[:, n-1]
    qdd = Qdd[:, n-1]
    Lam = LLam[:, n-1]
    R1n = 0
    JCond = 0
    Jiter = 0
    hopt=0
    Err = 2   # Criteria for accepting time step
    while Err > 1

        Phiq = mathfunction.PhiqEval(tnm, qnm, SJDT, par)
        Gamsq, Gamsqd = mathfunction.GamsqqdEval(tnm, qnm, qdnm, SJDT, par)
        P2 = mathfunction.P2Eval(tnm, qnm, qddnm, SJDT, par)
        M2 = mathfunction.M2Eval(qnm, qddnm, SMDT, par)
        QAsq, QAsqd = mathfunction.QAsqqd(tnm, qnm, qdnm, SMDT, STSDAT, par)
        Ssq, Ssqd = mathfunction.SqqdEval(qnm, qdnm, SMDT, par)
        P4 = mathfunction.P4Eval(tnm, qnm, Lamnm, SJDT, par)
        QA = mathfunction.QAEval(tnm, qnm, qdnm, SMDT, STSDAT, par)
        S = mathfunction.SEval(qnm, qdnm, SMDT, par)
        M = mathfunction.MEval(qnm, SMDT, par)

        # Term in HHT at tnm
        R1pnm = -(alpha / (1 + alpha)) * (Phiq' * Lamnm - QA - S)

        # Third Derivative calculation
        E = vcat(hcat(M, Phiq'), hcat(Phiq, zeros(nc, nc)))
        Rhs = vcat((-M2 - P4 + QAsq + Ssq) * qdnm + (QAsqd + Ssqd) * qddnm,
            (-P2 - Gamsq) * qdnm - Gamsqd * qddnm)
        
        df = DataFrame(E, :auto)
        # Write to CSV
        CSV.write("ii3_E.csv", df)
        x = E \ Rhs
        Pqdd = hcat(I, zeros(ngc, nc))
        PLam = hcat(zeros(nc, ngc), I)
        qdddnm = Pqdd * x
        Lamdnm = PLam * x

        # Estimate for qdd and Lam
        qdd = qddnm + h * qdddnm
        Lam = Lamnm + h * Lamdnm

        qd = qdnm + h * ((1 - delta) * qddnm + delta * qdd)
        q = qnm + h * qdnm + ((h^2) / 2) * ((1 - 2 * beta) * qddnm + 2 * beta * qdd)

        Phiq = mathfunction.PhiqEval(tn, q, SJDT, par)
        Gam = mathfunction.GamEval(tn, q, qd, SJDT, par)
        Gamsq, Gamsqd = mathfunction.GamsqqdEval(tn, q, qd, SJDT, par)
        P2 = mathfunction.P2Eval(tn, q, qdd, SJDT, par)
        M = mathfunction.MEval(q, SMDT, par)
        M2 = mathfunction.M2Eval(q, qdd, SMDT, par)
        QA = mathfunction.QAEval(tn, q, qd, SMDT, STSDAT, par)
        QAsq, QAsqd = mathfunction.QAsqqd(tn, q, qd, SMDT, STSDAT, par)
        S = mathfunction.SEval(q, qd, SMDT, par)
        Ssq, Ssqd = mathfunction.SqqdEval(q, qd, SMDT, par)
        P4 = mathfunction.P4Eval(tn, q, Lam, SJDT, par)

        # Jacobian Evaluation
        R1sqdd = (1 / (1 + alpha)) * M - h * delta * (QAsqd + Ssqd) + ((h^2) * beta) * ((1 / (1 + alpha)) * M2 + P4 - QAsq - Ssq)
        R2sqdd = Phiq
        J = vcat(hcat(R1sqdd, Phiq'), hcat(R2sqdd, zeros(nc, nc)))
        JCond = cond(J)

        # Solve Discretized Equations
        i = 1  # Set solution iteration counter
        err = intol + 1
        Pqdd = hcat(I, zeros(ngc, nc))
        PLam = hcat(zeros(nc, ngc), I)


        while err > intol
            # Residual Calculation
            Phiq = mathfunction.PhiqEval(tn, q, SJDT, par)
            M = mathfunction.MEval(q, SMDT, par)
            QA = mathfunction.QAEval(tn, q, qd, SMDT, STSDAT, par)
            S = mathfunction.SEval(q, qd, SMDT, par)
            R = vcat((1 / (1 + alpha)) * M * qdd + Phiq' * Lam - QA - S + R1pnm,
                (1 / (beta * h^2)) * mathfunction.PhiEval(tn, q, SJDT, par))

            if i == 1
                R1n = norm(R)
            end

            # Newton Correction
            x = -J \ R
            qdd = qdd + Pqdd * x
            Lam = Lam + PLam * x

            # Evaluate q and qd
            qd = qdnm + h * ((1 - delta) * qddnm + delta * qdd)
            q = qnm + h * qdnm + ((h^2) / 2) * ((1 - 2 * beta) * qddnm + 2 * beta * qdd)

            err = norm(R)
            i = i + 1
        end
        Jiter = i - 1

        # Higher order estimate of solution
        qe = qnm + h * qdnm + 0.5 * (h^2) * qddnm + ((h^3) / 6) * qdddnm
        qde = qdnm + h * qddnm + 0.5 * (h^2) * qdddnm

        # Evaluate error with sciq = Atol * (1 + abs(q[i])) and sciqd = Atol * (1 + abs(qd[i])),
        qdiff = q - qe
        qddiff = qd - qde
        Er = 0
        for i in 1:ngc
            sciq = Atol * (1 + abs(q[i]))
            sciqd = Atol * (1 + abs(qd[i]))
            Er += (qdiff[i] / sciq)^2  # Position error component
            Er += (qddiff[i] / sciqd)^2  # Velocity error component
        end
        Err = sqrt(Er / (2 * ngc))  # Factor of 2 removed in position error control

        hopt = h * (1 / Err)^(1 / 3)

        if hvar == 1  # variable step
            # Change step size
            if hopt < h
                h = 0.9 * hopt
                nch = n
            end

            if hopt > h
                if n > nch + 5
                    h = min(2 * h, 0.9 * hopt)
                    nch = n
                end
            end

            if h > hmax
                h = hmax
            end

            if h < 0.5 * 1e-4
                h = 0.5 * 1e-4
                Err = 0.5
            end
        end

        if hvar == 2
            Err = 0.5
        end
    end

    return q, qd, qdd, Lam, R1n, Jiter, JCond, h, nch, Err, hopt
end