function Ind0IC(q0, qd0, SMDT, SJDT, STSDAT, par, N)

    nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)

    # Evaluate Initial Acceleration and Lagrange Multipliers

    # Initial Parameterization
    Phiq = PhiqEval(0, q0, SJDT, par)
    U = Phiq'
    B = inv(U' * U)
    V = nullspace(U')
    v = zeros(nv)
    vd = V' * qd0
    u = zeros(nu)
    Vv = zeros(nv, 1)
    Vvd = zeros(nv, 1)
    Uu = zeros(nu, 1)
    Vv[:, 1] = v
    Vvd[:, 1] = vd
    Uu[:, 1] = u

    # Increment friction coefficients to obtain initial conditions on Lam and vdd
    w = 0
    Lam = zeros(nc)
    vdd = zeros(nv)
    Econd0 = zeros(N + 1)
    jodeiter0 = zeros(Int, N + 1)
    Qdd0 = zeros(size(q0, 1), N + 1)
    LLam0 = zeros(size(Lam, 1), N + 1)
    qdd = zeros(size(q0, 1))
    while w <= N
        vdd, Lam, jodeiter, ECond = FrODEfunct0w(q0, qd0, vdd, Lam, V, U, B, SJDT, SMDT, STSDAT, par, w, N)

        w += 1

        Gam = GamEval(0, q0, qd0, SJDT, par)
        qdd = V * vdd - U * B * Gam
        Econd0[w] = ECond
        jodeiter0[w] = jodeiter
        Qdd0[:, w] = qdd
        LLam0[:, w] = Lam
    end

    return qdd, Lam, Qdd0, LLam0, w
end
