function ExplicitInd0RKFN45(n, tn, Vv, Vvd, Vvdd, LLam, Uu, U, V, B, q0, h, hmax, par, SMDT, STSDAT, SJDT, nch)
    # Parameter Extraction
    nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)

    # RKFN45 Coefficients
    a1 = [0 0 0 0 0 0; 1/4 0 0 0 0 0; 3/32 9/32 0 0 0 0; 1932/2197 -7200/2197 7296/2197 0 0 0; 439/216 -8 3680/513 -845/4104 0 0; -8/27 2 -3544/2565 1859/4104 -11/40 0]
    A = a1 * a1
    c = [0, 1 / 4, 3 / 8, 12 / 13, 1, 1 / 2]
    b1 = [25 / 216, 0, 1408 / 2565, 2197 / 4104, -1 / 5, 0]'
    d1 = [16 / 135, 0, 6656 / 12825, 28561 / 56430, -9 / 50, 2 / 55]'
    B1 = b1 * a1
    D1 = d1 * a1
    v=[] 
    vd=[]
    vdd=[]
    Lam=[]
    ECond=0
    # Initialize variables
    tnm = tn - h
    vnm = Vv[:, n-1]
    vdnm = Vvd[:, n-1]
    vddnm = Vvdd[:, n-1]
    Lamnm = LLam[:, n-1]
    u = Uu[:, n-1]

    Err = 2  # Criteria for accepting time step

    # Main calculation loop
    while Err > 1
        # Calculate stage variables
        # First RKFN45 step calculation
        k1, Lam, jodeiter, ECond = FrODEfunct(n, tnm, vnm, vdnm, vddnm, Lamnm, u, q0, V, U, B, SJDT, SMDT, STSDAT, par)
        v2p = vnm + c[2] * h * vdnm
        vd2p = vdnm + h * a1[2, 1] * k1

        # Second RKFN45 step calculation
        k2, Lam, jodeiter, ECond = FrODEfunct(n, tnm + c[2] * h, v2p, vd2p, k1, Lam, u, q0, V, U, B, SJDT, SMDT, STSDAT, par)
        v3p = vnm + c[3] * h * vdnm + (h^2) * A[3, 1] * k1
        vd3p = vdnm + h * (a1[3, 1] * k1 + a1[3, 2] * k2)

        # Third RKFN45 step calculation
        k3, Lam, jodeiter, ECond = FrODEfunct(n, tnm + c[3] * h, v3p, vd3p, k2, Lam, u, q0, V, U, B, SJDT, SMDT, STSDAT, par)
        v4p = vnm + c[4] * h * vdnm + (h^2) * (A[4, 1] * k1 + A[4, 2] * k2)
        vd4p = vdnm + h * (a1[4, 1] * k1 + a1[4, 2] * k2 + a1[4, 3] * k3)

        # Fourth RKFN45 step calculation
        k4, Lam, jodeiter, ECond = FrODEfunct(n, tnm + c[4] * h, v4p, vd4p, k3, Lam, u, q0, V, U, B, SJDT, SMDT, STSDAT, par)
        v5p = vnm + c[5] * h * vdnm + (h^2) * (A[5, 1] * k1 + A[5, 2] * k2 + A[5, 3] * k3)
        vd5p = vdnm + h * (a1[5, 1] * k1 + a1[5, 2] * k2 + a1[5, 3] * k3 + a1[5, 4] * k4)

        # Fifth RKFN45 step calculation
        k5, Lam, jodeiter, ECond = FrODEfunct(n, tnm + c[5] * h, v5p, vd5p, k4, Lam, u, q0, V, U, B, SJDT, SMDT, STSDAT, par)
        v6p = vnm + c[6] * h * vdnm + (h^2) * (A[6, 1] * k1 + A[6, 2] * k2 + A[6, 3] * k3 + A[6, 4] * k4)
        vd6p = vdnm + h * (a1[6, 1] * k1 + a1[6, 2] * k2 + a1[6, 3] * k3 + a1[6, 4] * k4 + a1[6, 5] * k5)

        # Sixth RKFN45 step calculation
        k6, Lam, jodeiter, ECond = FrODEfunct(n, tnm + c[6] * h, v6p, vd6p, k5, Lam, u, q0, V, U, B, SJDT, SMDT, STSDAT, par)
        k = [k1 k2 k3 k4 k5 k6]

        # Evaluate Solution for v, vd, vdd
        v = vnm + h * vdnm + (h^2) * k * B1'
        vd = vdnm + h * k * b1'
        vdd, Lam, jodeiter, ECond = FrODEfunct(n, tnm + h, v, vd, k6, Lam, u, q0, V, U, B, SJDT, SMDT, STSDAT, par)

        if hvar == 1  # Variable step
            # Evaluate embedded solution for error estimate
            vnhat = vnm + h * vdnm + (h^2) * k * D1'
            vdnhat = vdnm + h * k * d1'

            # Error evaluation
            vdiff = v - vnhat
            vddiff = vd - vdnhat
            Er = 0
            for i in 1:ngc-nc
                sciv = Atol * (1 + abs(v[i]))
                scivd = Atol * (1 + abs(vd[i]))
                Er += (vdiff[i] / sciv)^2  # v-contribution
                Er += (vddiff[i] / scivd)^2  # vd-contribution
            end
            Err = sqrt(Er / (2 * (ngc - nc)))  # factor of 2 for displacement + velocity

            # Step size adjustment
            hopt = h * (1 / Err)^(1 / 6)

            if hopt < h
                h = 0.9 * hopt
                nch = n
                if h < 1e-5
                    h = 1e-5
                    Err = 0.5
                end
            elseif hopt > h && n > nch + 5
                h = min(2 * h, 0.9 * hopt)
                nch = n
            end

            if h > hmax
                h = hmax
            end
        elseif hvar == 2
            Err = 0.5
        end
    end

    return v, vd, vdd, Lam, ECond, h, nch
end