export PhiqEval
function PhiqEval(tn, q, SJDT, par)

    nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt=parPart(par)
    Phiq = zeros(Float64,nc, ngc)

    I3 = Matrix{Float64}(I, 3, 3)
    k = 1  # Joint number
    m = 0  # Constraint equation counter - 1

    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqdist(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        # Spherical Constraint
        if SJDT[1, k] == 2
            i, j, s1pr, s2pr = SphPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqsph(i, j, s1pr, s2pr, tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7 * (i - 1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7 * (j - 1))
            end

            m += 3
        end

        # Revolute Constraint
        if SJDT[1, k] == 4
            i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm = RevPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr
            Phiqcyl11,Phiqcyl12=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
            Phiqcyl21,Phiqcyl22=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
            Phiqcyl31,Phiqcyl32=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
            Phiqcyl41,Phiqcyl42=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
            Phiqrev51,Phiqrev52=bbPhiqdot2(i,j,vz2pr,s1pr,s2pr,tn,q,par);
        
            Phiq1k = vcat(Phiqcyl11, Phiqcyl21, Phiqcyl31, Phiqcyl41, Phiqrev51)
            Phiq2k = vcat(Phiqcyl12, Phiqcyl22, Phiqcyl32, Phiqcyl42, Phiqrev52)
        
            Phiq = add_constraint!(Phiq, Phiq1k, m, 7 * (i - 1))
        
            if j ≥ 1
                Phiq = add_constraint!(Phiq, Phiq2k, m, 7 * (j - 1))
            end
        
            m += 5
        end

        if SJDT[1, k] == 5
            i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,b,mus,mud,ms,nm = TranPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr
        
            Phiqcyl11,Phiqcyl12=bbPhiqdot2(i,j,vx2pr,s1pr,s2pr,tn,q,par);
            Phiqcyl21,Phiqcyl22=bbPhiqdot2(i,j,vy2pr,s1pr,s2pr,tn,q,par);
            Phiqcyl31,Phiqcyl32=bbPhiqdot1(i,j,vz1pr,vx2pr,tn,q,par);
            Phiqcyl41,Phiqcyl42=bbPhiqdot1(i,j,vz1pr,vy2pr,tn,q,par);
            Phiqtran51,Phiqtran52=bbPhiqdot1(i,j,vy1pr,vx2pr,tn,q,par);
        
            Phiq1k = vcat(Phiqcyl11, Phiqcyl21, Phiqcyl31, Phiqcyl41, Phiqtran51)
            Phiq2k = vcat(Phiqcyl12, Phiqcyl22, Phiqcyl32, Phiqcyl42, Phiqtran52)
        
            Phiq = add_constraint!(Phiq, Phiq1k, m, 7 * (i - 1))
        
            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2k, m, 7 * (j - 1))
            end
        
            m += 5
        end
        

        # Rotation Driver
        if SJDT[1, k] == 10
            i, j, vx1pr, vy1pr, vx2pr = RotDrPart(k, SJDT)
            Phiq1, Phiq2 = bbPhiqRotDr(i, j, vx1pr, vy1pr, vx2pr,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7 * (i - 1))
        
            if j ≥ 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7 * (j - 1))
            end
            m += 1
        end

        # ... (Continue with the rest of the constraints as per the given code)
        # fxc Constraint
        if SJDT[1, k] == 1010
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqfxc(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        if SJDT[1, k] == 1020
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqfy(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        if SJDT[1, k] == 1030
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqfz(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end

        if SJDT[1, k] == 1050
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqf(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        if SJDT[1, k] == 1060
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiq_Poly(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
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

    # Convert to DataFrame
    df = DataFrame(Phiq, :auto)
    # Write to CSV
    #CSV.write("Phiq.csv", df)

    return Phiq
end