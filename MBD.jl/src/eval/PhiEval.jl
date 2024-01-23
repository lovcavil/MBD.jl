export PhiEval
function PhiEval(tn, q, SJDT, par)
    nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt=parPart(par)
    Phi = zeros(nc)
    I3 = Matrix{Float64}(I, 3, 3)
    k = 1        # Joint No.
    m = 0        # Constraint Counter - 1
    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)
            Phik = bbPhidist(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end

        # Spherical Constraint
        if SJDT[1, k] == 2
            i, j, s1pr, s2pr = SphPart(k, SJDT)
            Phik = bbPhisph(i, j, s1pr, s2pr,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 3
        end

        # Cylindrical Constraint
        if SJDT[1, k] == 3
            i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm = CylPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr
            Phicyl1 = bbPhidot2(i,j,vx2pr,s1pr,s2pr,tn,q,par)
            Phicyl2 = bbPhidot2(i,j,vy2pr,s1pr,s2pr,tn,q,par)
            Phicyl3 = bbPhidot1(i,j,vz1pr,vx2pr,tn,q,par)
            Phicyl4 = bbPhidot1(i,j,vz1pr,vy2pr,tn,q,par)
            Phik = vcat(Phicyl1, Phicyl2, Phicyl3, Phicyl4)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 4
        end
        # Revolute Constraint
        if SJDT[1, k] == 4
            i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm = RevPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            Phicyl1 = bbPhidot2(i, j, vx2pr, s1pr, s2pr,tn, q, par)
            Phicyl2 = bbPhidot2(i, j, vy2pr, s1pr, s2pr,tn, q, par)
            Phicyl3 = bbPhidot1(i, j, vz1pr, vx2pr,tn, q, par)
            Phicyl4 = bbPhidot1(i, j, vz1pr, vy2pr,tn, q, par)
            Phirev5 = bbPhidot2(i, j, vz2pr, s1pr, s2pr,tn, q, par)

            Phik = vcat(Phicyl1, Phicyl2, Phicyl3, Phicyl4, Phirev5)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 5
        end
        # Check if it's a translational constraint
        if SJDT[1, k] == 5
            i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,b,mus,mud,ms,nm = TranPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            Phicyl1=bbPhidot2(i,j,vx2pr,s1pr,s2pr,tn,q,par)
            Phicyl2=bbPhidot2(i,j,vy2pr,s1pr,s2pr,tn,q,par)
            Phicyl3=bbPhidot1(i,j,vz1pr,vx2pr,tn,q,par)
            Phicyl4=bbPhidot1(i,j,vz1pr,vy2pr,tn,q,par)
            Phitran5=bbPhidot1(i,j,vy1pr,vx2pr,tn,q,par)
            Phik = vcat(Phicyl1, Phicyl2, Phicyl3, Phicyl4, Phitran5)
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
        # fxc Constraint
        if SJDT[1, k] == 1010
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)
            Phik = bbPhifxc(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end
        if SJDT[1, k] == 1020
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)
            Phik = bbPhify(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end
        if SJDT[1, k] == 1030
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)
            Phik = bbPhifz(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end

        if SJDT[1, k] == 1050
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)
            Phik = bbPhif(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end
        if SJDT[1, k] == 1060
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)
            Phik = bbPhi_Poly(i, j, s1pr, s2pr, d,tn, q, par)
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
    #println("Phi",Phi)
    return Phi
end