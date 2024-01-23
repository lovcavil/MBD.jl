export P2Eval
function P2Eval(tn,q,qd,SJDT,par)
    # Retrieve parameters from 'par'
    nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt=parPart(par)
    x=qd
    # Initialize P2 as a zero matrix of size ngc x ngc
    P2 = zeros(nc, ngc)

    # Initialize joint number and constraint counter
    k = 1
    m = 0

    # Loop over each joint/constraint
    while k <= nh
        constraintType = SJDT[1, k]
        # Process according to constraint type
        if constraintType == 1
            # Check if the constraint type is a Distance Constraint
            # Extract parameters for the Distance Constraint
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            # Compute P21 and P22 for the Distance Constraint
            P21, P22 = bbP2dist(i, j, s1pr, s2pr, d,tn, q, qd, par)

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
            P21, P22 = bbP2sph(i, j, s1pr, s2pr,tn, q, qd, par)

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
            i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm = CylPart(k, SJDT)
            uy1pr=atil(vz1pr)*vx1pr
            uy2pr=atil(vz2pr)*vx2pr
            P2cyl11,P2cyl12=bbP2dot2(i,j,vx2pr,s1pr,s2pr,tn,q,x,par);
            P2cyl21,P2cyl22=bbP2dot2(i,j,uy2pr,s1pr,s2pr,tn,q,x,par);
            P2cyl31,P2cyl32=bbP2dot1(i,j,vz1pr,vx2pr,tn,q,x,par);
            P2cyl41,P2cyl42=bbP2dot1(i,j,vz1pr,uy2pr,tn,q,x,par);

            P21k = vcat(P2cyl11, P2cyl21, P2cyl31, P2cyl41)
            P22k = vcat(P2cyl12, P2cyl22, P2cyl32, P2cyl42)

            P2 = add_constraint!(P2, P21k, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22k, m, 7 * (j - 1))
            end
            m += 4
        # Revolute Constraint
        elseif constraintType == 4
            i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,R,mus,mud,ms,nm = RevPart(k, SJDT)
            uy1pr=atil(vz1pr)*vx1pr 
            uy2pr=atil(vz2pr)*vx2pr
            
            P2cyl11,P2cyl12=bbP2dot2(i,j,vx2pr,s1pr,s2pr,tn,q,x,par)
            P2cyl21,P2cyl22=bbP2dot2(i,j,uy2pr,s1pr,s2pr,tn,q,x,par)
            P2cyl31,P2cyl32=bbP2dot1(i,j,vz1pr,vx2pr,tn,q,x,par)
            P2cyl41,P2cyl42=bbP2dot1(i,j,vz1pr,uy2pr,tn,q,x,par)
            P2rev51,P2rev52=bbP2dot2(i,j,vz2pr,s1pr,s2pr,tn,q,x,par)

            P21k = vcat(P2cyl11, P2cyl21, P2cyl31, P2cyl41, P2rev51)
            P22k = vcat(P2cyl12, P2cyl22, P2cyl32, P2cyl42, P2rev52)

            P2 = add_constraint!(P2, P21k, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22k, m, 7 * (j - 1))
            end
            m += 5

        # Translational Constraint
        elseif constraintType == 5
            i,j,s1pr,s2pr,vx1pr,vz1pr,vx2pr,vz2pr,a,b,mus,mud,ms,nm = TranPart(k, SJDT)
            uy1pr=atil(vz1pr)*vx1pr 
            uy2pr=atil(vz2pr)*vx2pr
            
            P2cyl11,P2cyl12=bbP2dot2(i,j,vx2pr,s1pr,s2pr,tn,q,x,par)
            P2cyl21,P2cyl22=bbP2dot2(i,j,uy2pr,s1pr,s2pr,tn,q,x,par)
            P2cyl31,P2cyl32=bbP2dot1(i,j,vz1pr,vx2pr,tn,q,x,par)
            P2cyl41,P2cyl42=bbP2dot1(i,j,vz1pr,uy2pr,tn,q,x,par)
            P2tran51,P2tran52=bbP2dot1(i,j,uy1pr,vx2pr,tn,q,x,par)

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
        elseif constraintType == 1010  
            # Check if the constraint type is a Distance Constraint
            # Extract parameters for the Distance Constraint
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            # Compute P21 and P22 for the Distance Constraint
            P21, P22 = bbP2fxc(i, j, s1pr, s2pr, d,tn, q, qd, par)

            # Add P21 to P2 at the appropriate location
            P2 = add_constraint!(P2, P21, m, 7 * (i - 1))

            # If j is not zero, add P22 as well
            if j >= 1
                P2 = add_constraint!(P2, P22, m, 7 * (j - 1))
            end

            # Increment the constraint counter
            m += 1
        elseif constraintType == 1020  
            # Check if the constraint type is a Distance Constraint
            # Extract parameters for the Distance Constraint
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            # Compute P21 and P22 for the Distance Constraint
            P21, P22 = bbP2fy(i, j, s1pr, s2pr, d,tn, q, qd, par)

            # Add P21 to P2 at the appropriate location
            P2 = add_constraint!(P2, P21, m, 7 * (i - 1))

            # If j is not zero, add P22 as well
            if j >= 1
                P2 = add_constraint!(P2, P22, m, 7 * (j - 1))
            end

            # Increment the constraint counter
            m += 1
        elseif constraintType == 1030  
            # Check if the constraint type is a Distance Constraint
            # Extract parameters for the Distance Constraint
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            # Compute P21 and P22 for the Distance Constraint
            P21, P22 = bbP2fz(i, j, s1pr, s2pr, d,tn, q, qd, par)

            # Add P21 to P2 at the appropriate location
            P2 = add_constraint!(P2, P21, m, 7 * (i - 1))

            # If j is not zero, add P22 as well
            if j >= 1
                P2 = add_constraint!(P2, P22, m, 7 * (j - 1))
            end

            # Increment the constraint counter
            m += 1   
        elseif constraintType == 1050  
            # Check if the constraint type is a Distance Constraint
            # Extract parameters for the Distance Constraint
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            # Compute P21 and P22 for the Distance Constraint
            P21, P22 = bbP2f(i, j, s1pr, s2pr, d,tn, q, qd, par)

            # Add P21 to P2 at the appropriate location
            P2 = add_constraint!(P2, P21, m, 7 * (i - 1))

            # If j is not zero, add P22 as well
            if j >= 1
                P2 = add_constraint!(P2, P22, m, 7 * (j - 1))
            end

            # Increment the constraint counter
            m += 1      
        elseif constraintType == 1060  
            # Check if the constraint type is a Distance Constraint
            # Extract parameters for the Distance Constraint
            i, j, s1pr, s2pr, d,ms,nm = DistPart(k, SJDT)

            # Compute P21 and P22 for the Distance Constraint
            P21, P22 = bbP2_Poly(i, j, s1pr, s2pr, d,tn, q, qd, par)

            # Add P21 to P2 at the appropriate location
            P2 = add_constraint!(P2, P21, m, 7 * (i - 1))

            # If j is not zero, add P22 as well
            if j >= 1
                P2 = add_constraint!(P2, P22, m, 7 * (j - 1))
            end

            # Increment the constraint counter
            m += 1               
        # Check if the constraint type is a Spherical Constraint
          
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