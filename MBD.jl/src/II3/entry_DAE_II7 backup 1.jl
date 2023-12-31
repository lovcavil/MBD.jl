include("mathfunction_II7.jl")
include("AppData_II7.jl")
include("ImplicitIndex1.jl")
include("ImplicitIndex3.jl")
include("ExplicitNystrom4.jl")
include("ExplicitRKFN45.jl")
using LinearAlgebra
using Plots
using CSV, DataFrames
using IterativeRefinement

function run(app,my_p)
    #println("app", app)
    # Integration and Error Control Parameters
    intol = 10^-6     # Tolerance in solving discretized equations of motion
    Atol = 10^-5      # Absolute error tolerance for variable step methods
    MaxIntiter = 8    # Limit on number of integration iterations
    MaxJcond = 200    # Limit on magnitude of Jcond
    R1nmax = 15.0       # Limit on residual in integration
    MaxECond = 10     # Limit on magnitude of ECond
    PosConstrMax = 10^-4  # Limit on position constraint error
    VelConstrMax = 10^-4  # Limit on velocity constraint error
    AccConstrMax = 10^20  # Limit on acceleration constraint error

    h0 = 0.002       # Initial time step
    h = h0
    hmax = 0.01
    hvar = 1          # hvar=1, variable h; hvar=2, constant h
  
    constrcor = 1     # constrcor=1, correct; constrcor=2, no correct
    Initialpositioncorrection = 1
    tfinal = 5.5

    # Integration Method
    integ = 5     
    # integ=1-ImplicitIndex1; alpha=0-Trapezoidal; alpha<0-HHT
    # integ=2-ImplicitIndex2; alpha=0-Trapezoidal; alpha<0-HHT
    # integ=3-ImplicitIndex3; alpha=0-Trapezoidal; alpha<0-HHT
    # integ=4-ExplicitNystrom4Index1
    # integ=5-ExplicitRKFN45Index1

    alpha = -0.3         # Enter -1/3 <= alpha < 0 for HHT
    # alpha not used in ExplicitNystrom4 and ExplicitRKFN45

    g = 9.81

    # Application Data         
    # Select application:
    # 1: Pendulum, Spherical to Ground
    # 2: Spin Stabilized Top, Spherical to Ground
    # 3: One Body Pendulum, Dist. to Ground
    # 4: Double Pendulum
    # 5: One Body Cylindrical with Spring
    # 6: Spatial Slider-Crank
    # 7: Transient Top, Spherical to Ground

    # Application Data Function
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(app)
    
    par = Any[nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA,my_p]


    # Initial condition calculation, if required

    ux = [1, 0, 0]
    uy = [0, 1, 0]
    uz = [0, 0, 1]
    zer = zeros(3)

    # Data Storage Arrays
    # Preallocate a large array with more columns than currently needed.
    # For example, 100 columns as a starting point:
    Q = zeros(ngc, 100)
    Qd = zeros(ngc, 100)
    Qdd = zeros(ngc, 100)
    LLam = zeros(nc, 100)
    ConPhi = zeros(nc, 100)
    PE = zeros(100)
    filled_columns = 1

    # Initial position correction
    z = zeros(ngc)
    nu = zeros(nc)
    Ingc = Matrix{Float64}(I, ngc, ngc)
    Pz = hcat(Ingc, zeros(ngc, nc))
    Inc = Matrix{Float64}(I, nc, nc)
    Pnu = hcat(zeros(nc, ngc), Inc)
    err = intol + 1
    i = 0
    w=[]
    if Initialpositioncorrection == 1
        println("Initial Position Correction",
                "\n--------------------------")
        while err > intol
            println("-------------i=", i)
            Phiq = mathfunction.PhiqEval(0, q + z, SJDT, par)
            Resid = vcat(z + Phiq' * nu, mathfunction.PhiEval(0, q + z, SJDT, par))
            # Concatenate the matrices
            JJ = vcat(hcat(I + mathfunction.P4Eval(0, q + z, nu, SJDT, par), Phiq'), hcat(Phiq, zeros(nc, nc)))

            # Convert to DataFrame
            df = DataFrame(JJ, :auto)

            # Write to CSV
            CSV.write("JJ.csv", df)
            # Convert to DataFrame
            #df = DataFrame(Resid, :auto)

            # Write to CSV
            #CSV.write("Resid.csv", df)
            #w = -JJ \ Resid

            # if norm(Resid)!=0
            #     #w = -JJ \ Resid
            # end
            #w, bnorm, bcomp = rfldiv(-JJ,Resid)
            #w = -pinv(JJ) * Resid

            println("rank", rank(JJ))
            try
                println("normal")
                w = -JJ \ Resid
                #w = -(JJ+0.01*I) \ Resid
            catch e
                println(e)
                println("w = -(JJ+0.0001*I) \\ Resid")
                w = -(JJ + 0.01 * I) \ Resid
            end

            z += Pz * w
            nu += Pnu * w
            err = norm(Resid)
            println("err", err)
            i += 1
        end
        q += z
    end
    
    # # Initial velocity correction
    # Phiq = mathfunction.PhiqEval(0, q, SJDT, par)
    # mu = pinv(Phiq * Phiq') * (Phiq * qd)
    # delqd = -Phiq' * mu
    # qd += delqd

    Q[:, 1] = q
    Qd[:, 1] = qd

    # Calculate Initial Acceleration and Lagrange Multipliers
    M = mathfunction.MEval(q, SMDT, par)
    Phiq = mathfunction.PhiqEval(0, q, SJDT, par)
    QA = mathfunction.QAEval(0, q, qd, SMDT, STSDAT, par)
    S = mathfunction.SEval(q, qd, SMDT, par)
    Gam = mathfunction.GamEval(0, q, qd, SJDT, par)
    # Concatenate the matrices
    EE = vcat(hcat(M, Phiq'), hcat(Phiq, zeros(nc, nc)))
    EEcond = cond(EE)
    RHS = [QA + S; -Gam]
    #x = EE \ RHS
    x = pinv(EE) * RHS
    Pqdd = hcat(I, zeros(ngc, nc))
    PLam = hcat(zeros(nc, ngc), I)
    qdd = Pqdd * x
    Lam = PLam * x

    Qdd[:, 1] = qdd
    LLam[:, 1] = Lam

    println("cor-q0",q)
    println("cor-qd0",qd)


    # Initialize Data For Integration
    n = 1
    t = [0.0]
    Jiterrpt = [0]
    Jiter = 0
    JCondrpt = [0.0]
    R1Normrpt = [R1nmax + 1]
    ECondrpt = [MaxECond + 1]
    Econd = [1.0]
    PosConstrNorm = [0.0]
    VelConstrNorm = [0.0]
    AccConstrNorm = [0.0]
    corvel = 0
    corpos = 0
    nch = 1.0
    Errrpt = [0.0]
    hrpt = [0.0]
    KE = [0.0]
    SE = [0.0]
    PE[1] = 0.0
    TE = [0.0]
    x1 = [0.0]
    y1 = [0.0]
    z1 = [0.0]
    x2 = [0.0]
    y2 = [0.0]
    z2 = [0.0]
    x3 = [0.0]
    y3 = [0.0]
    z3 = [0.0]

    xd1 = [0.0]
    yd1 = [0.0]
    zd1 = [0.0]
    xd2 = [0.0]
    yd2 = [0.0]
    zd2 = [0.0]
    xd3 = [0.0]
    yd3 = [0.0]
    zd3 = [0.0]

    corvelrpt = [0.0]
    corposrpt = [0]
    corpositer = [0]
    ECondrpt = [0.0]
    # Integration
    while t[n] < tfinal

        #global t, tn, n, filled_columns, Q, Qd, Qdd, LLam, ConPhi, h,nch,KE,PE,SE,TE
        #global q, qd, qdd, Lam, R1n, Jiter, JCond, h, nch, Err
        # When you need to add a new column:
        if filled_columns == size(Q, 2)
            # If the array is full, double its size.
            Q = hcat(Q, zeros(ngc, size(Q, 2)))
            Qd = hcat(Qd, zeros(ngc, size(Qd, 2)))
            Qdd = hcat(Qdd, zeros(ngc, size(Qdd, 2)))
            PE = hcat(PE, zeros(size(PE)))
            LLam = hcat(LLam, zeros(nc, size(LLam, 2)))
        end


        filled_columns += 1

        n += 1
        push!(t, t[n-1] + h)
        tn = t[n]
        tnm = t[n-1]

        # Integration
        if integ == 1
            q, qd, qdd, Lam, R1n, Jiter, JCond, h, nch, Err = ImplicitIndex1(n, tn, Q, Qd, Qdd, LLam, h, hmax, SMDT, STSDAT, SJDT, par, alpha, nch)

            push!(R1Normrpt, R1n)
            push!(Jiterrpt, Jiter)
            push!(JCondrpt, JCond)
            push!(Errrpt, Err)
            push!(hrpt, h)
        end

        if integ == 3
            q,qd,qdd,Lam,R1n,Jiter,JCond,h,nch,Err,hopt =ImplicitIndex3(n,tn,Q,Qd,Qdd,LLam,h,hmax,SMDT,STSDAT,SJDT,par,alpha,nch);
            
            push!(R1Normrpt, R1n)
            push!(Jiterrpt, Jiter)
            push!(JCondrpt, JCond)
            push!(Errrpt, Err)
            push!(hrpt, h)
        end

        if integ == 4
            q, qd, qdd, Lam, ECond = ExplicitNystrom4(n, tn, Q, Qd, h, SMDT, STSDAT, SJDT, par)
            push!(ECondrpt, ECond)
        end

        if integ == 5
            q, qd, qdd, ECond, h, nch = ExplicitRKFN45(n, tn, Q, Qd, h, hmax, par, SMDT, STSDAT, SJDT, nch)
            #println("q",q)
            push!(ECondrpt, ECond)
            push!(hrpt, h)
        end

        # Corrections if velocity or position errors exceed tolerances
        if constrcor == 1

            Phiq = mathfunction.PhiqEval(tn, q, SJDT, par)
            if norm(Phiq * qd) > VelConstrMax
                mu = pinv(Phiq * Phiq') * (Phiq * qd)
                delqd = -Phiq' * mu
                qd += delqd
                corvel += 1
                push!(corvelrpt, corvel)
                
            end
            # println("1norm(Phiq * qd)",norm(Phiq * qd))
            if norm(mathfunction.PhiEval(tn, q, SJDT, par)) > PosConstrMax
                z = zeros(ngc)
                nu = zeros(nc)
                Pz = hcat(I(ngc), zeros(ngc, nc))
                Pnu = hcat(zeros(nc, ngc), I(nc))
                err = intol + 1
                i = 0
                while err > intol
                    Phiq = mathfunction.PhiqEval(tn, q + z, SJDT, par)
                    Resid = vcat(z + Phiq' * nu, mathfunction.PhiEval(tn, q + z, SJDT, par))
                    JJ = vcat(hcat(I + mathfunction.P4Eval(0, q + z, nu, SJDT, par), Phiq'), hcat(Phiq, zeros(nc, nc)))
                    w = -inv(JJ) * Resid
                    z += 0.5*Pz * w
                    nu += 0.5*Pnu * w
                    err = norm(Resid)
                    i += 1
                end
                q += z
                corpos += 1
                push!(corposrpt, corpos)
                push!(corpositer, i)
            end
            if norm(Phiq * qd) > VelConstrMax
                mu = pinv(Phiq * Phiq') * (Phiq * qd)
                delqd = -Phiq' * mu
                qd += delqd
                corvel += 1
                push!(corvelrpt, corvel)

                
            end
            Phiq = mathfunction.PhiqEval(tn, q, SJDT, par)
            # println("2norm(Phiq * qd)",norm(Phiq * qd))
        end
        # End Corrections

        # Record Solution
        Q[:, n] = q
        Qd[:, n] = qd
        Qdd[:, n] = qdd
        LLam[:, n] = Lam

        # Calculate Total Energy
        M = mathfunction.MEval(q, SMDT, par)
        push!(KE, 0.5 * qd' * M * qd)

        i = 1
        PE[n] = 0
        while i <= nb
            PE[n] += SMDT[1, i] * g * q[7*(i-1)+3]
            i += 1
        end

        push!(SE, 0)
        T = 1
        while T <= NTSDA
            i, j, s1pr, s2pr, K, C, el0, F = STSDATPart(STSDAT, T)
            r1, p1 = qPart(q, i)
            r2 = [0, 0, 0]
            p2 = [1, 0, 0, 0]
            r1d = [0, 0, 0]
            p1d = zeros(4)
            if j >= 1
                r2, p2 = qPart(q, j)
            end
            A1 = ATran(p1)
            A2 = ATran(p2)
            d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
            el = sqrt(d12' * d12)
            SET = 0.5 * K * (el - el0)^2
            SE[n] += SET
            T += 1
        end
        push!(TE, 0.0)
        TE[n] = KE[n] + PE[n] + SE[n]

        # Data of Interest (Enter for each application)
        if app == 1 || app == 101 || app == 102 || app == 103 || app == 201|| app == 204|| app == 205|| app == 206|| app == 207
            push!(x1, q[1])
            push!(y1, q[2])
            push!(z1, q[3])
        end
        if app == 2
            push!(x1, q[1])
            push!(y1, q[2])
            push!(z1, q[3])
        end
        if app == 6
            push!(x1, q[1])
            push!(y1, q[2])
            push!(z1, q[3])
            push!(x2, q[8])
            push!(y2, q[9])
            push!(z2, q[10])

        end
        if app == 202 || app == 203
            push!(x1, q[1])
            push!(y1, q[2])
            push!(z1, q[3])
            push!(x2, q[8])
            push!(y2, q[9])
            push!(z2, q[10])
        end
        target=[]
        if app == 301
            temp=[0.0]
            # Using an array of arrays
            target = [x1, y1, z1,temp,temp,temp,temp, x2, y2, z2,temp,temp,temp,temp]
            target = [x1, y1, z1]
            flags = [true, true, true,false,false,false,false,true, true, true,false,false,false,false]
            flags = [true, true, true,false,false,false,false]
            # Apply the function
            process_vector(q, flags, target, append_to_vector)
        end
        if app == 302||app==303||app==304||app==305
            temp=[0.0]
            # Using an array of arrays
            target = [x1, y1, z1,temp,temp,temp,temp, x2, y2, z2,temp,temp,temp,temp, x3, y3, z3,temp,temp,temp,temp]

            flags = [true, true, true,false,false,false,false,true, true, true,false,false,false,false,true, true, true,false,false,false,false]

            # Apply the function
            process_vector(q, flags, target, append_to_vector)
           
            target = [xd1, yd1, zd1,temp,temp,temp,temp,     xd2, yd2, zd2,temp,temp,temp,temp,       xd3, yd3, zd3,temp,temp,temp,temp]

            flags = [true, true, true,false,false,false,false,true, true, true,false,false,false,false,true, true, true,false,false,false,false]

            # Apply the function
            process_vector(qd, flags, target, append_to_vector)
        end

        # Calculate constraint error
        Phi = mathfunction.PhiEval(tn, q, SJDT, par)
        Phiq = mathfunction.PhiqEval(tn, q, SJDT, par)
        Gam = mathfunction.GamEval(tn, q, qd, SJDT, par)
        #println("3norm(Phiq * qd)",norm(Phiq * qd))
        push!(PosConstrNorm, norm(Phi))
        push!(VelConstrNorm, norm(Phiq * qd))
        push!(AccConstrNorm, norm(Phiq * qdd + Gam))

    end
    Plots.default(show = true)
    if app == 1 || app == 101 || app == 102 || app == 103 || app == 201|| app == 204|| app == 205|| app == 206|| app == 207
        #f = plot3d(x1, y1, z1, xlabel="X-axis", ylabel="Y-axis", zlabel="Z-axis", title="3D Plot")
        #display(f)
        f0 = plot(t, x1, label="x", title="Output Data Plot", xlabel="t", ylabel="x", legend=:topright)
        display(f0)
        f1 = plot(t, y1, label="y", title="Output Data Plot", xlabel="t", ylabel="y", legend=:topright)
        display(f1)
        f2 = plot(t, z1, label="z", title="Output Data Plot", xlabel="t", ylabel="z", legend=:topright)
        display(f2)
        f = plot3d(x1, y1, z1, xlabel="X-axis", ylabel="Y-axis", zlabel="Z-axis", title="3D Plot")
        display(f)
        f3 = plot(t, TE, label="TE", title="Output Data Plot", xlabel="t", ylabel="TE", legend=:topright)
        display(f3)
        f4 = plot(t, TE, label="TE", title="Output Data Plot", xlabel="t", ylabel="TE", legend=:topright)
        display(f4)
        f5 = plot(t, PosConstrNorm, label="PosConstrNorm", title="Output Data Plot", xlabel="t", ylabel="PosConstrNorm", legend=:topright)
        display(f5)
        f6 = plot(t, VelConstrNorm, label="VelConstrNorm", title="Output Data Plot", xlabel="t", ylabel="VelConstrNorm", legend=:topright)
        display(f6)
        f7 = plot(t, AccConstrNorm, label="AccConstrNorm", title="Output Data Plot", xlabel="t", ylabel="AccConstrNorm", legend=:topright)
        display(f7)
        #f8 = plot(t, Errrpt, label="Errrpt", title="Output Data Plot", xlabel="t", ylabel="Errrpt", legend=:topright)
        #display(f8)
        #f9 = plot(t, hrpt, label="hrpt", title="Output Data Plot", xlabel="t", ylabel="hrpt", legend=:topright)
        #display(f9)
        println(corvelrpt)
        println(corposrpt)
        println(corpositer)

        #f41 = plot(t, corvelrpt, label="corvelrpt", title="Output Data Plot", xlabel="t", ylabel="corvelrpt", legend=:topright)
        #display(f41)
        #f42 = plot(t, corposrpt, label="corposrpt", title="Output Data Plot", xlabel="t", ylabel="corposrpt", legend=:topright)
        #display(f42)
        #f43 = plot(t, corpositer, label="corpositer", title="Output Data Plot", xlabel="t", ylabel="corpositer", legend=:topright)
        #display(f43)
    end
    if app == 6
        f1 = plot(t, z2, label="z2", title="Output Data Plot", xlabel="t", ylabel="z", legend=:topright)
        display(f1)

        folder_path = "demo6"
        if !isdir(folder_path)
            mkdir(folder_path)
        end
        #plot()  # Initialize an empty plot
        target = [(x1, "x1"), (y1, "y1"), (z1, "z1"), (x2, "x2"), (y2, "y2"), (z2, "z2"),
            (PosConstrNorm, "PosConstrNorm"), (VelConstrNorm, "VelConstrNorm"), (AccConstrNorm, "AccConstrNorm")]

        df = DataFrame(t=t, x1=x1, y1=y1, z1=z1, x2=x2, y2=y2, z2=z2, PosConstrNorm=PosConstrNorm,
            VelConstrNorm=VelConstrNorm, AccConstrNorm=AccConstrNorm)
        filename = "data.csv"
        CSV.write(joinpath(folder_path, filename), df)
    end
    if app == 205

        folder_path = "demo205"
        if !isdir(folder_path)
            mkdir(folder_path)
        end
        #plot()  # Initialize an empty plot
        target = [(x1, "x1"), (y1, "y1"), (z1, "z1"), 
            (PosConstrNorm, "PosConstrNorm"), (VelConstrNorm, "VelConstrNorm"), (AccConstrNorm, "AccConstrNorm")]

        df = DataFrame(t=t, x1=x1, y1=y1, z1=z1,  PosConstrNorm=PosConstrNorm,
            VelConstrNorm=VelConstrNorm, AccConstrNorm=AccConstrNorm)
        filename = "data.csv"
        CSV.write(joinpath(folder_path, filename), df)
    end
    if app == 202 || app == 203
        f5 = plot(t, x1, label="x", title="Output Data Plot", xlabel="t", ylabel="x", legend=:topright)
        display(f5)
        f6 = plot(t, y1, label="y", title="Output Data Plot", xlabel="t", ylabel="y", legend=:topright)
        display(f6)
        f7 = plot(t, z1, label="z", title="Output Data Plot", xlabel="t", ylabel="z", legend=:topright)
        display(f7)
        f = plot3d(x1, y1, z1, xlabel="X-axis", ylabel="Y-axis", zlabel="Z-axis", title="3D Plot")
        display(f)
        f2 = plot3d(x2, y2, z2, xlabel="X-axis", ylabel="Y-axis", zlabel="Z-axis", title="3D Plot2")
        display(f2)
    end
    
    if app==301
        # Sample dictionary
        my_dict = par[11]

        # Initialize an empty string
        concatenated_string = ""

        # Iterate over the values in the dictionary
        for value in values(my_dict)
            # Concatenate each value to the string
            concatenated_string *= string(value)
        end

        # Display the concatenated string
        println(concatenated_string)  # Output will be "12.02.0"
        folder_path = concatenated_string
        mkdir(folder_path)
        
        #plot()  # Initialize an empty plot
        target = [(x1, "x1"),(y1,"y1"),(z1,"z1"),(x2, "x2"),(y2,"y2"),(z2,"z2"),(PosConstrNorm,"PosConstrNorm"),(VelConstrNorm,"VelConstrNorm"),(AccConstrNorm,"AccConstrNorm")]
        target = [(x1, "x1"),(y1,"y1"),(z1,"z1"),(PosConstrNorm,"PosConstrNorm"),(VelConstrNorm,"VelConstrNorm"),(AccConstrNorm,"AccConstrNorm")]
        pyplot()
        for (vec, label) in target
            #f=plot(t, vec, label=label, title="Output Data Plot", xlabel="t", ylabel=label, legend=:topright)  # Add each vector to the plot
            #display(f)
            filename = "$label.png"
            #savefig(f,joinpath(folder_path, filename))
        end
        #display(plot())  # Display the plot
        #f = plot3d(x1, y1, z1, xlabel="X-axis", ylabel="Y-axis", zlabel="Z-axis", title="3D Plot")
        #display(f)
        filename = "3d.png"
        #savefig(f,joinpath(folder_path, filename))



        df = DataFrame(t=t,x1= x1,y1=y1,z1=z1,PosConstrNorm=PosConstrNorm,
        VelConstrNorm=VelConstrNorm,AccConstrNorm=AccConstrNorm)
        filename="data.csv"

        CSV.write(joinpath(folder_path, filename), df)
        


    end

    if app==302||app==303||app==304||app==305
        # Sample dictionary
        my_dict = par[11]

        # Initialize an empty string
        concatenated_string = ""

        # Iterate over the values in the dictionary
        for value in values(my_dict)
            # Concatenate each value to the string
            concatenated_string *= string(value)
        end

        # Display the concatenated string
        println(concatenated_string)  # Output will be "12.02.0"
        folder_path = concatenated_string
        if !isdir(folder_path)
            mkdir(folder_path)
        end
        
        #plot()  # Initialize an empty plot
        target = [(x1, "x1"),(y1,"y1"),(z1,"z1"),(x2, "x2"),(y2,"y2"),(z2,"z2"),(x3, "x3"),(y3,"y3"),(z3,"z3"),
        (PosConstrNorm,"PosConstrNorm"),(VelConstrNorm,"VelConstrNorm"),(AccConstrNorm,"AccConstrNorm")]
        
        pyplot()
        for (vec, label) in target
            #f=plot(t, vec, label=label, title="Output Data Plot", xlabel="t", ylabel=label, legend=:topright)  # Add each vector to the plot
            #display(f)
            filename = "$label.png"
            #savefig(f,joinpath(folder_path, filename))
        end
        #display(plot())  # Display the plot
        #f = plot3d(x1, y1, z1, xlabel="X-axis", ylabel="Y-axis", zlabel="Z-axis", title="3D Plot")
        #display(f)
        filename = "3d.png"
        #savefig(f,joinpath(folder_path, filename))

        df = DataFrame(t=t,x1= x1,y1=y1,z1=z1,x2= x2,y2=y2,z2=z2,x3= x3,y3=y3,z3=z3,
        xd1= xd1,yd1=yd1,zd1=zd1,xd2= xd2,yd2=yd2,zd2=zd2,xd3= xd3,yd3=yd3,zd3=zd3,PosConstrNorm=PosConstrNorm,
        VelConstrNorm=VelConstrNorm,AccConstrNorm=AccConstrNorm)
        filename="data.csv"

        CSV.write(joinpath(folder_path, filename), df)
        
    end
end

#run(101)
#run(102)
#run(103)
#run(301)
my_p = Dict("p1" => 1, "p2" => 1,"p3" => 1, "p4" => 1)
#run(302, my_p)  # Call function from file_b.jl
#run(6, my_p)  # Call function from file_b.jl

#run(303, my_p)  # Call function from file_b.jl
run(305, my_p)  # Call function from file_b.jl