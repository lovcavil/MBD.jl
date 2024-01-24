include("solver_II8.jl")
include("../problem/AppData_II8.jl")
using LinearAlgebra
using DifferentialEquations
using OrdinaryDiffEq
using Sundials
using Plots
using CSV, DataFrames

function test_II1()

    h0 = 0.001       # Initial time step
    hvar = 0
    g = 9.81

    # Application Data Function
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(201)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Integration P
    t = 0.0
    tspan = (0.0, 1.0)
    println("START")
    qdd, Lam, ECond = mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀ = vcat(q, Lam, qd)
    du₀ = vcat(qd, zeros(nc), qdd)

    differential_vars = create_bool_vector(nb, nc)
    simulation_params_dict = Dict(
        :g => 9.81,
        :SMDT => SMDT,
        :STSDAT => STSDAT,
        :SJDT => SJDT,
        :par => par
    )
    println("u₀=", u₀)
    println("du₀=", du₀)
    println("differential_vars=", differential_vars)  # 输出 differential_vars 的值

    prob = DAEProblem(implicitdae, du₀, u₀, tspan, simulation_params_dict, differential_vars=differential_vars)

    # Integration
    default_solve_kwargs = Dict(:alg => IDA(), :reltol => 1e-6, :abstol => 1e-6, :progress => true, :dtmin => 0.00001)
    sol = solve(prob; default_solve_kwargs...)

    return sol
end

function test_II3()

    h0 = 0.001       # Initial time step
    hvar = 0
    g = 9.81

    # Application Data Function
    nb, ngc, nh, nc, NTSDA, SJDT, SMDT, STSDAT, q, qd = AppData_II7(201)
    par = Any[nb, ngc, nh, nc, g, 0, 0, h0, hvar, NTSDA]

    # Integration P
    t = 0.0
    tspan = (0.0, 1.0)
    println("START")
    qdd, Lam, ECond = mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    u₀ = vcat(q, Lam, qd)
    du₀ = vcat(qd, zeros(nc), qdd)

    differential_vars = create_bool_vector(nb, nc)
    simulation_params_dict = Dict(
        :g => 9.81,
        :SMDT => SMDT,
        :STSDAT => STSDAT,
        :SJDT => SJDT,
        :par => par
    )
    println("u₀=", u₀)
    println("du₀=", du₀)
    println("differential_vars=", differential_vars)  # 输出 differential_vars 的值

    prob = DAEProblem(implicit3dae, du₀, u₀, tspan, simulation_params_dict, differential_vars=differential_vars)

    # Integration
    default_solve_kwargs = Dict(:alg => IDA(), :reltol => 1e-6, :abstol => 1e-6, :progress => true, :dtmin => 0.00001)
    sol = solve(prob; default_solve_kwargs...)

    return sol
end

function test_EI0(; app::Int=6, tspan::Tuple{Float64,Float64}=(0.0, 5.0))
    # CAKD initialization
    # & Application Data Function
    h0 = 0.01       # Initial time step
    h = h0
    hvar = 0
    g = 9.8
    utol = 10^-6
    Btol = 10^-6
    intol = 10^-6
    Atol = 10^-5
    Vtol = 10^-4
    Maxv = 0.8#       %Limit on magnitude of v
    MaxIntiter = 12 #      %Limit on number of integration iterations
    MaxUiter = 8 #      %Limit on number of U iterations
    MaxJcond = 100 #   %Limit on magnitude of Econd in Trap
    R1nmax = 30000 #     %Limit on residual in integrationt
    MaxECond = 2000.0 #  %Limit on magnitude of ECond
    MaxBCond = 6000.0 #  %Limit on magnitude of BCond
    MaxBnormRat = 10.0 #  %Maximum ratio of Bnorm to that after parameterizationm
    MaxBCondRat = 10.0 #     %Maximum ratio of BCond to BCond0 at parameterization

    vt = 100 * h0
    integ = 3
    tfinal = 1.0
    nb, ngc, nh, nc, nv, nu, NTSDA, SJDT, SMDT, STSDAT, q0, qd0 = AppData_II8(app)
    par = Any[nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt]

    # Initialization
    println("START of Explicit Index 0")

    # Correction
    # t=tspan[1]
    # qdd, Lam, ECond=mathfunction.ODEfunct(t, q, qd, SMDT, STSDAT, SJDT, par)
    # u₀  = vcat(q,Lam,qd)
    # du₀ = vcat(qd,zeros(nc),qdd)

    # p=Any[SMDT,STSDAT,SJDT,par]

    # Initial condition calculation, if required

    ux = [1; 0; 0]
    uy = [0; 1; 0]
    uz = [0; 0; 1]
    zer = zeros(3)

    Q = zeros(ngc, 10)
    Qd = zeros(ngc, 10)
    Qdd = zeros(ngc, 10)
    Vv = zeros(nv, 10)
    Vvd = zeros(nv, 10)
    Vvdd = zeros(nv, 10)
    Uu = zeros(nu, 10)
    Uud = zeros(nu, 10)
    Uudd = zeros(nu, 10)
    LLam = zeros(nc, 10)
    LLam0 = zeros(nc, 10)
    QQA = zeros(ngc,10)

    Q[:, 1] = q0
    Qd[:, 1] = qd0

    # Calculation of initial values of qdd and Lam to start integration
    # Coefficients of friction are indexed by w/N, w=0,1,...N in Function 
    # Ind0IC to calculate a sequence of qdd and Lam that converge to the
    # desired coefficients of friction. The user may select N.

    N = 10
    qdd, Lam, Qdd0, LLam0, w = Ind0IC(q0, qd0, SMDT, SJDT, STSDAT, par, N)

    Qdd[:, 1] = qdd
    LLam[:, 1] = Lam

    # Initialize Data For Integration

    n = 1
    t = zeros(10)
    vnorm = [Maxv + 1]  # Assuming Maxv is defined
    jRepar = 0  # Counter for Reparameterization
    Cr = 2
    Jiterrpt = [0]
    Jiter = 0
    Uterurpt = [0]
    JCondrpt = [0]
    R1Normrpt = [R1nmax + 1]  # Assuming R1nmax is defined
    vnormrpt = [Maxv + 1]
    vdnormrpt = [0]
    vddnormrpt = [0]
    ECondrpt = [MaxECond + 1]  # Assuming MaxECond is defined
    Econd = [1]
    uiterrpt = [1]
    nch = 1
    VelConstrNorm = [0]
    AccConstrNorm = [0]
    nres = 1
    jnres = 0
    B0norm = 1
    Bnorm=[0]
    BnormRat = [1]
    BCond = [MaxBCond + 1]  # Assuming MaxBCond is defined
    Bcond = [MaxBCond + 1] 
    BCondRat = [1]
    Jinv = I(nv)  # Assuming nv is defined
    Crrpt = [0]
    jReparrpt = [0]
    jodeiterRpt = [0]
    Iterurpt=[0]
    hrpt = [0]
    unorm=[0]
    qnormrpt=[0]
    Biterrpt=[0]
    Binvnorm=[0]
    qdnormrpt=[0]
    qddnormrpt=[0]
    TE=[0]
    PE=[0]
    SE=[0]
    KE=[0]
    # Integration
    filled_columns=1
    V=[]
    U=[]
    B=[]
    npar=0
    BCond0=0
    # app=4 initial
    omegaz1 = [0.0]
    theta1 = [0.0]
    dely2pr = [0.0]

    while t[n] < tfinal

        if filled_columns == size(Q, 2)
            #println("double")
            # If the array is full, double its size.
            t = vcat(t, zeros(size(t)))
            Q = hcat(Q, zeros(ngc, size(Q, 2)))
            Qd = hcat(Qd, zeros(ngc, size(Qd, 2)))
            Qdd = hcat(Qdd, zeros(ngc, size(Qdd, 2)))

            LLam = hcat(LLam, zeros(nc, size(LLam, 2)))
            LLam0 = hcat(LLam0, zeros(nc, size(LLam0, 2)))

            Vv = hcat(Vv, zeros(nv, size(Vv, 2)))
            Vvd = hcat(Vvd, zeros(nv, size(Vvd, 2)))
            Vvdd = hcat(Vvdd, zeros(nv, size(Vvdd, 2)))

            Uu = hcat(Uu, zeros(nu, size(Uu, 2)))
            Uud = hcat(Uud, zeros(nu, size(Uud, 2)))
            Uudd = hcat(Uudd, zeros(nu, size(Uudd, 2)))
            QQA = hcat(QQA, zeros(ngc, size(QQA, 2)))

        end
        filled_columns += 1

        n += 1
        t[n] = t[n-1] + h
        tn = t[n]
        tnm = t[n-1]

        # Start implicit reparameterization criteria
        if integ < 3
            if vnormrpt[n-1] > Maxv
                Cr = 2
            end

            if Jiterrpt[n-1] > MaxIntiter
                Cr += 2
            end

            if JCondrpt[n-1] > MaxJcond
                Cr += 2
            end

            if R1Normrpt[n-1] > R1nmax
                Cr += 2
            end

            if BnormRat[n-1] > MaxBnormRat
                Cr += 2
            end
        end  # end implicit reparameterization criteria

        # Start explicit reparameterization criteria
        if integ > 2
            if vnormrpt[n-1] > Maxv
                Cr = 2
            end

            if ECondrpt[n-1] > MaxECond
                Cr += 2
            end

            if BnormRat[n-1] > MaxBnormRat
                Cr += 2
            end
        end  # end explicit reparameterization criteria

        if Cr > 1
            npar = n - 1
            Crrpt = vcat(Crrpt, Cr)
            nch = 1

            # Parameterization
            tnm = t[n-1]
            v, vd, vdd, q0, U, V, B, jRepar = Param(n, tnm, Q, Qd, Qdd, SJDT, par, jRepar)
            println("Cr>1 V=",V)
            u = zeros(nu)
            Uu[:, n-1] = u
            Vv[:, n-1] = v
            Vvd[:, n-1] = vd
            Vvdd[:, n-1] = vdd

            B0norm = norm(B)
            BCond0 = B0norm * norm(U' * U)

            jReparrpt = vcat(jReparrpt, jRepar)

            Cr = 1
        end
        # Integration methods
        if integ == 1  # Implicit Ind0 Trap
            v, vd, vdd, Lam, R1n, Jiter, JCond, h, nch = ImplicitInd0Trap(n, tn, npar, Vv, Vvd, Vvdd, LLam, Uu, q0, V, U, B, h, hmax, nch, SMDT, SJDT, STSDAT, par)

            R1Normrpt[n] = R1n
            JCondrpt[n] = JCond
            Jiterrpt[n] = Jiter
            hrpt = vcat(hrpt, h)
            LLam[:, n] = Lam
        end

        if integ == 2  # Implicit Ind0 SDIRK54
            v, vd, vdd, Lam, Jiter, R1n, J, JCond, h, nch = ImplicitInd0SDIRK54(n, tn, Vv, Vvd, Vvdd, LLam, Uu, q0, U, V, B, par, h, hmax, nch, npar, SJDT, SMDT, STSDAT)

            R1Normrpt[n] = R1n
            JCondrpt[n] = JCond
            Jiterrpt[n] = Jiter
            hrpt = vcat(hrpt, h)
            LLam[:, n] = Lam
        end
        #println(V)
        if integ == 3  # Explicit Ind0 Nystrom4
            v, vd, vdd, Lam, ECond, jodeiter1 = ExplicitInd0Nystrom4(n, tn, Vv, Vvd, Vvdd, LLam, Uu, V, U, B, q0, h, npar, SJDT, SMDT, STSDAT, par)

            push!(ECondrpt, ECond)
            LLam[:, n] = Lam
            jodeiterRpt = vcat(jodeiterRpt, jodeiter1)
        end

        if integ == 4  # Explicit Ind0 RKFN45
            v, vd, vdd, Lam, ECond, h, nch = ExplicitInd0RKFN45(n, tn, Vv, Vvd, Vvdd, LLam, Uu, U, V, B, q0, h, hmax, par, SMDT, STSDAT, SJDT, nch)

            push!(ECondrpt, ECond)
            LLam[:, n] = Lam
            hrpt = vcat(hrpt, h)
        end

        # Process Results
        Vv[:, n] = v
        Vvd[:, n] = vd
        Vvdd[:, n] = vdd
        vnormrpt = vcat(vnormrpt, norm(v))
        vdnormrpt = vcat(vdnormrpt, norm(vd))
        vddnormrpt = vcat(vddnormrpt, norm(vdd))

        # Evaluate q
        u = Uu[:, n-1]
        u, Iteru = usolv(tn, u, v, q0, SJDT, V, U, B, par)
        Iterurpt = vcat(Iterurpt, Iteru)
        Uu[:, n] = u
        unorm = vcat(unorm, norm(u))
        q = q0 + V * v - U * u
        Q[:, n] = q
        qnormrpt = vcat(qnormrpt, norm(q))

        # Update B and Evaluate qd
        B, Biter = CorrectB(tn, q, B, U, SJDT, par)
        Biterrpt = vcat(Biterrpt, Biter)
        Bnorm = vcat(Bnorm, norm(B))
        BnormRat = vcat(BnormRat, Bnorm[n] / B0norm)
        Phiq = PhiqEval(tn, q, SJDT, par)
        BCond = vcat(BCond, Bnorm[n] * norm(Phiq * U))
        BCondRat = vcat(BCondRat, BCond[n] / BCond0)
        Binvnorm = vcat(Binvnorm, norm(Phiq * U))
        Bcond = vcat(Bcond, Bnorm[n] * Binvnorm[n])
        D = (I - U * B * Phiq) * V
        Pst, Pstt, Pstq, Psttq = P5Eval(tn, q, par)
        qd = D * vd - U * B * Pst
        Qd[:, n] = qd
        qdnormrpt = vcat(qdnormrpt, norm(qd))

        # Evaluate qdd
        Gam = GamEval(tn, q, qd, SJDT, par)
        qdd = D * vdd - U * B * Gam
        Qdd[:, n] = qdd
        qddnormrpt = vcat(qddnormrpt, norm(qdd))
        LLam[:, n] = Lam

        # End Tangent Space Integration

        # Calculate Total Energy
        M = MEval(q, SMDT, par)
        KE = vcat(KE, 0.5 * qd' * M * qd)

        i = 1
        PE = vcat(PE, 0)
        PEn=0.0
        while i <= nb
            PEn += SMDT[1, i] * g * q[7*(i-1)+3]
            i += 1
        end
        PE = vcat(PE, PEn)
        SE = vcat(SE, 0)
        TS = 1
        SEn=0
        while TS <= NTSDA
            i, j, s1pr, s2pr, K, C, el0, F = STSDATPart(STSDAT, TS)
            r1, p1 = qPart(q, i)
            r2 = [0; 0; 0]
            p2 = [1; 0; 0; 0]
            if j >= 1
                r2, p2 = qPart(q, j)
            end
            A1 = ATran(p1)
            A2 = ATran(p2)
            d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
            el = sqrt(d12' * d12)
            SET = 0.5 * K * (el - el0)^2
            SEn += SET
            TS += 1
        end
        SE = vcat(SE, SEn)
        TE = vcat(TE, KE[n] + PE[n] + SE[n])
        # println(n)
        println(tn)
        if app == 4       # Rotating Disk with Translating Body
            theta1[1] = 0
            p1 = [q[4], q[5], q[6], q[7]]
            p1d = [qd[4], qd[5], qd[6], qd[7]]
            r2 = [q[8], q[9], q[10]]
            E1 = EEval(p1)
            A1 = ATran(p1)
            push!(omegaz1, 2 * uz' * E1 * p1d)
            push!(theta1, theta1[n - 1] + h * omegaz1[n])
            push!(dely2pr, (A1' * (r2 - A1 * ux))' * uy)
        end


    end
    if app==4
        df = DataFrame()
        println(size(t))
        println(size(theta1))
        println(size(omegaz1))
        df[!, "t"] = t[1:length(theta1)]
        df[!, "theta1"] = theta1
        df[!, "omegaz1"] = omegaz1
        CSV.write("jl_solver_app4.csv", df)
    end

    if app==6
        df = DataFrame()
        df[!, "t"] = t
        df[!, "x"] = Q[1, :]
        CSV.write("jl_solver.csv", df)
    end



    # println("u₀=", u₀)
    # println("du₀=", du₀)

    # prob2 = ODEProblem(odequation, u₀, tspan, p)

    # # Integration
    # default_solve_kwargs = Dict(:alg => Tsit5(), :reltol => 1e-6, :abstol => 1e-6, :progress => true)
    # sol = solve(prob2; default_solve_kwargs...)
    # return sol
end

function draw(sol)
    f01 = plot(sol, vars=1:3)
    display(f01)
    f02 = plot(sol, vars=8:11)
    display(f02)
    f03 = plot(sol, vars=13:15)
    display(f03)
end

function save(sol)
    col_names = ["x1", "y1", "z1", "p1_1", "p1_2", "p1_3", "p1_4"]

    # Create a DataFrame using the column names
    df = DataFrame()
    df[!, "t"] = sol.t
    df[!, col_names[1]] = sol[1, :]
    df[!, col_names[2]] = sol[2, :]
    df[!, col_names[3]] = sol[3, :]

    CSV.write("jl_solver.csv", df)
end

using ProfileView
#@ProfileView.profview test_EI0(app=4, tspan=(0.0, 0.1))  # run once to trigger compilation (ignore this one)
#@ProfileView.profview test_EI0(app=4, tspan=(0.0, 1.0))
#sol = test_EI0(app=6, tspan=(0.0, 1.0))
test_EI0(app=4, tspan=(0.0, 1.0))
#draw(sol)
#save(sol)