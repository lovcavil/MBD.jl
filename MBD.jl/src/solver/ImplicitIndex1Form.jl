using LinearAlgebra
using Plots
using DifferentialEquations
using OrdinaryDiffEq, ProgressLogging
using Sundials
using Printf
using LSODA
export implicitdae

# tspan = (0.0, 1.0)
# "₁₂₃₄₅₆₇₈₉₀"ᵃ ᵇ ᶜ ᵈ ᵉ ᶠ ᵍ ʰ ⁱ ʲ ᵏ ˡ ᵐ ⁿ ᵒ ᵖ ʳ ˢ ᵗ ᵘ ᵛ ʷ ˣ ʸ ᶻ ⁰ ¹ ² ³ ⁴ ⁵ ⁶ ⁷ ⁸ ⁹ ⁺ ⁻ ⁼ ⁽ ⁾

function rhs_II1(x,tnm,qnm,qdnm,SMDT, STSDAT, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = mathfunction.parPart(par)
    Phiq = mathfunction.PhiqEval(tnm, qnm, SJDT, par)
    #Gamsq, Gamsqd = mathfunction.GamsqqdEval(tnm, qnm, qdnm, SJDT, par)
    #P2 = mathfunction.P2Eval(tnm, qnm, qddnm, SJDT, par)
    #M2 = mathfunction.M2Eval(qnm, qddnm, SMDT, par)
    #QAsq, QAsqd = mathfunction.QAsqqd(tnm, qnm, qdnm, SMDT, STSDAT, par)
    #Ssq, Ssqd = mathfunction.SqqdEval(qnm, qdnm, SMDT, par)
    #P4 = mathfunction.P4Eval(tnm, qnm, Lamnm, SJDT, par)
    QA = mathfunction.QAEval(tnm, qnm, qdnm, SMDT, STSDAT, par)
    #S = mathfunction.SEval(qnm, qdnm, SMDT, par)
    M = mathfunction.MEval(qnm, SMDT, par)
    Gam = mathfunction.GamEval(tnm, qnm, qdnm, SJDT, par)

    # A=[m1 0 0 q[1];
    #    0 m1 0 q[2] ;
    #    0 0 m1 q[3];
    #    q[1] q[2] q[3] 0]
    # Qa = [0, 0, -m1 * params.g]
    # b=vcat(Qa, -q_v[1]^2 - q_v[2]^2 - q_v[3]^2)

    A=vcat(hcat(M,Phiq'),
    hcat(Phiq,zeros(nc,nc)))
    b=vcat(QA, -Gam)
    #println("A*x-b=",A*x-b)
    return A*x-b
end


function implicitdae(out, du, u, p, t)
    @unpack SMDT, STSDAT, SJDT, par = p
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = mathfunction.parPart(par)
    
    sec1=1:nb*7
    sec2=nb*7+1:nb*7+nc
    sec3=nb*7+nc+1:nb*14+nc
    ddq=du[sec3]
    dq=du[sec1]
    v=u[sec3]
    q=u[sec1]
    l=u[sec2]

    x=vcat(ddq,l)
    res=rhs_II1(x,t,q,dq,SMDT, STSDAT, SJDT, par)
    out[sec1]=res[sec1]
    out[sec2]=res[sec2]
    out[sec3] = v - dq
    #@printf("residual=")
    #for element in out
    #    @printf("%8.2e ", element)
    #end
    #println("")
    #println("t=",t)
    println(out)
end


function sol_I(out, p, t)
    #equation(out, du₀, u₀, p, t)
    #println(out)
    tspan = (0.0, 1.0)
    du₀, u₀=p.du₀,p.u₀ 
    prob = DAEProblem(p.equation, du₀, u₀, tspan, p, differential_vars=p.differential_vars)

    #sol = solve(prob, DFBDF() , reltol=1e-6, abstol=1e-6, progress=true)
    sol = solve(prob, IDA(max_num_steps_ic=10), progress=true)
    #sol = solve(prob, lsoda(), progress=true)
    # # Create a plot(init_all=false)
    # p = plot(title="Solution to the linear ODE with a thick line", xaxis="Time (t)", yaxis="u(t) (in μm)")

    # # Plot only the first 5 curves
    # for i in 1:min(3, length(sol))
    #     plot!(p, sol.t, sol[i, :], linewidth=3, label="Curve $i")
    # end

    # # Display the plot
    # display(p)
    return sol
end

