using LinearAlgebra
using Plots
using DifferentialEquations
using OrdinaryDiffEq, ProgressLogging
using Sundials
using Printf
using LSODA
export odequation



# tspan = (0.0, 1.0)
# "₁₂₃₄₅₆₇₈₉₀"

function ode_rhs(tnm,qnm,qdnm,SMDT, STSDAT, SJDT, par)
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
    #println("A\\b=",A\b)
    return A\b
end


function odequation(du, u, p, t)
    SMDT, STSDAT, SJDT, par = p
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = mathfunction.parPart(par)
    
    sec1=1:nb*7
    sec2=nb*7+1:nb*7+nc
    sec3=nb*7+nc+1:nb*14+nc

    dq=u[sec3]
    q=u[sec1]
    l=u[sec2]

    res=ode_rhs(t,q,dq,SMDT, STSDAT, SJDT, par)
    du[sec1]=u[sec3]
    du[sec2]=res[sec2]
    du[sec3] = res[sec1]

end


function sol_O(out, p, t)

    tspan = (0.0, 1.0)
    du₀, u₀=p.du₀,p.u₀ 
    prob = ODEProblem(odequation,  u₀, tspan, p)
    sol = solve(prob, Tsit5(), reltol=1e-15, abstol=1e-15, progress=true)

    return sol
end