using LinearAlgebra
using Plots
using DifferentialEquations
using OrdinaryDiffEq, ProgressLogging
using Sundials
using Printf
using LSODA
export implicit3dae

# tspan = (0.0, 1.0)
# "₁₂₃₄₅₆₇₈₉₀"

function rhs(x,tnm,qnm,qdnm,SMDT, STSDAT, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = mathfunction.parPart(par)
    Phiq = mathfunction.PhiqEval(tnm, qnm, SJDT, par)
    QA = mathfunction.QAEval(tnm, qnm, qdnm, SMDT, STSDAT, par)
    M = mathfunction.MEval(qnm, SMDT, par)
    A=hcat(M,Phiq')
    #println("A*x-b=",A*x-b)
    return A*x-QA
end


function implicit3dae(out, du, u, p, t)
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
    Phi=mathfunction.PhiEval(t, q, SJDT, par)
    x=vcat(ddq,l)
    res=rhs(x,t,q,dq,SMDT, STSDAT, SJDT, par)
    out[sec1]=res[sec1]
    out[sec2]=Phi
    out[sec3] = v - dq
    println(out)
end




