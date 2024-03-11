using LinearAlgebra
using Plots
using DifferentialEquations
using OrdinaryDiffEq, ProgressLogging
#using Sundials
using Printf
#using LSODA
export odequation

function odequation(du, u, p, t)
    # println("t $t")
    SMDT, STSDAT, SJDT, par, p_contact = p
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = mathfunction.parPart(par)

    sec1=1:nb*7
    sec2=nb*7+1:nb*7+nc
    sec3=nb*7+nc+1:nb*14+nc

    dq=u[sec3]
    q=u[sec1]
    l=u[sec2]

    Phiq = mathfunction.PhiqEval(t, q, SJDT, par)
    QA = mathfunction.QAEval(t, q, dq, SMDT, STSDAT, par)
    QAC = mathfunction.QACEval(t, q, dq, SMDT, STSDAT, par,p_contact)
    M = mathfunction.MEval(q, SMDT, par)
    Gam = mathfunction.GamEval(t, q, dq, SJDT, par)
    A=vcat( hcat(M,         Phiq'           ),
            hcat(Phiq,      zeros(nc,nc)    )   )
    b=vcat(QA + QAC, -Gam)
    res= A \ b
    du[sec1]=u[sec3]
    du[sec2]=res[sec2]
    du[sec3] = res[sec1]
end
