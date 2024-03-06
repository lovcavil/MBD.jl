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
    # A=[m1 0 0 q[1];
    #    0 m1 0 q[2] ;
    #    0 0 m1 q[3];
    #    q[1] q[2] q[3] 0]
    # Qa = [0, 0, -m1 * params.g]
    # b=vcat(Qa, -q_v[1]^2 - q_v[2]^2 - q_v[3]^2)


function odequation(du, u, p, t)
    SMDT, STSDAT, SJDT, par = p
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = mathfunction.parPart(par)

    sec1=1:nb*7
    sec2=nb*7+1:nb*7+nc
    sec3=nb*7+nc+1:nb*14+nc

    dq=u[sec3]
    q=u[sec1]
    l=u[sec2]

    Phiq = mathfunction.PhiqEval(t, q, SJDT, par)
    QA = mathfunction.QAEval(t, q, dq, SMDT, STSDAT, par)
    M = mathfunction.MEval(q, SMDT, par)
    Gam = mathfunction.GamEval(t, q, dq, SJDT, par)
    A=vcat( hcat(M,         Phiq'           ),
            hcat(Phiq,      zeros(nc,nc)    )   )
    b=vcat(QA, -Gam)
    res= A \ b
    du[sec1]=u[sec3]
    du[sec2]=res[sec2]
    du[sec3] = res[sec1]

end
