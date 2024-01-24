export parPart,RevPart,TranPart,CylPart,DistPart,SphPartparPart
println("AbstractArray")
function parPart(par)
    nb = par[1]
    ngc = par[2]
    nh = par[3]
    nc = par[4]
    nv=par[5];
    nu=par[6];
    g=par[7];
    utol=par[8];
    Btol=par[9];
    intol=par[10];
    Atol=par[11];
    Vtol=par[12];
    hvar=par[13];
    NTSDA=par[14];
    vt=par[15];
    
    return nb,ngc,nh,nc,nv,nu,g,utol,Btol,intol,Atol,Vtol,hvar,NTSDA,vt
end

function RevPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    vx1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    vz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    vx2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]
    vz2pr = [SJDT[20, k], SJDT[21, k], SJDT[22, k]]
    a = SJDT[23, k]
    R = SJDT[24, k]
    mus = SJDT[25, k]
    mud = SJDT[26, k]
    ms = SJDT[27, k]
    nm = SJDT[28, k]

    return i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr, a, R, mus, mud, ms, nm
end


function TranPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = Vector{Float64}([SJDT[4, k], SJDT[5, k], SJDT[6, k]])
    s2pr = Vector{Float64}([SJDT[7, k], SJDT[8, k], SJDT[9, k]])
    vx1pr = Vector{Float64}([SJDT[11, k], SJDT[12, k], SJDT[13, k]])
    vz1pr = Vector{Float64}([SJDT[14, k], SJDT[15, k], SJDT[16, k]])
    vx2pr = Vector{Float64}([SJDT[17, k], SJDT[18, k], SJDT[19, k]])
    vz2pr = Vector{Float64}([SJDT[20, k], SJDT[21, k], SJDT[22, k]])
    a = SJDT[23, k]
    b = SJDT[24, k]
    mus = SJDT[25, k]
    mud = SJDT[26, k]
    ms = SJDT[27, k]
    nm = SJDT[28, k]

    return i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr, a, b, mus, mud, ms, nm
end
function CylPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = Vector{Float64}([SJDT[4, k], SJDT[5, k], SJDT[6, k]])
    s2pr = Vector{Float64}([SJDT[7, k], SJDT[8, k], SJDT[9, k]])
    vx1pr = Vector{Float64}([SJDT[11, k], SJDT[12, k], SJDT[13, k]])
    vz1pr = Vector{Float64}([SJDT[14, k], SJDT[15, k], SJDT[16, k]])
    vx2pr = Vector{Float64}([SJDT[17, k], SJDT[18, k], SJDT[19, k]])
    vz2pr = Vector{Float64}([SJDT[20, k], SJDT[21, k], SJDT[22, k]])
    a = SJDT[23, k]
    R = SJDT[24, k]
    mus = SJDT[25, k]
    mud = SJDT[26, k]
    ms = SJDT[27, k]
    nm = SJDT[28, k]

    return i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr, a, R, mus, mud, ms, nm
end

function DistPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = Vector{Float64}([SJDT[4, k], SJDT[5, k], SJDT[6, k]])
    s2pr = Vector{Float64}([SJDT[7, k], SJDT[8, k], SJDT[9, k]])
    d = SJDT[10, k]
    ms = SJDT[25, k]
    nm = SJDT[26, k]

    return i, j, s1pr, s2pr, d, ms, nm
end

function SphPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]

    return i, j, s1pr, s2pr
end

function qPart(q, i)
    r = q[7*(i-1)+1:7*(i-1)+3]
    p = q[7*(i-1)+4:7*(i-1)+7]
    return r, p
end

function STSDATPart(STSDAT, T)
    i = STSDAT[1, T]
    j = STSDAT[2, T]
    sipr = [STSDAT[3, T], STSDAT[4, T], STSDAT[5, T]]
    sjpr = [STSDAT[6, T], STSDAT[7, T], STSDAT[8, T]]
    K = STSDAT[9, T]
    C = STSDAT[10, T]
    el0 = STSDAT[11, T]
    F = STSDAT[12, T]

    return i, j, sipr, sjpr, K, C, el0, F
end
