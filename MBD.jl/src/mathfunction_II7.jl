# mathfunction.jl
module mathfunction
using Test
using LinearAlgebra
using BlockDiagonals
using CSV, DataFrames
include("./constraint/constraint_fPy.jl")
include("./constraint/constraint_fPz.jl")
include("./constraint/constraint_fP.jl")
include("./constraint/constraint_Poly.jl")
include("./eval/QAC.jl")
export add_constraint!
export atilde
export ATran,BTran,qPart,bbP2dist,bbP2dot1,bbP2dot2,bbP2RotDr,bbP2sph
export bbPhidist,bbPhiqdist,bbPhidot1,bbPhidot2,bbPhiqdot1,bbPhiqdot2
export bbPhiqRotDr,bbPhiqsph,InitConfig,PhiqEval,PhiEval#,einv
# function einv(A,b)
#     x=[]
#     # try
#     #     x = inv(A) * b  # 尝试执行此操作
#     # catch e
#     #     if e isa SingularException
#     #         println("发生1 SingularException：", e)
#     #         x = inv(A+0.01*I) * b  
#     #         # 处理 SingularException 的代码
#     #     elseif e isa MethodError
#     #         x = pinv(A) * b  # 如果上面的代码失败，执行这行代码
#     #         println("发生2 MethodError：", e)
#     #         # 处理 MethodError 的代码
#     #     else
#     #         x = pinv(A) * b  # 如果上面的代码失败，执行这行代码
#     #         println("发生异常3，错误信息：", e)  # 打印异常信息
#     #     end

#     # end

#     println("rank",rank(A))
#     try
#         println("eI")
#         x = (A+0.001*I) \ b
        
#     catch e
#         println(e)
#         println("w = -(JJ+0.0001*I) \\ Resid")
#         x = (A+0.01*I) \ b
#     end
#     try
#         println("err_\\",norm(A \ b))
#     catch e
#     end
#     try
#         println("err_\\I",norm((A+0.001*I) \ b))
#     catch e
#     end
#     try
#         println("err_inv",norm(inv(A) * b))
#     catch e
#     end
#     try
#         println("err_pinv",norm(pinv(A) * b))
#     catch e
#     end
#     return x
# end
# Define the function from the previous translation
function add_constraint!(A, B, m, n)
    r, s = size(B, 1), size(B, 2)
    for i in 1:r
        for j in 1:s
            A[m+i, n+j] += B[i, j]
        end
    end
    return A
end


function atilde(a)
    return [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
end

function ATran(p)
    e0 = p[1]
    e = p[2:4]
    I3 = I  # The identity matrix
    etil = atilde(e)
    AT = (e0^2 - dot(e, e)) * I3 + 2 * e * e' + 2 * e0 * etil
    return AT
end

#= p = [1, 2, 3,0]
at = ATran(p)
println(at) =#

function BTran(p, apr)
    e0 = p[1]
    e = p[2:4]
    I3 = I  # Identity matrix in Julia is represented by I
    etil = atilde(e)
    BT = 2 * hcat((e0 * I3 + etil) * apr, e * apr' - (e0 * I3 + etil) * atilde(apr))
    return BT
end

#= # Example usage:
p = [1, 2, 3, 4]  # Example quaternion components
apr = [5, 6, 7]   # Example vector
BT_matrix = BTran(p, apr)
println("The BTran matrix is:")
println(BT_matrix) =#

function qPart(q, i)
    r = q[7*(i-1)+1:7*(i-1)+3]
    p = q[7*(i-1)+4:7*(i-1)+7]
    return r, p
end

# Example usage:
#q = rand(14) # Example vector of length 14, which can hold data for two 'i' indices
#i = 2        # We want to extract the second set of (r, p)

#= r, p = qPart(q, i)
println("Vector r: ", r)
println("Quaternion p: ", p) =#

function parPart(par)
    nb = par[1]
    ngc = par[2]
    nh = par[3]
    nc = par[4]
    g = par[5]
    intol = par[6]
    Atol = par[7]
    h0=par[8]
    hvar=par[9]
    NTSDA=par[10]

    return nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA
end


function bbP2dist(i, j, s1pr, s2pr, d,tn, q, qd, par)
    # Assuming par is a structure or a dictionary with the following keys.
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1x = BTran(xp1, s1pr)
    a1 = (xr1' + xp1' * BT1')

    P21 = zeros(1, 7)
    P22 = zeros(1, 7)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P21 = [a1, a1 * BT1 - d12' * BT1x]
    elseif j >= 1
        r2, p2 = qPart(q, j)
        xr2, xp2 = qPart(qd, j)
        A2 = ATran(p2)
        BT2 = BTran(p2, s2pr)
        BT2x = BTran(xp2, s2pr)
        a2 = (xr2' + xp2' * BT2')
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P21 = hcat((a1 - a2), (a1 - a2) * BT1 - d12' * BT1x)
        P22 = hcat((a2 - a1), (a2 - a1) * BT2 + d12' * BT2x)
    end

    return P21, P22
end

function bbP2dot1(i, j, a1pr, a2pr,tn, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    
    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, a1pr)
    BTx1 = BTran(xp1, a1pr)
    
    P21 = zeros(1, 7)
    P22 = zeros(1, 7)
    
    if j == 0
        P21 = hcat(zeros(1, 3), a2pr' * BTx1)
    elseif j >= 1
        r2, p2 = qPart(q, j)
        xr2, xp2 = qPart(qd, j)
        A2 = ATran(p2)
        BT2 = BTran(p2, a2pr)
        BTx2 = BTran(xp2, a2pr)
        P21 = hcat(zeros(1, 3), a2pr' * A2' * BTx1 + xp2' * BT2' * BT1)
        P22 = hcat(zeros(1, 3), a1pr' * A1' * BTx2 + xp1' * BT1' * BT2)
    end
    
    return P21, P22
end


function bbP2dot2(i, j, a2pr, s1pr, s2pr,tn, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    BT1s1 = BTran(p1, s1pr)
    xr1, xp1 = qPart(qd, i)
    BTx1s1 = BTran(xp1, s1pr)
    
    P21 = zeros(1, 7)
    P22 = zeros(1, 7)
    
    if j == 0
        P21 = -hcat(zeros(1, 3), a2pr' * BTx1s1)
    elseif j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        xr2, xp2 = qPart(qd, j)
        BT2s2 = BTran(p2, s2pr)
        BT2a2 = BTran(p2, a2pr)
        BTx2s2 = BTran(xp2, s2pr)
        BTx2a2 = BTran(xp2, a2pr)
        c = xp2' * BT2a2' * BT2s2 + d12' * BTx2a2 + a2pr' * A2' * BTx2s2 + (xr2' + xp2' * BT2s2') * BT2a2
        P21 = -hcat(xp2' * BT2a2', a2pr' * A2' * BTx1s1 + xp2' * BT2a2' * BT1s1)
        P22 = hcat(xp2' * BT2a2', c - (xr1' + xp1' * BT1s1') * BT2a2)
    end
    
    return P21, P22
end

function bbP2fxc(i, j, s1pr, s2pr, d,tn, q, qd, par)
    # Assuming par is a structure or a dictionary with the following keys.
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    E11=[1 0 0;0 0 0;0 0 0]
    r1, p1 = qPart(q, i)
    xr1, xp1 = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1x = BTran(xp1, s1pr)
    a1 = (xr1' + xp1' * BT1')
    a1E11 = a1 * E11
    P21 = zeros(1, 7)
    P22 = zeros(1, 7)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P21 = hcat(a1, a1 * BT1 - d12' * BT1x)
    elseif j >= 1
        r2, p2 = qPart(q, j)
        xr2, xp2 = qPart(qd, j)
        A2 = ATran(p2)
        BT2 = BTran(p2, s2pr)
        BT2x = BTran(xp2, s2pr)
        a2 = (xr2' + xp2' * BT2')
        a2E11 = a2 * E11
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P21 = hcat((a1E11 - a2E11), (a1E11 - a2E11) * BT1 - d12' *E11* BT1x)
        P22 = hcat((a2E11 - a1E11), (a2E11 - a1E11) * BT2 + d12' *E11* BT2x)
    end
    # Special
    #println("xr1=",xr1')
    P21 = hcat(xr1[1],zeros(1, 6))
    return P21, P22
end



function bbP2RotDr(i, j, vx1pr, vy1pr, vx2pr, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    s = bbPhidot1(i, j, vy1pr, vx2pr,tn, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr,tn, q, par)

    if abs(c) >= abs(s)
        P21, P22 = bbP2dot1(i, j, vy1pr, vx2pr,tn, q, qd, par)
    else
        P21, P22 = bbP2dot1(i, j, vx1pr, vx2pr,tn, q, qd, par)
    end

    return P21, P22
end

# Example usage:
# You would need to define the variables `i`, `j`, `vx1pr`, `vy1pr`, `vx2pr`, `q`, `qd`, and `par`
# appropriately before calling this function.
# P21, P22 = bbP2RotDr(i, j, vx1pr, vy1pr, vx2pr, q, qd, par)
function bbP2sph(i, j, s1pr, s2pr,tn, q, qd, par)
    # Unpacking parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    xr1, xp1 = qPart(qd, i)
    BT1x = BTran(xp1, s1pr)

    P21 = zeros(3, 7)
    P22 = zeros(3, 7)
    
    if j == 0
        P21[:, 4:7] = -BT1x
    elseif j >= 1
        xr2, xp2 = qPart(qd, j)
        BT2x = BTran(xp2, s2pr)
        P21[:, 4:7] = -BT1x
        P22[:, 4:7] = BT2x
    end

    return P21, P22
end

function bbP3dist(i, j, s1pr, s2pr, d, tn, q, qd, par)
    I3 = Matrix{Float64}(I, 3, 3)  # Equivalent of eye(3) in MATLAB

    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    BT1 = BTran(p1, s1pr)
    BT1d = BTran(p1d, s1pr)
    a1bar = r1d + BT1 * p1d
    if j == 0
        P31 = hcat(p1d' * BT1d', 2 * (a1bar)' * BT1d + p1d' * BT1d' * BT1)
        P32 = zeros(1, 7)
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        r2d, p2d = qPart(qd, j)
        BT2 = BTran(p2, s2pr)
        BT2d = BTran(p2d, s2pr)
        a2bar = r2d + BT2 * p2d
        P31 = hcat(p1d' * BT1d' - p2d' * BT2d', 2 * (a1bar - a2bar)' * BT1d +
               (p1d' * BT1d' - p2d' * BT2d') * BT1)
        P32 = hcat(p2d' * BT2d' - p1d' * BT1d', 2 * (a2bar - a1bar)' * BT2d +
               (p2d' * BT2d' - p1d' * BT1d') * BT2)
    end

    return P31, P32
end


function bbP3dot1(i, j, a1pr, a2pr, tn, q, qd, par)
    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    BT1 = BTran(p1, a1pr)
    BT1d = BTran(p1d, a1pr)

    P31 = zeros(1, 7)
    P32 = zeros(1, 7)

    if j >= 1
        r2, p2 = qPart(q, j)
        r2d, p2d = qPart(qd, j)
        BT2 = BTran(p2, a2pr)
        BT2d = BTran(p2d, a2pr)
        P31 = hcat(0, 0, 0, 2 * (p2d'* BT2' * BT1d) + (p2d'* BT2d' * BT1))# from source code error
        P32 = hcat(0, 0, 0, 2 * (p1d'* BT1' * BT2d) + (p1d'* BT1d' * BT2))
    end 

    return P31, P32
end
function bbP3dot2(i, j, a2pr, s1pr, s2pr, tn, q, qd, par)
    # nb, ngc, nh, nc, nv, nu, g, utol, Btol, intol, Atol, Vtol, hvar, NTSDA, vt = parPart(par)

    I3 = Matrix{Float64}(I, 3, 3)
    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)
    BT1d = BTran(p1d, s1pr)

    P31 = zeros(1, 7)
    P32 = zeros(1, 7)

    if j >= 1
        r2, p2 = qPart(q, j)
        r2d, p2d = qPart(qd, j)
        A2 = ATran(p2)
        BT2s2 = BTran(p2, s2pr)
        BT2a2 = BTran(p2, a2pr)
        BT2s2d = BTran(p2d, s2pr)
        BT2a2d = BTran(p2d, a2pr)
        d = p2d' * BT2a2d' * BT2s2 + p2d' * BT2s2d' * BT2a2 + 2 * r2d' * BT2a2d +
            2 * p2d' * BT2s2' * BT2a2d + 2 * p2d' * BT2a2' * BT2s2d
        P31 = -hcat(p2d' * BT2a2', 2 * p2d' * BT2a2' * BT1d + p2d' * BT2a2d' * BT1)
        P32 = hcat(p2d' * BT2a2', -p1d' * BT1d' * BT2a2 - 2 * (r1d' + p1d' * BT1') * BT2a2d + d)
    end

    return P31, P32
end

function bbP3fxc(i, j, s1pr, s2pr, d, tn, q, qd, par)
    I3 = Matrix{Float64}(I, 3, 3)  # Equivalent of eye(3) in MATLAB

    r1, p1 = qPart(q, i)
    r1d, p1d = qPart(qd, i)
    BT1 = BTran(p1, s1pr)
    BT1d = BTran(p1d, s1pr)
    a1bar = r1d + BT1 * p1d
    E11=[1 0 0;0 0 0;0 0 0]
    a1barE11=a1bar'*E11
    if j == 0
        P31 = hcat(p1d' * BT1d'*E11, 2 * (a1barE11) * BT1d + p1d'* BT1d' *E11 * BT1)
        P32 = zeros(1, 7)
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        r2d, p2d = qPart(qd, j)
        BT2 = BTran(p2, s2pr)
        BT2d = BTran(p2d, s2pr)
        a2bar = r2d + BT2 * p2d
        a2barE11=a2bar*E11
        P31 = hcat(p1d' * BT1d'*E11 - p2d' * BT2d'*E11, 2 * (a1barE11 - a2barE11)' * BT1d +
               (p1d' * BT1d' - p2d' * BT2d') *E11* BT1)
        P32 = hcat(p2d' * BT2d'*E11 - p1d' * BT1d'*E11, 2 * (a2barE11 - a1barE11)' * BT2d +
               (p2d'* BT2d' - p1d' * BT1d') *E11 * BT2)
    end
    # speical
    P31 = zeros(1, 7)
    println("nnnnoooo")
    return P31, P32
end

function bbP4dist(i, j, s1pr, s2pr, d, tn, q, etak, par)
    I3 = Matrix{Float64}(I, 3, 3)

    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3,BT1), hcat(BT1',BT1' * BT1 - KEval(s1pr, d12)))
        P412 = zeros(7, 7)
        P422 = zeros(7, 7)
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        BT2 = BTran(p2, s2pr)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3, BT1), hcat(BT1', BT1' * BT1 - KEval(s1pr, d12)))
        P412 = -etak * vcat(hcat(I3, BT2), hcat(BT1', BT1' * BT2))
        P422 = etak * vcat(hcat(I3, BT2), hcat(BT2', BT2' * BT2 + KEval(s2pr, d12))    )
    end

    return P411, P412, P422
end

function bbP4dot1(i, j, a1pr, a2pr, tn, q, etak, par)
    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, a1pr)
    A1 = ATran(p1)

    P411 = zeros(7, 7)
    P412 = zeros(7, 7)
    P422 = zeros(7, 7)

    if j == 0
        P411 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3), KEval(a1pr, a2pr)))
    elseif j >= 1
        r2, p2 = qPart(q, j)
        BT2 = BTran(p2, a2pr)
        A2 = ATran(p2)

        P411 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3), KEval(a1pr, A2 * a2pr)))
        P412 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3), BT1' * BT2))
        P422 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3), KEval(a2pr, A1 * a1pr)))
    end

    return P411, P412, P422
end
function bbP4dot2(i, j, a2pr, s1pr, s2pr, tn, q, etak, par)
    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)

    P411 = zeros(7, 7)
    P412 = zeros(7, 7)
    P422 = zeros(7, 7)

    if j == 0 #to ch8
        P411 = -etak * vcat(zeros(3, 7), hcat(zeros(4, 3), KEval(s1pr, a2pr)))
    elseif j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        BT2a2 = BTran(p2, a2pr)
        BT2s2 = BTran(p2, s2pr)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        e = KEval(s2pr, A2 * a2pr) + KEval(a2pr, d12) + BT2s2' * BT2a2 + BT2a2' * BT2s2
        P411 = -etak * vcat(zeros(3, 7), hcat(zeros(4, 3),KEval(s1pr, A2 * a2pr)))#to ch8
        P412 = -etak * vcat(hcat(zeros(3, 3), BT2a2), hcat(zeros(4, 3), BT1' * BT2a2))
        P422 = etak * vcat(hcat(zeros(3, 3), BT2a2), hcat(BT2a2', e))
    end

    return P411, P412, P422
end
function bbP4fxc(i, j, s1pr, s2pr, d, tn, q, etak, par)
    I3 = Matrix{Float64}(I, 3, 3)

    r1, p1 = qPart(q, i)
    BT1 = BTran(p1, s1pr)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3*E11,E11*BT1), hcat(BT1'*E11,BT1'*E11 * BT1 - KEval(s1pr, E11*d12)))
        P412 = zeros(7, 7)
        P422 = zeros(7, 7)
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        BT2 = BTran(p2, s2pr)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        P411 = etak * vcat(hcat(I3*E11, E11*BT1), hcat(BT1'*E11, BT1' *E11* BT1 - KEval(s1pr, E11*d12)))
        P412 = -etak * vcat(hcat(I3*E11, E11*BT2), hcat(BT1'*E11, BT1' *E11* BT2))
        P422 = etak * vcat(hcat(I3*E11, E11*BT2), hcat(BT2'*E11, BT2' *E11* BT2 + KEval(s2pr, E11*d12))    )
    end
    P411 = etak * vcat(hcat(I3*E11,zeros(3,4)), hcat(zeros(4,3),- KEval(s1pr,  vcat(r1[1],zeros(2,1)) )))
    println("nnnnoooo4")
    return P411, P412, P422
end

function bbP4sph(i, j, s1pr, s2pr, tn, q, etak, par)
    # P412 is identically zero
    P412 = zeros(7, 7)

    if j == 0
        P411 = vcat(zeros(3, 7), hcat(zeros(4, 3),-KEval(s1pr, etak)))
        P422 = zeros(7, 7)
    end

    if j >= 1
        P411 = vcat(zeros(3, 7), hcat(zeros(4, 3),-KEval(s1pr, etak)))
        P422 = vcat(zeros(3, 7), hcat(zeros(4, 3),KEval(s2pr, etak)))
    end

    return P411, P412, P422
end

# Example usage:
# You would need to define the variables `i`, `j`, `s1pr`, `s2pr`, `q`, `qd`, and `par`
# appropriately before calling this function.
# P21, P22 = bbP2sph(i, j, s1pr, s2pr, q, qd, par)
function bbPhidist(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phi = (dot(d12, d12) - d^2) / 2
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phi = (dot(d12, d12) - d^2) / 2
    end
    return Phi
end

function bbPhiqdist(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * hcat(I3, BTran(p1, s1pr))
        Phiq2 = zeros(1, 7)
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * hcat(I3, BTran(p1, s1pr))
        Phiq2 = d12' * hcat(I3, BTran(p2, s2pr))
    end
    return Phiq1, Phiq2
end

function bbPhifxc(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phi = ((d12'*E11* d12) - d^2) / 2
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phi = ((d12'*E11* d12)  - d^2) / 2
    end
    Phi=(r1[1]*r1[1]-1) / 2
    return Phi
end

function bbPhiqfxc(i, j, s1pr, s2pr, d,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    I3 = Matrix{Float64}(I, 3, 3) # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    E11=[1 0 0;0 0 0;0 0 0]
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * E11*hcat(I3, BTran(p1, s1pr))
        Phiq2 = zeros(1, 7)
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phiq1 = -d12' * E11* hcat(I3, BTran(p1, s1pr))
        Phiq2 = d12' * E11* hcat(I3, BTran(p2, s2pr))
    end
    #println("bbPhiqfxc-r1",r1)
    Phiq1 = r1'*hcat(I3, zeros(3,4))
    Phiq1 = hcat(r1[1],zeros(1,2), zeros(1,4))
    #println("bbPhiqfxc",Phiq1)
    return Phiq1, Phiq2
end

function bbPhidot1(i, j, a1pr, a2pr,tn, q, par)
    # Unpack parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Get position and orientation of body i
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)

    Phi = 0.0
    if j == 0
        Phi = a1pr' * A1' * a2pr
    elseif j >= 1
        # Get position and orientation of body j
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        Phi = a1pr' * A1' * A2 * a2pr
    end

    return Phi
end

function bbPhidot2(i, j, a2pr, s1pr, s2pr,tn, q, par)
    # Unpack parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Get position and orientation of body i
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)

    Phi = 0.0
    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phi = a2pr' * d12
    elseif j >= 1
        # Get position and orientation of body j
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phi = a2pr' * A2' * d12
    end

    return Phi
end
function bbPhiqdot1(i, j, a1pr, a2pr,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, a1pr)

    Phiq1 = zeros(1, 7)
    Phiq2 = zeros(1, 7)

    if j == 0
        Phiq1[4:7] = (a2pr' * BT1)'
    elseif j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        BT2 = BTran(p2, a2pr)
        Phiq1[4:7] = (a2pr' * A2' * BT1)'
        Phiq2[4:7] = (a1pr' * A1' * BT2)'
    end

    return Phiq1, Phiq2
end
function bbPhiqdot2(i, j, a2pr, s1pr, s2pr,tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix in Julia
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)
    BT1 = BTran(p1, s1pr)

    Phiq1 = zeros(1, 7)
    Phiq2 = zeros(1, 7)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phiq1[1:7] = (-a2pr' * hcat(I3, BT1))'
    elseif j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        BT2s2 = BTran(p2, s2pr)
        BT2a2 = BTran(p2, a2pr)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phiq1[1:7] = (-a2pr' * A2' * hcat(I3, BT1))'
        Phiq2[1:7] = (a2pr' * A2' * hcat(I3, BT2s2) + d12' * hcat(zeros(3, 3), BT2a2))'
    end

    return Phiq1, Phiq2
end
function bbPhiRotDr(i, j, vx1pr, vy1pr, vx2pr, q, par)
    # Unpack parameters; parPart is assumed to be a function that extracts these values
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Compute s and c using the bbPhidot1 function which needs to be defined in Julia
    s = bbPhidot1(i, j, vy1pr, vx2pr,tn, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr,tn, q, par)

    # Use conditional statement to decide which value to return
    Phi = abs(c) >= abs(s) ? s : c

    return Phi
end

function bbPhiqRotDr(i, j, vx1pr, vy1pr, vx2pr,tn, q, par)
    # Unpack parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Compute scalar values of the constraint derivatives
    s = bbPhidot1(i, j, vy1pr, vx2pr,tn, q, par)
    c = bbPhidot1(i, j, vx1pr, vx2pr,tn, q, par)

    # Decide which derivative to use based on the magnitude comparison
    if abs(c) >= abs(s)
        return bbPhiqdot1(i, j, vy1pr, vx2pr,tn, q, par)
    else
        return bbPhiqdot1(i, j, vx1pr, vx2pr,tn, q, par)
    end
end
function bbPhisph(i, j, s1pr, s2pr, tn, q, par)
    r1, p1 = qPart(q, i)
    A1 = ATran(p1)

    if j == 0
        d12 = s2pr - r1 - A1 * s1pr
        Phi = d12
    end

    if j >= 1
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        Phi = d12
    end

    return Phi
end


function bbPhiqsph(i, j, s1pr, s2pr,tn, q, par)
    # Unpack parameters
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Identity matrix in Julia
    I3 = Matrix{Float64}(I, 3, 3)

    # Constraint Jacobians
    Phiq1 = zeros(3, 7)
    Phiq2 = zeros(3, 7)

    # Get the first part of the generalized coordinate for body i
    r1, p1 = qPart(q, i)

    if j == 0
        Phiq1[:, 1:7] = -hcat(I3, BTran(p1, s1pr))
        # Phiq2 remains zero as initialized
    elseif j >= 1
        r2, p2 = qPart(q, j)
        Phiq1[:, 1:7] = -hcat(I3, BTran(p1, s1pr))
        Phiq2[:, 1:7] = hcat(I3, BTran(p2, s2pr))
    end

    return Phiq1, Phiq2
end

function CylPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    ux1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    uz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    ux2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]
    uz2pr = [SJDT[20, k], SJDT[21, k], SJDT[22, k]]

    return i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr
end
function DistDrPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    return i, j, s1pr, s2pr
end
function DistPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    d = SJDT[10, k]
    return i, j, s1pr, s2pr, d
end

function SphPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]

    return i, j, s1pr, s2pr
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


function EEval(p)
    e0 = p[1]
    e = p[2:4]
    Ebar = [-e atil(e) + e0 * I]  # I in Julia is the identity matrix
    return Ebar
end

# Helper function 'atil' that was used in EEval
function atil(a)
    return [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
end

function GamEval(tn, q, qd, SJDT, par)
    Pst, Pstt, Pstq, Psttq = P5Eval(tn, q, par)
    P2 = P2Eval(tn, q, qd, SJDT, par)
    Gam = P2 * qd + Pstt

    return Gam
end


function GEval(x)
    e0 = x[1]
    e = x[2:4]
    etil = atil(e)
    G = hcat(-e,-etil + e0 * I)
    return G
end
function GamsqqdEval(tn, q, qd, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    P2 = P2Eval(tn, q, qd, SJDT, par)
    P3 = P3Eval(tn, q, qd, SJDT, par)
    Gamsq = P3
    Gamsqd = 2 * P2
    return Gamsq, Gamsqd
end

function InitConfig(rO, rP, rQ)
    # Define A using point definition of Sect 2.5.5.
    f = (1 / norm(rP - rO)) * (rP - rO)
    h = (1 / norm(atil(f) * (rQ - rO))) * atil(f) * (rQ - rO)
    g = -atil(f) * h
    A = [f g h]

    # Compute Euler parameter vector p in Section 2.5.4
    if norm(A' - A) > 1e-10
        trA = A[1, 1] + A[2, 2] + A[3, 3]
        e0 = 0.5 * sqrt(trA + 1)
        e = (1 / (2 * sqrt(trA + 1))) * [A[3, 2] - A[2, 3], A[1, 3] - A[3, 1], A[2, 1] - A[1, 2]]
    end
    if norm(A' - A) ≤ 1e-10
        e0 = 1
        e = zeros(3)
    end        
    p = [e0; e]

    # Evaluate q
    q = [rO; p]
    return q, p, A
end

function KEval(apr, b)
    # Assuming `atil` function exists or is replaced by an equivalent Julia function
    K = 2 * vcat(hcat(dot(apr, b),     apr' * atil(b)),
             hcat(atil(apr) * b,   apr * b' + b * apr' - dot(apr, b) * I))
    return K
end

function MEval(q, SMDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    M = zeros(ngc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix

    for i in 1:nb
        m = SMDT[1, i]
        J = Diagonal([SMDT[2, i], SMDT[3, i], SMDT[4, i]])
        r, p = qPart(q, i)
        G = GEval(p)
        Mi = BlockDiagonal([m * I3, 4 * G' * J * G])  # blockdiag function to create block diagonal matrix
        M = add_constraint!(M, Mi, 7 * (i - 1), 7 * (i - 1))
    end

    return M
end

function M2Eval(q, mu, SMDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    M2 = zeros(ngc, ngc)

    for i in 1:nb
        r, p = qPart(q, i)
        mur, mup = qPart(mu, i)
        m = SMDT[1, i]
        J = Diagonal([SMDT[2, i], SMDT[3, i], SMDT[4, i]])
        G = GEval(p)
        Gmu = GEval(mup)
        M2i = vcat(zeros(3, 7), hcat(zeros(4, 3), TEval(4 * J * G * mup) - 4 * G' * J * Gmu))
        M2 = add_constraint!(M2, M2i, 7 * (i - 1), 7 * (i - 1))
    end

    return M2
end


function QAEval(tn, q, qd, SMDT, STSDAT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    uz = [0; 0; 1]
    uy = [0; 1; 0]

    QA = zeros(ngc)

    # Account for gravitational force in negative y direction
    for i in 1:nb
        mi = SMDT[1, i]
        QAGi = vcat(-mi * g * uz, zeros(4))
        QA = add_constraint!(QA, QAGi, 7 * (i - 1), 0)
    end

    # Account for TSDA forces
    for T in 1:NTSDA
        i, j, s1pr, s2pr, K, C, el0, F = STSDATPart(STSDAT, T)
        r1, p1 = qPart(q, i)
        r1d, p1d = qPart(qd, i)
        r2, p2, r2d, p2d = [0; 0; 0], [1; 0; 0; 0], [0; 0; 0], zeros(4)

        if j >= 1
            r2, p2 = qPart(q, j)
            r2d, p2d = qPart(qd, j)
        end

        A1 = ATran(p1)
        A2 = ATran(p2)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        BT1 = BTran(p1, s1pr)
        BT2 = BTran(p2, s2pr)
        el = sqrt(dot(d12, d12))
        eld = (1 / el) * dot(d12, r2d + BT2 * p2d - r1d - BT1 * p1d)
        f = K * (el - el0) + C * eld + F  # User insert F(el,eld) if needed
        QA1 = (f / el) * vcat(d12, BT1' * d12)
        QA = add_constraint!(QA, QA1, 7 * (i - 1), 0)

        if j >= 1
            QA2 = -(f / el) * vcat(d12, BT2' * d12)
            QA = add_constraint!(QA, QA2, 7 * (j - 1), 0)
        end
    end

    return QA
end


function QAsqqd(tn, q, qd, SMDT, STSDAT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    #uy = [0, 0, 1]

    QAsq = zeros(ngc, ngc)
    QAsqd = zeros(ngc, ngc)

    I3 = Matrix{Float64}(I, 3, 3)
    Z3 = zeros(3, 1)

    T = 1
    while T <= NTSDA
        # Evaluate QAsq
        i, j, s1pr, s2pr, K, C, el0, F = STSDATPart(STSDAT, T)
        r1, p1 = qPart(q, i)
        r1d, p1d = qPart(qd, i)
        r2 = zeros(3)
        p2 = [1, 0, 0, 0]
        r2d = zeros(3)
        p2d = zeros(4)
        if j >= 1
            r2, p2 = qPart(q, j)
            r2d, p2d = qPart(qd, j)
        end
        A1 = ATran(p1)
        A2 = ATran(p2)
        BT1 = BTran(p1, s1pr)
        BT2 = BTran(p2, s2pr)
        BT1d = BTran(p1d, s1pr)
        BT2d = BTran(p2d, s2pr)
        d12 = r2 + A2 * s2pr - r1 - A1 * s1pr
        el = sqrt(d12' * d12)
        a = r2d + BT2 * p2d - r1d - BT1 * p1d
        eld = (1 / el) * d12' * a
        f = K * (el - el0) + C * eld + F  # User insert F(el,eld) if needed
        # ... continue with the rest of the code
        # Remaining calculations
        elsq1 = -(1 / el) * d12' * hcat(I3, BT1)
        eldsq1 = -(1 / el^2) * d12' * a * elsq1 - (1 / el) * a' * hcat(I3, BT1) + (1 / el) * d12' * hcat(Z3, BT1d)
        b1 = (1 / el) * (K * elsq1 + C * eldsq1) - (f / el^2) * elsq1
        Q1sq1 = vcat(d12, BT1' * d12) * b1 + (f / el) * vcat(hcat(-I3, BT1'), hcat(BT1', BT1' * BT1 + KEval(s1pr, d12)))
        QAsq = add_constraint!(QAsq, Q1sq1, 7 * (i - 1), 7 * (i - 1))

        if j >= 1
            elsq2 = (1 / el) * d12' * hcat(I3, BT2)
            eldsq2 = -(1 / el^2) * d12' * a * elsq2 + (1 / el) * a' * hcat(I3, BT2) + (1 / el) * d12' * hcat(Z3, BT2d)
            b2 = (1 / el) * (K * elsq2 + C * eldsq2) - (f / el^2) * elsq2
            Q1sq2 = vcat(d12, BT1'* d12)  * b2 + (f / el) * vcat(hcat(I3, BT2), BT1'* hcat(I3, BT2))
            QAsq = add_constraint!(QAsq, Q1sq2, 7 * (i - 1), 7 * (j - 1))
            Q2sq1 = -vcat(d12, BT2' * d12) * b1 - (f / el) * vcat(hcat(-I3, BT1'), BT2'*hcat(-I3, -BT1'))
            QAsq = add_constraint!(QAsq, Q2sq1, 7 * (j - 1), 7 * (i - 1))
            Q2sq2 = -vcat(d12, BT2' * d12) * b2 - (f / el) * vcat(hcat(I3, BT2), hcat(BT2', BT2' * BT2 + KEval(s2pr, d12)))
            QAsq = add_constraint!(QAsq, Q2sq2, 7 * (j - 1), 7 * (j - 1))
        end

        # Evaluate QAsqd
        eldsq1d = (1 / el) * d12' * hcat(-I3, -BT1)
        Q1sq1d = (C / el) * vcat(d12, BT1' * d12) * eldsq1d
        QAsqd = add_constraint!(QAsqd, Q1sq1d, 7 * (i - 1), 7 * (i - 1))

        if j >= 1
            eldsq2d = (1 / el) * d12' * hcat(I3, BT2)
            Q1sq2d = (C / el) * vcat(d12, BT1' * d12) * eldsq2d
            QAsqd = Aadd_constraint!dd(QAsqd, Q1sq2d, 7 * (i - 1), 7 * (j - 1))
            Q2sq1d = -(C / el) * vcat(d12, BT2' * d12) * eldsq1d
            QAsqd = add_constraint!(QAsqd, Q2sq1d, 7 * (j - 1), 7 * (i - 1))
            Q2sq2d = -(C / el) * vcat(d12, BT2' * d12) * eldsq2d
            QAsqd = add_constraint!(QAsqd, Q2sq2d, 7 * (j - 1), 7 * (j - 1))
        end

        T = T + 1
    end

    return QAsq, QAsqd
end


function SEval(q, qd, SMDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    S = zeros(ngc)
    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix in Julia
    for i in 1:nb
        m = SMDT[1, i]
        J = diagm([SMDT[2, i], SMDT[3, i], SMDT[4, i]])
        r, p = qPart(q, i)
        rd, pd = qPart(qd, i)
        G = GEval(p)
        Gd = GEval(pd)
        Si = vcat([0; 0; 0], 8 * Gd' * J * Gd * p)  # Adjust if necessary for correct matrix/vector sizes
        # Add Si to S (Assuming 'Add' is defined elsewhere)
        S = add_constraint!(S, Si, 7 * (i - 1), 0)  # Adjust this line if the 'Add' function works differently in Julia
    end
    return S
end


function SqqdEval(q, qd, SMDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    Ssq = zeros(ngc, ngc)
    Ssqd = zeros(ngc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix of size 3x3

    for i = 1:nb
        m = SMDT[1, i]
        J = Diagonal([SMDT[2, i], SMDT[3, i], SMDT[4, i]])
        r, p = qPart(q, i)
        rd, pd = qPart(qd, i)
        G = GEval(p)
        Gd = GEval(pd)
        Ssqi = vcat(zeros(3, 7), hcat(zeros(4, 3),8 * Gd' * J * Gd))
        Ssq = add_constraint!(Ssq, Ssqi, 7 * (i - 1), 7 * (i - 1))
        Ssqdi = vcat(zeros(3, 7), hcat(zeros(4, 3),8 * TEval(J * Gd * p) - 8 * Gd' * J * G))
        Ssqd = add_constraint!(Ssqd, Ssqdi, 7 * (i - 1), 7 * (i - 1))
    end

    return Ssq, Ssqd
end


function TEval(a)
    T =vcat( hcat(0,-a'),hcat( a,-atil(a)))
    return T
end


function ODEfunct(tn, q, qd, SMDT, STSDAT, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    Gam = GamEval(tn, q, qd, SJDT, par)
    QA = QAEval(tn, q, qd, SMDT, STSDAT, par)
    S = SEval(q, qd, SMDT, par)
    RHS = vcat(QA + S, -Gam)
    M = MEval(q, SMDT, par)
    Phiq = PhiqEval(tn, q, SJDT, par)
    E = vcat(hcat(M,Phiq'), hcat(Phiq,zeros(nc, nc)))
    ECond = cond(E)
    println("cond",ECond)
    println("r",rank(E))
    # x = E \ RHS
    x = pinv(E) * RHS

    qdd = zeros(ngc)
    for i in 1:ngc
        qdd[i] = x[i]
    end
    Lam = zeros(nc)
    for i in 1:nc
        Lam[i] = x[ngc + i]
    end

    return qdd, Lam, ECond
end


function P3Eval(tn, q, qd, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    P3 = zeros(nc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)  # Identity matrix of size 3x3
    k = 1   # Joint number
    m = 0   # Constraint counter - 1

    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            P31, P32 = bbP3dist(i, j, s1pr, s2pr, d, tn, q, qd, par)
            P3 = add_constraint!(P3, P31, m, 7 * (i - 1))

            if j >= 1
                P3 = add_constraint!(P3, P32, m, 7 * (j - 1))
            end

            m = m + 1
        end

        # Spherical Constraint
        if SJDT[1, k] == 2
            i, j, s1pr, s2pr = SphPart(k, SJDT)
            # P3 spherical contribution is zero, no nonzero terms to add

            m = m + 3
        end

        # Rev Constraint
        if SJDT[1, k] == 4
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = RevPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            P3cyl11, P3cyl12 = bbP3dot2(i, j, ux2pr, s1pr, s2pr, tn, q, qd, par)
            P3cyl21, P3cyl22 = bbP3dot2(i, j, uy2pr, s1pr, s2pr, tn, q, qd, par)
            P3cyl31, P3cyl32 = bbP3dot1(i, j, uz1pr, ux2pr, tn, q, qd, par)
            P3cyl41, P3cyl42 = bbP3dot1(i, j, uz1pr, uy2pr, tn, q, qd, par)
            P3rev51, P3rev52 = bbP3dot2(i, j, uz2pr, s1pr, s2pr, tn, q, qd, par)
        
            P31k = [P3cyl11; P3cyl21; P3cyl31; P3cyl41; P3rev51]
            P32k = [P3cyl12; P3cyl22; P3cyl32; P3cyl42; P3rev52]
        
            P3 = add_constraint!(P3, P31k, m, 7 * (i - 1))
        
            if j >= 1
                P3 = add_constraint!(P3, P32k, m, 7 * (j - 1))
            end
        
            m += 5
        end
        # Translational Constraint
        if SJDT[1, k] == 5
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = TranPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            P3cyl11, P3cyl12 = bbP3dot2(i, j, ux2pr, s1pr, s2pr, tn, q, qd, par)
            P3cyl21, P3cyl22 = bbP3dot2(i, j, uy2pr, s1pr, s2pr, tn, q, qd, par)
            P3cyl31, P3cyl32 = bbP3dot1(i, j, uz1pr, ux2pr, tn, q, qd, par)
            P3cyl41, P3cyl42 = bbP3dot1(i, j, uz1pr, uy2pr, tn, q, qd, par)
            P3tran51, P3tran52 = bbP3dot1(i, j, uy1pr, ux2pr, tn, q, qd, par)
        
            P31k = [P3cyl11; P3cyl21; P3cyl31; P3cyl41; P3tran51]
            P32k = [P3cyl12; P3cyl22; P3cyl32; P3cyl42; P3tran52]
        
            P3 = add_constraint!(P3, P31k, m, 7 * (i - 1))
        
            if j >= 1
                P3 = add_constraint!(P3, P32k, m, 7 * (j - 1))
            end
        
            m += 5
        end
        # fxc Constraint
        if SJDT[1, k] == 1010
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            P31, P32 = bbP3fxc(i, j, s1pr, s2pr, d, tn, q, qd, par)
            P3 = add_constraint!(P3, P31, m, 7 * (i - 1))

            if j >= 1
                P3 = add_constraint!(P3, P32, m, 7 * (j - 1))
            end

            m = m + 1
        end                
        if SJDT[1, k] == 1020
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            P31, P32 = bbP3fy(i, j, s1pr, s2pr, d, tn, q, qd, par)
            P3 = add_constraint!(P3, P31, m, 7 * (i - 1))

            if j >= 1
                P3 = add_constraint!(P3, P32, m, 7 * (j - 1))
            end

            m = m + 1
        end           
        k = k + 1
    end

    return P3
end

function P2Eval(tn,q,qd,SJDT,par)
    # x=qd
    # Retrieve parameters from 'par'
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

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
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

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
            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = CylPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            P2cyl11, P2cyl12 = bbP2dot2(i, j, vx2pr, s1pr, s2pr,tn, q, qd, par)
            P2cyl21, P2cyl22 = bbP2dot2(i, j, vy2pr, s1pr, s2pr,tn, q, qd, par)
            P2cyl31, P2cyl32 = bbP2dot1(i, j, vz1pr, vx2pr, q,tn, qd, par)
            P2cyl41, P2cyl42 = bbP2dot1(i, j, vz1pr, vy2pr, q,tn, qd, par)

            P21k = vcat(P2cyl11, P2cyl21, P2cyl31, P2cyl41)
            P22k = vcat(P2cyl12, P2cyl22, P2cyl32, P2cyl42)

            P2 = add_constraint!(P2, P21k, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22k, m, 7 * (j - 1))
            end
            m += 4
        # Revolute Constraint
        elseif constraintType == 4
            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = RevPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            P2cyl11, P2cyl12 = bbP2dot2(i, j, vx2pr, s1pr, s2pr,tn, q, qd, par)
            P2cyl21, P2cyl22 = bbP2dot2(i, j, vy2pr, s1pr, s2pr,tn, q, qd, par)
            P2cyl31, P2cyl32 = bbP2dot1(i, j, vz1pr, vx2pr,tn, q, qd, par)
            P2cyl41, P2cyl42 = bbP2dot1(i, j, vz1pr, vy2pr,tn, q, qd, par)
            P2rev51, P2rev52 = bbP2dot2(i, j, vz2pr, s1pr, s2pr,tn, q, qd, par)

            P21k = vcat(P2cyl11, P2cyl21, P2cyl31, P2cyl41, P2rev51)
            P22k = vcat(P2cyl12, P2cyl22, P2cyl32, P2cyl42, P2rev52)

            P2 = add_constraint!(P2, P21k, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22k, m, 7 * (j - 1))
            end
            m += 5

        # Translational Constraint
        elseif constraintType == 5
            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = TranPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            P2cyl11, P2cyl12 = bbP2dot2(i, j, vx2pr, s1pr, s2pr,tn, q, qd, par)
            P2cyl21, P2cyl22 = bbP2dot2(i, j, vy2pr, s1pr, s2pr,tn, q, qd, par)
            P2cyl31, P2cyl32 = bbP2dot1(i, j, vz1pr, vx2pr,tn, q, qd, par)
            P2cyl41, P2cyl42 = bbP2dot1(i, j, vz1pr, vy2pr,tn, q, qd, par)
            P2tran51, P2tran52 = bbP2dot1(i, j, vy1pr, vx2pr,tn, q, qd, par)

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
        elseif constraintType == 20

            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = TranPart(k, SJDT)
            vy1pr = atil(vz1pr) * vx1pr
            vy2pr = atil(vz2pr) * vx2pr

            P2cyl11, P2cyl12 = bbP2dot2(i, j, vx2pr, s1pr, s2pr,tn, q, qd, par)
            P2cyl21, P2cyl22 = bbP2dot2(i, j, vy2pr, s1pr, s2pr,tn, q, qd, par)
            P2cyl31, P2cyl32 = bbP2dot1(i, j, vz1pr, vx2pr,tn, q, qd, par)
            P2cyl41, P2cyl42 = bbP2dot1(i, j, vz1pr, vy2pr,tn, q, qd, par)
            P2rev51, P2rev52 = bbP2dot2(i, j, vz2pr, s1pr, s2pr,tn, q, qd, par)
            P2tran61, P2tran62 = bbP2dot1(i, j, vy1pr, vx2pr,tn, q, qd, par)

            P21k = vcat(P2cyl11, P2cyl21, P2cyl31, P2cyl41,P2rev51, P2tran61)
            P22k = vcat(P2cyl12, P2cyl22, P2cyl32, P2cyl42,P2rev52, P2tran62)

            P2 = add_constraint!(P2, P21k, m, 7 * (i - 1))
            if j >= 1
                P2 = add_constraint!(P2, P22k, m, 7 * (j - 1))
            end
            m += 6           
        elseif constraintType == 1010  
            # Check if the constraint type is a Distance Constraint
            # Extract parameters for the Distance Constraint
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

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
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

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
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

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
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

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
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

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

function P4Eval(tn, q, eta, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    P4 = zeros(ngc, ngc)
    I3 = Matrix{Float64}(I, 3, 3)
    Z1 = zeros(3, 1)
    Z3 = zeros(3, 3)
    k = 1    # Joint number
    m = 0    # Address in vector eta

    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            etak = eta[m + 1]

            P411, P412, P422 = bbP4dist(i, j, s1pr, s2pr, d, tn, q, etak, par)
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))

            if j >= 1
                P421 = transpose(P412)
                P4 = add_constraint!(P4, P412, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 1
        end

        # Spherical Constraint
        if SJDT[1, k] == 2
            i, j, s1pr, s2pr = SphPart(k, SJDT)
            etak = eta[m + 1 : m + 3]

            P411, P412, P422 = bbP4sph(i, j, s1pr, s2pr, tn, q, etak, par)    
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))    

            if j >= 1
                # P412 and P421 are zero
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 3
        end
        # Rev
        if SJDT[1, k] == 4
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = RevPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            P4cyl111, P4cyl112, P4cyl122 = bbP4dot2(i, j, ux2pr, s1pr, s2pr, tn, q, eta[m + 1], par)
            P4cyl211, P4cyl212, P4cyl222 = bbP4dot2(i, j, uy2pr, s1pr, s2pr, tn, q, eta[m + 2], par)
            P4cyl311, P4cyl312, P4cyl322 = bbP4dot1(i, j, uz1pr, ux2pr, tn, q, eta[m + 3], par)
            P4cyl411, P4cyl412, P4cyl422 = bbP4dot1(i, j, uz1pr, uy2pr, tn, q, eta[m + 4], par)
            P4rev511, P4rev512, P4rev522 = bbP4dot2(i, j, uz2pr, s1pr, s2pr, tn, q, eta[m + 5], par)
        
            P411k = P4cyl111 + P4cyl211 + P4cyl311 + P4cyl411 + P4rev511
            P412k = P4cyl112 + P4cyl212 + P4cyl312 + P4cyl412 + P4rev512
            P421k = transpose(P412k)
            P422k = P4cyl122 + P4cyl222 + P4cyl322 + P4cyl422 + P4rev522
        
            P4 = add_constraint!(P4, P411k, 7 * (i - 1), 7 * (i - 1))
        
            if j >= 1
                P4 = add_constraint!(P4, P412k, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421k, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422k, 7 * (j - 1), 7 * (j - 1))
            end
        
            m += 5
        end
        # Tra
        if SJDT[1, k] == 5
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = TranPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            P4cyl111, P4cyl112, P4cyl122 = bbP4dot2(i, j, ux2pr, s1pr, s2pr, tn, q, eta[m + 1], par)
            P4cyl211, P4cyl212, P4cyl222 = bbP4dot2(i, j, uy2pr, s1pr, s2pr, tn, q, eta[m + 2], par)
            P4cyl311, P4cyl312, P4cyl322 = bbP4dot1(i, j, uz1pr, ux2pr, tn, q, eta[m + 3], par)
            P4cyl411, P4cyl412, P4cyl422 = bbP4dot1(i, j, uz1pr, uy2pr, tn, q, eta[m + 4], par)
            P4tran511, P4tran512, P4tran522 = bbP4dot1(i, j, uy1pr, ux2pr, tn, q, eta[m + 3], par)
        
            P411k = P4cyl111 + P4cyl211 + P4cyl311 + P4cyl411 + P4tran511
            P412k = P4cyl112 + P4cyl212 + P4cyl312 + P4cyl412 + P4tran512
            P421k = transpose(P412k)
            P422k = P4cyl122 + P4cyl222 + P4cyl322 + P4cyl422 + P4tran522
        
            P4 = add_constraint!(P4, P411k, 7 * (i - 1), 7 * (i - 1))
        
            if j >= 1
                P4 = add_constraint!(P4, P412k, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421k, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422k, 7 * (j - 1), 7 * (j - 1))
            end
        
            m += 5
        end
        # fxc Constraint
        if SJDT[1, k] == 1010
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            etak = eta[m + 1]

            P411, P412, P422 = bbP4fxc(i, j, s1pr, s2pr, d, tn, q, etak, par)
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))

            if j >= 1
                P421 = transpose(P412)
                P4 = add_constraint!(P4, P412, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 1
        end        
        if SJDT[1, k] == 1020
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            etak = eta[m + 1]

            P411, P412, P422 = bbP4fy(i, j, s1pr, s2pr, d, tn, q, etak, par)
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))

            if j >= 1
                P421 = transpose(P412)
                P4 = add_constraint!(P4, P412, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 1
        end 
        if SJDT[1, k] == 1030
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            etak = eta[m + 1]

            P411, P412, P422 = bbP4fz(i, j, s1pr, s2pr, d, tn, q, etak, par)
            P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))

            if j >= 1
                P421 = transpose(P412)
                P4 = add_constraint!(P4, P412, 7 * (i - 1), 7 * (j - 1))
                P4 = add_constraint!(P4, P421, 7 * (j - 1), 7 * (i - 1))
                P4 = add_constraint!(P4, P422, 7 * (j - 1), 7 * (j - 1))
            end

            m += 1
        end       
        k += 1
    end

    # Euler Parameter Normalization Constraint
    i = 1
    while i <= nb
        etak = eta[m + 1]
        P411 = etak * vcat(zeros(3, 7), hcat(zeros(4, 3),I))
        P4 = add_constraint!(P4, P411, 7 * (i - 1), 7 * (i - 1))
        i += 1
        m += 1
    end

    return P4
end

function P5Eval(tn, q, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)

    # Enter Constraint t derivatives of P(t, q, qd, par); Default is Zeros
    Pst = zeros(nc)
    Pstt = zeros(nc)
    Pstq = zeros(nc, ngc)
    Psttq = zeros(nc, ngc)

    return Pst, Pstt, Pstq, Psttq
end

function computeDriverInputs(u1pr, v1pr, u2pr, A1, theta, thetad, thetadd, j, q)
    c = s = PD = PDd = PDdd = 0
    if j == 0
        c = dot(u1pr, A1' * u2pr)
        s = dot(v1pr, A1' * u2pr)
    else
        r2, p2 = qPart(q, j)
        A2 = ATran(p2)
        c = dot(u1pr, A1' * A2 * u2pr)
        s = dot(v1pr, A1' * A2 * u2pr)
    end

    if abs(c) >= abs(s)
        PD = -sin(theta)
        PDd = -thetad * cos(theta)
        PDdd = -thetadd * cos(theta) + thetad^2 * sin(theta)
    else
        PD = -cos(theta)
        PDd = thetad * sin(theta)
        PDdd = thetadd * sin(theta) + thetad^2 * cos(theta)
    end

    return PD, PDd, PDdd
end

function PhiEval(tn, q, SJDT, par)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    Phi = zeros(nc)
    I3 = Matrix{Float64}(I, 3, 3)
    k = 1        # Joint No.
    m = 0        # Constraint Counter - 1
    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
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
            i, j, s1pr, s2pr, u1pr, vz1pr, u2pr, vz2pr = CylPart(k, SJDT)
            v1pr = atil(vz1pr) * u1pr
            vy2pr = atil(vz2pr) * u2pr
            Phicyl1 = bbPhidot2(i, j, u2pr, s1pr, s2pr,tn, q, par)
            Phicyl2 = bbPhidot2(i, j, vy2pr, s1pr, s2pr,tn, q, par)
            Phicyl3 = bbPhidot1(i, j, vz1pr, u2pr,tn, q, par)
            Phicyl4 = bbPhidot1(i, j, vz1pr, vy2pr,tn, q, par)
            Phik = vcat(Phicyl1, Phicyl2, Phicyl3, Phicyl4)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 4
        end
        # Revolute Constraint
        if SJDT[1, k] == 4
            i, j, s1pr, s2pr, u1pr, vz1pr, u2pr, vz2pr = RevPart(k, SJDT)
            v1pr = atil(vz1pr) * u1pr
            vy2pr = atil(vz2pr) * u2pr

            Phicyl1 = bbPhidot2(i, j, u2pr, s1pr, s2pr,tn, q, par)
            Phicyl2 = bbPhidot2(i, j, vy2pr, s1pr, s2pr,tn, q, par)
            Phicyl3 = bbPhidot1(i, j, vz1pr, u2pr,tn, q, par)
            Phicyl4 = bbPhidot1(i, j, vz1pr, vy2pr,tn, q, par)
            Phirev5 = bbPhidot2(i, j, vz2pr, s1pr, s2pr,tn, q, par)

            Phik = vcat(Phicyl1, Phicyl2, Phicyl3, Phicyl4, Phirev5)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 5
        end
        # Check if it's a translational constraint
        if SJDT[1, k] == 5
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = TranPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr

            Phicyl1 = bbPhidot2(i, j, ux2pr, s1pr, s2pr, tn, q, par)
            Phicyl2 = bbPhidot2(i, j, uy2pr, s1pr, s2pr, tn, q, par)
            Phicyl3 = bbPhidot1(i, j, uz1pr, ux2pr, tn, q, par)
            Phicyl4 = bbPhidot1(i, j, uz1pr, uy2pr, tn, q, par)
            Phitran5 = bbPhidot1(i, j, uy1pr, ux2pr, tn, q, par)

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
        if SJDT[1, k] == 20# fix body
            
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = TranPart(k, SJDT)

            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr

            Phicyl1 = bbPhidot2(i, j, ux2pr, s1pr, s2pr, tn, q, par)
            Phicyl2 = bbPhidot2(i, j, uy2pr, s1pr, s2pr, tn, q, par)
            Phicyl3 = bbPhidot1(i, j, uz1pr, ux2pr, tn, q, par)
            Phicyl4 = bbPhidot1(i, j, uz1pr, uy2pr, tn, q, par)
            Phirev5 = bbPhidot2(i, j, uz2pr, s1pr, s2pr, tn, q, par)
            Phitran6 = bbPhidot1(i, j, uy1pr, ux2pr, tn, q, par)

            Phik = vcat(Phicyl1, Phicyl2, Phicyl3, Phicyl4, Phirev5, Phitran6)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 6
        end
        # fxc Constraint
        if SJDT[1, k] == 1010
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            Phik = bbPhifxc(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end
        if SJDT[1, k] == 1020
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            Phik = bbPhify(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end
        if SJDT[1, k] == 1030
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            Phik = bbPhifz(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end

        if SJDT[1, k] == 1050
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
            Phik = bbPhif(i, j, s1pr, s2pr, d,tn, q, par)
            Phi = add_constraint!(Phi, Phik, m, 0)
            m += 1
        end
        if SJDT[1, k] == 1060
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)
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

function PhiqEval(tn, q, SJDT, par)

    nb,ngc,nh,nc,g,intol,Atol,h0,hvar,NTSDA = parPart(par)
    Phiq = zeros(Float64,nc, ngc)

    I3 = Matrix{Float64}(I, 3, 3)
    k = 1  # Joint number
    m = 0  # Constraint equation counter - 1

    while k <= nh
        # Distance Constraint
        if SJDT[1, k] == 1
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqdist(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        # Spherical Constraint
        if SJDT[1, k] == 2
            i, j, s1pr, s2pr = SphPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqsph(i, j, s1pr, s2pr, tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7 * (i - 1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7 * (j - 1))
            end

            m += 3
        end

        # Revolute Constraint
        if SJDT[1, k] == 4
            i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr = RevPart(k, SJDT)
            vy1pr = cross(vz1pr, vx1pr)
            vy2pr = cross(vz2pr, vx2pr)
        
            Phiqcyl11, Phiqcyl12 = bbPhiqdot2(i, j, vx2pr, s1pr, s2pr,tn, q, par)
            Phiqcyl21, Phiqcyl22 = bbPhiqdot2(i, j, vy2pr, s1pr, s2pr,tn, q, par)
            Phiqcyl31, Phiqcyl32 = bbPhiqdot1(i, j, vz1pr, vx2pr,tn, q, par)
            Phiqcyl41, Phiqcyl42 = bbPhiqdot1(i, j, vz1pr, vy2pr,tn, q, par)
            Phiqrev51, Phiqrev52 = bbPhiqdot2(i, j, vz2pr, s1pr, s2pr,tn, q, par)
        
            Phiq1k = vcat(Phiqcyl11, Phiqcyl21, Phiqcyl31, Phiqcyl41, Phiqrev51)
            Phiq2k = vcat(Phiqcyl12, Phiqcyl22, Phiqcyl32, Phiqcyl42, Phiqrev52)
        
            Phiq = add_constraint!(Phiq, Phiq1k, m, 7 * (i - 1))
        
            if j ≥ 1
                Phiq = add_constraint!(Phiq, Phiq2k, m, 7 * (j - 1))
            end
        
            m += 5
        end

        if SJDT[1, k] == 5
            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = TranPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            Phiqcyl11, Phiqcyl12 = bbPhiqdot2(i, j, ux2pr, s1pr, s2pr, tn, q, par)
            Phiqcyl21, Phiqcyl22 = bbPhiqdot2(i, j, uy2pr, s1pr, s2pr, tn, q, par)
            Phiqcyl31, Phiqcyl32 = bbPhiqdot1(i, j, uz1pr, ux2pr, tn, q, par)
            Phiqcyl41, Phiqcyl42 = bbPhiqdot1(i, j, uz1pr, uy2pr, tn, q, par)
            Phiqtran51, Phiqtran52 = bbPhiqdot1(i, j, uy1pr, ux2pr, tn, q, par)
        
            Phiq1k = vcat(Phiqcyl11, Phiqcyl21, Phiqcyl31, Phiqcyl41, Phiqtran51)
            Phiq2k = vcat(Phiqcyl12, Phiqcyl22, Phiqcyl32, Phiqcyl42, Phiqtran52)
        
            Phiq = add_constraint!(Phiq, Phiq1k, m, 7 * (i - 1))
        
            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2k, m, 7 * (j - 1))
            end
        
            m += 5
        end
        

        # Rotation Driver
        if SJDT[1, k] == 10
            i, j, vx1pr, vy1pr, vx2pr = RotDrPart(k, SJDT)
            Phiq1, Phiq2 = bbPhiqRotDr(i, j, vx1pr, vy1pr, vx2pr,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7 * (i - 1))
        
            if j ≥ 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7 * (j - 1))
            end
            m += 1
        end

                # Revolute Constraint

        # Revolute Constraint
        if SJDT[1, k] == 20

            i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr = TranPart(k, SJDT)
            uy1pr = atil(uz1pr) * ux1pr
            uy2pr = atil(uz2pr) * ux2pr
        
            Phiqcyl11, Phiqcyl12 = bbPhiqdot2(i, j, ux2pr, s1pr, s2pr, tn, q, par)
            Phiqcyl21, Phiqcyl22 = bbPhiqdot2(i, j, uy2pr, s1pr, s2pr, tn, q, par)
            Phiqcyl31, Phiqcyl32 = bbPhiqdot1(i, j, uz1pr, ux2pr, tn, q, par)
            Phiqcyl41, Phiqcyl42 = bbPhiqdot1(i, j, uz1pr, uy2pr, tn, q, par)
            Phiqrev51, Phiqrev52 = bbPhiqdot2(i, j, uz2pr, s1pr, s2pr,tn, q, par)
            Phiqtran61, Phiqtran62 = bbPhiqdot1(i, j, uy1pr, ux2pr, tn, q, par)
        
            Phiq1k = vcat(Phiqcyl11, Phiqcyl21, Phiqcyl31, Phiqcyl41,Phiqrev51, Phiqtran61)
            Phiq2k = vcat(Phiqcyl12, Phiqcyl22, Phiqcyl32, Phiqcyl42,Phiqrev52, Phiqtran62)
        
            Phiq = add_constraint!(Phiq, Phiq1k, m, 7 * (i - 1))
        
            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2k, m, 7 * (j - 1))
            end
        
            m += 6
        end

        # ... (Continue with the rest of the constraints as per the given code)
        # fxc Constraint
        if SJDT[1, k] == 1010
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqfxc(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        if SJDT[1, k] == 1020
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqfy(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        if SJDT[1, k] == 1030
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqfz(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end

        if SJDT[1, k] == 1050
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiqf(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        if SJDT[1, k] == 1060
            i, j, s1pr, s2pr, d = DistPart(k, SJDT)

            Phiq1, Phiq2 = bbPhiq_Poly(i, j, s1pr, s2pr, d,tn, q, par)
            Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))

            if j >= 1
                Phiq = add_constraint!(Phiq, Phiq2, m, 7*(j-1))
            end

            m += 1
        end
        k += 1
    end

    # Euler Parameter Normalization Constraints
    i = 1
    while i <= nb
        r1, p1 = qPart(q, i)
        Phiq1 = hcat(zeros(1, 3),p1')
        Phiq = add_constraint!(Phiq, Phiq1, m, 7*(i-1))
        m += 1
        i += 1
    end

    # Convert to DataFrame
    df = DataFrame(Phiq, :auto)
    # Write to CSV
    #CSV.write("Phiq.csv", df)

    return Phiq
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

    return i, j, s1pr, s2pr, vx1pr, vz1pr, vx2pr, vz2pr
end

function TranPart(k, SJDT)
    i = SJDT[2, k]
    j = SJDT[3, k]
    s1pr = [SJDT[4, k], SJDT[5, k], SJDT[6, k]]
    s2pr = [SJDT[7, k], SJDT[8, k], SJDT[9, k]]
    ux1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    uz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    ux2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]
    uz2pr = [SJDT[20, k], SJDT[21, k], SJDT[22, k]]

    return i, j, s1pr, s2pr, ux1pr, uz1pr, ux2pr, uz2pr
end


function RotDrPart(k, SJDT)
    # Data derived from the host revolute or cylindrical joint data table
    i = SJDT[2, k]
    j = SJDT[3, k]
    vx1pr = [SJDT[11, k], SJDT[12, k], SJDT[13, k]]
    vz1pr = [SJDT[14, k], SJDT[15, k], SJDT[16, k]]
    vy1pr = cross(vz1pr, vx1pr)  # using cross product to get the perpendicular vector
    vx2pr = [SJDT[17, k], SJDT[18, k], SJDT[19, k]]

    return i, j, vx1pr, vy1pr, vx2pr
end


end # module mathfunction