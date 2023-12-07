using LinearAlgebra
using Plots
function run()
    println("start")
    # Given parameters

    m1 = 30

    g = 9.81
    l1 = 50
    l2 = 220
    w = 5 * π
    dt = 0.001
    t = 0:dt:5
    n = length(t)
    q = zeros(5, n)
    q_v = zeros(5, n)
    q_ac = zeros(5, n)
    q[:, 1] = [1, 1, 0, 0, 0]
    q_v[:, 1] = zeros(5)
    PosConstrNorm = zeros(1, n)
    PosConstrNorm[1, 1] = phi(q[:, 1])
    for i in 1:n
        q1 = q[1, i]
        q2 = q[2, i]
        q3 = q[3, i]
        qv1 = q_v[1, i]
        qv2 = q_v[2, i]
        qv3 = q_v[3, i]

        # Choose one of the A matrix definitions based on your scenario
        # Sphere and plane y=0 constraint
        A = [m1 0 0 q1 0;
            0 m1 0 q2 q2;
            0 0 m1 q3 0;
            q1 q2 q3 0 0;
            0 q2 0 0 0]

        # Calculate Qa
        Qa = [0, 0, -m1 * g]

        # Sphere and plane y^2/2=0 constraint
        garma = [-qv1^2 - qv2^2 - qv3^2, -qv2^2]
        # garma = [-qv1^2 - qv2^2 - qv3^2]
        B = vcat(Qa, garma)

        temp = A \ B
        q_ac[1:5, i] = temp[1:5]

        println(i)
        if i == n
            break
        end

        if i == 1
            q_v[:, i+1] = q_v[:, i] + dt * q_ac[:, i]
            q[:, i+1] = q[:, i] + dt * q_v[:, i]
        else
            q_v[:, i+1] = q_v[:, i] + dt * (q_ac[:, i-1] + q_ac[:, i]) / 2
            q[:, i+1] = q[:, i] + dt * (q_v[:, i-1] + q_v[:, i]) / 2
        end
        phi1 = phi(q) 
        #phi2=q2^2/2
        PosConstrNorm[1, i+1] = norm([phi1])
    end
    println("end")
    gr()
    f = plot(t, q[1, :], lw=2, ls=:solid, color=:red, label="x")
    f = plot!(t, q[2, :], lw=2, ls=:solid, color=:green, label="y")
    f = plot!(t, q[3, :], lw=2, ls=:solid, color=:blue, label="z")
    display(f)
    f2 = f = plot(t, PosConstrNorm[1, :], lw=3, ls=:solid, color=:blue, label="PosConstrNorm")
    display(f2)
    println("end2")
end

function phi(q)
   return (q[1]^2 + q[2]^2 + q[3]^2) / 2
end

run()