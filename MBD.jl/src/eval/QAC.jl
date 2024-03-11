

function QACEval(tn, q, qd, SMDT, STSDAT, par, p_contact)
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    ld_damper, ld_contact = p_contact
    uz = [0; 0; 1]
    uy = [0; 1; 0]

    QAC = zeros(ngc)
    for d_damper in ld_damper
        # Account for gravitational force in negative y direction
        b=d_damper["b"]
        damp=d_damper["damp"]
            # mi = SMDT[1, i]
        index = 7 * (b - 1) + 3
        QACi = vcat(-qd[index] * damp, -qd[index] * damp, -qd[index] * damp, zeros(4))
        QAC = add_constraint!(QAC, QACi, 7 * (b - 1), 0)
    end
    # if q[3] < -0.5
    #     delta_v = q[3] + 0.5
    #     vel_v = qd[3]
    #     updater = create_updater()
    #     init_vel_v = updater(-delta_v, vel_v)
    #     fz = calculate_F_flore_Abs(abs(delta_v), abs(vel_v), checkSign(delta_v, vel_v), 1e9, 0.97, 0.6, init_vel_v)
    #     println("d$delta_v v$vel_v inv$init_vel_v f$fz")
    #     QACi = vcat(0, 0, fz, zeros(4))
    #     QAC = add_constraint!(QAC, QACi, 7 * (1 - 1), 0)
    # end
    for d_contact in ld_contact
        b=d_contact["b"]
        threshold=d_contact["pos"]
        index = 7 * (b - 1) + 3
        if b == 4
            println("b $b th $threshold ind $index")
        end
        if q[index] < threshold
            delta_v = q[index] - (threshold)
            vel_v = qd[index]
            updater = create_updater()
            init_vel_v = updater(-delta_v, vel_v)

            fz = calculate_F_flore_Abs(abs(delta_v), abs(vel_v), checkSign(delta_v, vel_v), 1e10, 0.97, 0.6, init_vel_v)
            if b == 4
                println("q[index]$(q[index]) d$delta_v v$vel_v inv$init_vel_v f$fz")
            end
            QACi = vcat(0, 0, fz, zeros(4))
            QAC = add_constraint!(QAC, QACi, 7 * (b - 1), 0)
        end
    end
    return QAC
end

function calculate_F_flore_Abs(q::Float64, v::Float64, flag_move_dir::Int, Eeq::Float64, F2::Float64, E2::Float64, init_vel::Float64)::Float64
    # if flag_move_dir == 0 || abs(init_vel) < 10.0
    #     init_vel = 0.0
    # elseif abs(v) > 1.1 * abs(init_vel) && abs(init_vel) > 10.0
    #     v = 1.1 * abs(init_vel) * flag_move_dir
    # else
    #     v = v * flag_move_dir
    # end
    v = v * flag_move_dir
    init_vel = -abs(init_vel)
    init_vel=0.0
    Fd = Eeq * (q /1000/ F2)^1.5
    if init_vel != 0.0
        Fd *= (1.0 + 8.0 * (1.0 - E2) * v / (5.0 * E2 * init_vel))
    end
    return abs(Fd)
end

function checkSign(a::Float64, b::Float64)::Int
    if a == 0 || b == 0
        return 0 # If either number is zero, return 0
    elseif a * b > 0
        return 1 # If the product is positive, they have the same sign
    else
        return -1 # If the product is negative, they have opposite signs
    end
end

function create_updater()
    initial_v = 0.0 # 初始值
    prev_q = 0.0    # 上一次的q值

    # 定义并返回闭包，该闭包接受q和v作为输入
    return (q, v) -> begin
        #println("pq $prev_q q $q")
        # 如果prev_q <= 0 且当前q > 0，则更新initial_v
        if prev_q <= 0.0 && q > 0.0
            # println("update iv")
            initial_v = v
        else
            # println("no update iv")
        end
        prev_q = q # 更新prev_q为当前的q值

        # 返回当前有效的initial_v
        return initial_v
    end
end