include("./contact.jl")
function QACEval(tn, q, qd, SMDT, STSDAT, par, p_contact)
    println("t=$tn----------------------------------------------------------------")
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    ld_damper, ld_contact,PUSH = p_contact
    uz = [0; 0; 1]
    uy = [0; 1; 0]
    QAC = zeros(ngc)
    # if tn>=0.1 && tn<=0.5
    F=0.0
    #println(PUSH)
    if PUSH!=0
        F=PUSH(tn)
        QACi = vcat(-F,0.,0., zeros(4))
        QAC = add_constraint!(QAC, QACi, 7 * (11 - 1), 0)
    end

    # end
    # QACi = vcat(0.,-300.,0., zeros(4))
    # QAC = add_constraint!(QAC, QACi, 7 * (1 - 1), 0)
    for d_damper in ld_damper
        # Account for gravitational force in negative y direction
        b=d_damper["b"]
        damp=d_damper["damp"]
            # mi = SMDT[1, i]
        index1 = 7 * (b - 1) + 1
        index2 = 7 * (b - 1) + 2
        index3 = 7 * (b - 1) + 3
        QACi = vcat(-qd[index1] * damp[1], -qd[index2] * damp[2], -qd[index3] * damp[3], zeros(4))
        QAC = add_constraint!(QAC, QACi, 7 * (b - 1), 0)
    end

    for d_contact in ld_contact
        type = d_contact["type"]
        QACi = vcat(0, 0, 0, zeros(4))
        if type == "guide"
            b = d_contact["b"]
            index = 7 * (b - 1) + 2
            # fy=calculate_Fy_prepare(d_contact,q,qd)
            fx, fy = calculate_contact_geo(d_contact, q, qd)
            QACi = vcat(fx, fy, 0, zeros(4))
        elseif type == "pos"
            fz = calculate_F_prepare(d_contact, q, qd)
            b = d_contact["b"]
            index = 7 * (b - 1) + 3
            # if b == 4
            #   println("q[index]$(q[index]) d$delta_v v$vel_v inv$init_vel_v f$fz")
            # end
            QACi = vcat(0, 0, fz, zeros(4))
        end

        QAC = add_constraint!(QAC, QACi, 7 * (b - 1), 0)
    end
    return QAC

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