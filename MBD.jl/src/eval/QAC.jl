include("./contact.jl")
function QACEval(tn, q, qd, SMDT, STSDAT, par, p_contact)
    #println("t=$tn-QACEval--------------------------------------------------------")
    nb, ngc, nh, nc, g, intol, Atol, h0, hvar, NTSDA = parPart(par)
    ld_damper, ld_contact,PUSH,start_time = p_contact
    uz = [0; 0; 1]
    uy = [0; 1; 0]
    QAC = zeros(ngc)

    F=0.0
    #println(PUSH)
    if PUSH!=0
        if tn>=start_time
            F=PUSH(tn-start_time)
        end
        QACi = vcat(-F,0.,0., zeros(4))
        QAC = add_constraint!(QAC, QACi, 7 * (11 - 1), 0) # to the handle
    end

    for d_damper in ld_damper
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
            fx, fy, debug_dict = calculate_contact_geo(d_contact, q, qd)
            QACi = vcat(fx, fy, 0, zeros(4))
        elseif type == "pos"
            fx, fy, fz,debugDict = calculate_F_plus(d_contact, q, qd)

            b = d_contact["b"]
            index = 7 * (b - 1) + 3
            QACi = vcat(fx, fy,  fz, zeros(4))
        elseif type == "simple_y"
            fy=0.0
            thres=0.00000
            b = d_contact["b"]
            index = 7 * (b - 1) + 2
            y=q[index]-0.02
            vel_v = qd[index]
            vel_v=0.0
            if y < (-thres)
                delta_v = y - (thres)

                Eeq = d_contact["p"]["Eeq"]
                n = d_contact["p"]["n_damper"]
                init_vel_v = d_contact["p"]["B5"]
                Cr=0.5
                eta = 8 * (1.0 - Cr) / (5.0 * Cr * init_vel_v);
                fy = Eeq * (-delta_v)^ (1.5) * (1.0 + eta * vel_v);

            end
            if y > (thres)
                delta_v = y - (thres)

                Eeq = d_contact["p"]["Eeq"]
                n = d_contact["p"]["n_damper"]
                init_vel_v = d_contact["p"]["B5"]
                Cr=0.5
                eta = 8.0 * (1.0 - Cr) / (5.0 * Cr * -init_vel_v);
                fy = -Eeq * (delta_v)^ (1.5) * (1.0 + eta * -vel_v);

            end
            #println("y=$(y) fy=$fy")
            index = 7 * (b - 1) + 3
            QACi = vcat(0.0, fy,  0.0, zeros(4))    
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