function calculate_F_prepare(d_contact,q,qd)
    b=d_contact["b"]
    threshold=d_contact["pos"]
    index = 7 * (b - 1) + 3
    # if b == 4
    #     println("b $b th $threshold ind $index")
    # end
    fz=0.0
    if q[index] < threshold
        delta_v = q[index] - (threshold)
        vel_v = qd[index]
        init_vel_v=-0.1
        # println("delta_v $delta_v vel_v $vel_v")
        # fz = calculate_F_flore_Abs(abs(delta_v), abs(vel_v), checkSign(delta_v, vel_v), 1e10, 0.97, 0.6, init_vel_v)
        fz = calculate_F_flore_Abs_modified((delta_v), (vel_v), checkSign(delta_v, vel_v), 1e9, 0.97, 0.6, init_vel_v)
    end
    # println("fz $fz")
    return fz>1e20 ? 1e20 : fz
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
    #init_vel=0.0
    Fd = Eeq * (q / F2)^1.5
    if init_vel != 0.0
        Fd *= (1.0 + 8.0 * (1.0 - E2) * v / (5.0 * E2 * init_vel))
    end
    return abs(Fd)
end

function calculate_F_flore_Abs_modified(q::Float64, v::Float64, flag_move_dir::Int, Eeq::Float64, F2::Float64, E2::Float64, init_vel::Float64)::Float64

    # Conditional calculation of Fd based on the sign of q
    if q < 0
        Fd = Eeq * (abs(q) / F2)^1.5
        #Fd *= (1.0 + 8.0 * (1.0 - E2) * v / (5.0 * E2 * (init_vel))) 
    else
        Fd = 0.0
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