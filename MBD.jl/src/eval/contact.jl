include("contact_geo.jl")

function calculate_F_plus(d_contact,q,qd)
    b=d_contact["b"]
    threshold=d_contact["pos"]
    spl = d_contact["guide"]
    Eeq = d_contact["p"]["Eeq"]
    n = d_contact["p"]["n_damper"]
    B5 = d_contact["p"]["B5"]
    fric_force_mul = d_contact["p"]["fric_force_mul"]
    
    index_x = 7 * (b - 1) + 1
    index_y = 7 * (b - 1) + 2
    index = 7 * (b - 1) + 3
    pos_part = [q[index_x], q[index_y]]
    vel = [qd[index_x], qd[index_y]]

    fz=0.0

    if q[index] < threshold
        delta_v = q[index] - (threshold)
        vel_v = qd[index]
        init_vel_v=B5
        fz = calculate_F_D_Abs((delta_v), (vel_v), checkSign(delta_v, vel_v), Eeq, n)
    end

    #######################################################################################
    # calculate_contact_geo
    #######################################################################################
    pen_pos_delta_local, r_init_vel, ub = calculate_delta_spline(pos_part, spl, 0)
    
    #######################################################################################
    # friction
    #######################################################################################
    fx=0.0
    fy=0.0
    fric=0.0
    vel_slide=0.0
    miu=0.0
    if d_contact["p"]["f_type"]=="Coulomb"
        v1=d_contact["p"]["f_v_sec"]["1"]
        v2=d_contact["p"]["f_v_sec"]["2"]
        v3=d_contact["p"]["f_v_sec"]["3"]
        v4=d_contact["p"]["f_v_sec"]["4"]

        flag = 0
        if ub[2] > 0.0
            flag = 1
        elseif ub[2] < 0.0
            flag = -1
        else
        end
        vel_slide = project_vector_and_return_value(vel, ub)
        miu=abs(smoothstep(abs(vel_slide), v1, v2))-0.1*abs(smoothstep(abs(vel_slide), v3, v4))
        fric=miu * sign(vel[2]) * ub * fz *0.3
        fx, fy = fric_force_mul * fric
    end
    if d_contact["p"]["f_type"]=="fix"
        start_v=d_contact["p"]["f_start_v"]
        flag = 0
        if ub[2] > 0.0
            flag = 1
        elseif ub[2] < 0.0
            flag = -1
        else
        end
        vel_slide = project_vector_and_return_value(vel, ub)
        miu=abs(smoothstep(abs(vel_slide), 0.0002, start_v))
        fric=miu * sign(vel[2]) * ub * fz *0.3
        fx, fy = fric_force_mul * fric
    end
    if d_contact["p"]["f_type"]=="no"
        start_v=d_contact["p"]["f_start_v"]
        flag = 0
        if ub[2] > 0.0
            flag = 1
        elseif ub[2] < 0.0
            flag = -1
        else
        end
        vel_slide = project_vector_and_return_value(vel, ub)
        miu=abs(smoothstep(abs(vel_slide), 0.0002, start_v))
    end
    debug_dir=0.0
    test=fy*vel[2]
    if test>0.0
        debug_dir=1.0
    end
    if test<0.0
        debug_dir=-1.0
    end

    debugDict=Dict("fx"=>fx, "fy"=>fy,"fz"=>fz,"ub_x"=>ub[1], "ub_y"=>ub[2],"vx"=>vel[1], "vy"=>vel[2],
    "vel_slide"=>vel_slide,"miu"=>miu,"debug_dir"=>debug_dir)


    return fx,fy,(fz>1e20 ? 1e20 : fz) ,debugDict
end

function calculate_Fy_prepare(d_contact,q,qd)
    b=d_contact["b"]
    guide=d_contact["guide"]
    index_x = 7 * (b - 1) + 1
    threshold=guide(q[index_x])
    index = 7 * (b - 1) + 2
    # if b == 4
    #     println("b $b th $threshold ind $index")
    # end
    fz=0.0
    # Calculate the difference between the current value and the threshold
    delta_v = q[index] - threshold
    println("x=$(q[index_x]) ysp=$threshold y=$(q[index]) delta_v $delta_v")
    # Define initial velocity
    init_vel_v=-0.1
    
    # Check if q[index] is less than the threshold
    if delta_v < 0
        vel_v = qd[index]
        fz = calculate_F_flore_Abs_modified2(delta_v, vel_v, checkSign(delta_v, vel_v), 1e9, 0.97, 0.6, init_vel_v)
    # Check if q[index] is greater than the threshold
    elseif delta_v > 0
        # For q[index] > threshold, you may adjust the parameters accordingly
        # Assuming the calculation method is symmetric or needs to be adjusted for this condition
        # This example uses the same function, but parameters may be adjusted based on your specific needs
        vel_v = qd[index]
        fz = -calculate_F_flore_Abs_modified2(delta_v, vel_v, checkSign(delta_v, vel_v), 1e9, 0.97, 0.6, init_vel_v)
    end
    
    # Conditional check to limit the maximum value of fz
    fz = fz > 1e20 ? 1e20 : fz
    
    # Print or return fz as needed
    # println("fz $fz")
    return fz
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

function calculate_F_D_Abs(q::Float64, v::Float64, flag_move_dir::Int, Eeq::Float64, n_damper::Float64)::Float64

    Fd = Eeq * (-q)^1.5 - v * Eeq *  n_damper

    return Fd
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
function calculate_F_flore_Abs_modified2(q::Float64, v::Float64, flag_move_dir::Int, Eeq::Float64, F2::Float64, E2::Float64, init_vel::Float64)::Float64

    # Conditional calculation of Fd based on the sign of q
        Fd = Eeq * (abs(q) / F2)^1.5
        #Fd *= (1.0 + 8.0 * (1.0 - E2) * v / (5.0 * E2 * (init_vel))) 


    return abs(Fd)
end

function calculate_g_fric_force0(ub, vel, damper, start_v, max_f)
    flag = 0
    if ub[2] > 0.0
        flag = 1
    elseif ub[2] < 0.0
        flag = -1
    else
    end
    vel_slide = project_vector_and_return_value(vel, ub)

    # if vel_slide<0.0005
    #     return  abs(smoothstep(abs(vel_slide), 0., start_v))*sign(vel[2])*ub*flag*0.
    # else
    return abs(smoothstep(abs(vel_slide), 0.001, start_v)) * sign(vel[2]) * ub * flag * max_f
    # end

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