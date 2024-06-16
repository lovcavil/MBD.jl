include("contact_geo.jl")

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
        fz = calculate_F_flore_Abs_modified((delta_v), (vel_v), checkSign(delta_v, vel_v), 1e9, 0.97, 0.6, init_vel_v)
    end
    debugDict=Dict("fz"=>fz)
    return fz>1e20 ? 1e20 : fz ,debugDict
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
function checkSign(a::Float64, b::Float64)::Int
    if a == 0 || b == 0
        return 0 # If either number is zero, return 0
    elseif a * b > 0
        return 1 # If the product is positive, they have the same sign
    else
        return -1 # If the product is negative, they have opposite signs
    end
end