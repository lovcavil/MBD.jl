function calculate_resultant_force(ub, delta, pen_vel, init_vel, k, f1, AF)
    q = abs(delta)
    v = abs(pen_vel)
    flag_dir = delta >= 0 ? 1 : -1
    init_vel = delta < 0 ? -abs(init_vel) : init_vel
    r = calculate_flore_force_magnitude(q, v, flag_dir, k, f1, AF, init_vel)
    return [-r * ub[1], -r * ub[2]]
end

function calculate_delta_spline(pos_part1, pos_part2, vel, jm, init_vel_id, clearance)
    x_diff = pos_part2[1] - pos_part1[1]
    y_diff = pos_part2[2] - pos_part1[2]
    cl = clearance
    nva = x_diff
    while true
        nvb = nva
        # Placeholder for c_cubspl function
        vva, vvb, vvc = 0, 0, 0  # Placeholder values
        ha = (x_diff - nvb + (y_diff - vva) * vvb)
        hb = (-1 - vvb * vvb + (y_diff - vva) * vvc)
        nva = nvb - ha / hb
        eer = nva - nvb
        if abs(eer) <= 1e-6
            break
        end
    end
    x_closest_point = nva
    y_closest_point = 0  # Placeholder value
    xb = pos_part2[1] - x_closest_point
    yb = pos_part2[2] - y_closest_point
    rr = sqrt(xb^2 + yb^2)
    uxb = xb / rr
    uyb = yb / rr
    dr = rr - cl
    r_delta = dr < 0 ? 0 : dr * sign(yb)
    vel2 = copy(vel)
    r_ub = [uxb, uyb]
    r_pen_vel = [0, 0]  # Placeholder value
    r_init_vel = [0]  # Placeholder value
    return r_delta, r_pen_vel, r_init_vel, r_ub
end

function compute_and_assign_result(p0, disp, vel, jm, init_vel_id, clearance, ub, cnf_type, k, f1, AF, b)
    delta, pen_vel, init_vel, r_ub = calculate_delta_spline(p0, disp, vel, jm, init_vel_id, clearance, ub)
    if cnf_type == 1
        result = cnf_hz(ub, delta, pen_vel, init_vel, k, f1, AF)
    elseif cnf_type == 2
        result = cnf_hz_d(ub, delta, pen_vel, init_vel, k, b)
    elseif cnf_type == 3
        result = cnf_flore(ub, delta, pen_vel, init_vel, k, f1, AF)
    end
    return [result[1], result[2], 0]
end

function calculate_flore_force_magnitude(q, v, flag_move_dir, Eeq, F2, E2, init_vel)
    if flag_move_dir == 0 || abs(init_vel) < 10.0
        init_vel = 0.0
    else
        v = 1.1 * abs(init_vel) * flag_move_dir if abs(v) > 1.1 * abs(init_vel) && abs(init_vel) > 10.0 else v * flag_move_dir
    end
    
    if init_vel == 0.0
        Fd = Eeq * (q / 1000.0 / F2) ^ 1.5
    else
        Fd = Eeq * (q / 1000.0 / F2) ^ 1.5 * (1.0 + 8 * (1.0 - E2) * v / (5.0 * E2 * init_vel))
    end
    return abs(Fd)
end
