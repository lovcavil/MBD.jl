using DifferentialEquations
using Plots
using Dierckx


using LinearAlgebra

function calculate_delta_spline(pos_part, vel_glo, spline, clearance)
    println("calculate_delta_spline")
    x= pos_part[1]
    y= pos_part[2]
    println("x: ", x, " y: ", y)
    # println("vx: ", vel_glo[1], " vy: ", vel_glo[2])
    nva = x
    while true
        nvb = nva
        vva = spline(x)  
        vvb = derivative(spline,x,1)
        vvc = derivative(spline,x,2)
        ha = (x - nvb + (y - vva) * vvb)
        hb = (-1 - vvb * vvb + (y - vva) * vvc)
        nva = nvb - ha / hb
        eer = nva - nvb
        if abs(eer) <= 1e-6
            break
        end
    end
    x_closest_point = nva
    y_closest_point = spline(x_closest_point)
    # println("x_closest_point: ", x_closest_point, " y_closest_point: ", y_closest_point)
    xb = pos_part[1] - x_closest_point
    yb = pos_part[2] - y_closest_point
    # println("xb: ", xb, " yb: ", yb)
    rr = sqrt(xb^2 + yb^2)
    if rr>0
        uxb = xb / rr
        uyb = yb / rr
        dr = rr - clearance
    else
        uxb=0
        uyb=0
        dr=0
    end
    pen_pos_delta_local = -(dr < 0 ? 0 : dr )
    r_ub = [uxb, uyb]
    pen_vel_glo=project_vector_and_return_value(vel_glo, r_ub)
    r_init_vel = 0 
    return pen_pos_delta_local, pen_vel_glo, r_init_vel, r_ub
end

function calculate_resultant_force(ub, pen_pos_delta_local, pen_vel_glo, init_vel, k, f1, AF,damper)

    #*flag
    v_local=0.0
    iv=0.0
    if ub[2]>0
        v_local=-pen_vel_glo
        iv = abs(init_vel)
    elseif ub[2]<0
        v_local=pen_vel_glo
        iv = -abs(init_vel)
    else
        v_local=0.0
        iv = -1
    end
    # println("pen_pos_delta_local:$pen_pos_delta_local, pen_vel_glo:$pen_vel_glo, iv:$iv")
    println("*-pen_pos_delta_local:$pen_pos_delta_local, k:$k")
    # r = ContactMag(pen_pos_delta_local, v_local, k, f1, AF, iv)
    r = ContactMagBase(pen_pos_delta_local, k)
    if r>1e5
        r=1e5
    end
    fdMag=0.0
    fd=[0.,0.]
    # if r>0.0
    #     fdMag=(damper*pen_vel_glo)
    #     if ub[2]>0.0
    #         fd=[fdMag*ub[1],fdMag*ub[2]]
    #     elseif ub[2]<0.0
    #         fd=[-fdMag*ub[1],-fdMag*ub[2]]
    #     else
    #     end
    # end
    fd=[0.,0.]
    println("*------------------------------------------------ContactMag: ", r)
    return [-r * ub[1]+fd[1], -r * ub[2]+fd[2]]
end

function calculate_cont_damper_force(ub, vel,damper,start_v,max_f,pen_pos_delta_local,start_delta)
    flag=0
    if ub[2]>0.0
        flag=1
    elseif ub[2]<0.0
        flag=-1
    else
    end
    vel_slide=project_vector_and_return_value(vel, ub)
    
    return step_f(abs(pen_pos_delta_local),start_delta,1)*vel_slide*damper*sign(vel[2])*ub*flag
end

function calculate_load_carry_fric_force(ub, vel,damper,start_v,max_f)
    flag=0
    if ub[2]>0.0
        flag=1
    elseif ub[2]<0.0
        flag=-1
    else
    end
    vel_slide=project_vector_and_return_value(vel, ub)
    # println("vel:$vel, vel_slide:$vel_slide")
    aa = abs(step_f2(vel_slide, start_v, max_f)) * sign(vel[2]) * ub * flag
    println("calculate_load_carry_fric_force abs(step_f(vel_slide, start_v, max_f)):$(abs(step_f2(vel_slide, start_v, max_f))),sign(vel[2]):$(sign(vel[2])),damper:$aa")
    # println("ub * flag:$(ub * flag)")
    return abs(step_f2(vel_slide, start_v, max_f))*sign(vel[2])*ub*flag
    #return damper*sign(vel[2])*ub*flag
end

function step_f(v,start_v,max_f)
    if v > start_v
        return max_f
    elseif v < -start_v
        return -max_f
    else
        # Calculate proportion of v in the range (-start_v, start_v)
        # and scale it to the range (-max_f, max_f)
        return (v / start_v) * max_f
    end
end
function step_f2(v, start_v, max_f, zero_v=1e-1)
    if v > start_v
        return max_f
    elseif v < -start_v
        return -max_f
    elseif v > -zero_v && v < zero_v
        return 0
    else
        # Adjust calculation to ignore the zero_v zone and scale proportionally
        # First, adjust v to start the scale immediately outside the zero_v zone
        adjusted_v = v > 0 ? v - zero_v : v + zero_v
        # Then, adjust start_v to reduce its range by the zero_v zone on both sides
        adjusted_start_v = start_v - zero_v
        # Now, calculate the proportion of adjusted_v in the new range (-adjusted_start_v, adjusted_start_v)
        # and scale it to the range (-max_f, max_f)
        return (adjusted_v / adjusted_start_v) * max_f
    end
end
function ContactMag(delta_pos_local, v_local, Eeq, F2, E2, init_vel_local)
    # """delta_pos_local<0, v_local<,>,init_vel_localL<0,"""
    # println("-ContactMag-delta_pos_local:$delta_pos_local,v_local:$v_local,init_vel_local:$init_vel_local")
    # r=delta_pos_local < 0 ? ( Eeq * (abs(delta_pos_local)/F2)^1.5 * (1 + 8 .* (1-E2) .* v_local ./ 5 ./ E2 ./ init_vel_local ) ) : 0.0
    r=delta_pos_local < 0 ? ( Eeq * (abs(delta_pos_local))^1.5) : 0.0
    return r
end

function ContactMagBase(delta_pos_local, Eeq)
    # """delta_pos_local<0, v_local<,>,init_vel_localL<0,"""
    # println("-ContactMag-delta_pos_local:$delta_pos_local,v_local:$v_local,init_vel_local:$init_vel_local")
    # r=delta_pos_local < 0 ? ( Eeq * (abs(delta_pos_local)/F2)^1.5 * (1 + 8 .* (1-E2) .* v_local ./ 5 ./ E2 ./ init_vel_local ) ) : 0.0
    r=delta_pos_local < 0 ? ( Eeq * (abs(delta_pos_local))^1.5) : 0.0
    return r
end


function project_vector_and_return_value(p, u)
    dot_product = dot(p, u)
    mag_u = norm(u)
    if mag_u != 0.0
        return dot_product / mag_u
    else
        return 0.0
    end
end


function calculate_contact_geo(d_contact,q,qd)
    f1,AF,B5,Eeq,damper,start_v,max_f,start_delta,faaa = d_contact["p"]
    spl=d_contact["guide"]
    b=d_contact["b"]

    index_x = 7 * (b - 1) + 1
    index_y = 7 * (b - 1) + 2
    pos_part=[q[index_x],q[index_y]]
    vel=[qd[index_x],qd[index_y]]
    init_vel = B5
    clearance = 0.25
    pen_pos_delta_local, pen_vel_glo, r_init_vel, ub = calculate_delta_spline(pos_part, vel, spl, clearance)

    Fx , Fy = calculate_resultant_force(ub, pen_pos_delta_local, pen_vel_glo, init_vel, Eeq, f1, AF,damper)
    Fdx , Fdy=0.,0.
    if Fx!=0.0
        Fdx , Fdy = calculate_cont_damper_force(ub, vel, damper, start_v,0,pen_pos_delta_local,start_delta)
    end
    Ffx , Ffy =  faaa*calculate_load_carry_fric_force(ub, vel,damper,start_v,max_f)
    # Ffx=0.
    Fdx , Fdy=0.,0.
    # Ffx , Ffy=0.,0.
    println("Fx: ", Fx, "Fy: ", Fy)
    println("Fdx: ", -Fdx, "Fdy: ", -Fdy)
    println("Ffx: ", -Ffx, "Ffy: ", -Ffy)
    

    # return Fx , Fy
    return Fx-Fdx-Ffx , Fy-Fdy-Ffy
end



# # 构造数据点
# x = [0.0, 100.0, 200.0, 300.0, 400.0,600,800,1000]
# y = [1.0, 1.0, 1.0, 2.0, 1.0,1,1,1]

# # 使用 Dierckx.jl 中的 `Spline1D` 函数创建样条曲线
# spl = Spline1D(x, y)

# ub=[1.0,0.0]
# f1=0.97
# AF=0.6
# M = 20.0 # 示例值
# g = 9.81 # 重力加速度
# B5=-1.1#*testdir
# damper=0.005
# p=Any[ub,f1,AF,M,g,B5,spl,damper]
# u0=[200.0,1.0,1000.0,10.0] # 初始穿透深度和初始速度，需要根据实际情况调整

# # Time span
# tspan = (0.0, 0.8) # From 0 to 10 seconds

# # Define the problem
# prob = ODEProblem(myfunhertz!, u0, tspan, p)

# # Solve the problem
# sol = solve(prob,Tsit5(), reltol=1e-6, abstol=1e-8, maxiters=1000,dtmin=1e-10, dtmax=1e-2)

