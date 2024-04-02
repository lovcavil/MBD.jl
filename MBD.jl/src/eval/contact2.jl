using DifferentialEquations
using Plots
using Dierckx
include("contact1.jl")

using LinearAlgebra

function calculate_delta_spline(pos_part, vel_glo, spline, clearance)
    println("calculate_delta_spline")
    x= pos_part[1]
    y= pos_part[2]
    println("x: ", x, " y: ", y)
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
    println("x_closest_point: ", x_closest_point, " y_closest_point: ", y_closest_point)
    xb = pos_part[1] - x_closest_point
    yb = pos_part[2] - y_closest_point
    println("xb: ", xb, " yb: ", yb)
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
    println("pen_pos_delta_local:$pen_pos_delta_local, pen_vel_glo:$pen_vel_glo, iv:$iv")
    println("pen_pos_delta_local:$pen_pos_delta_local, v_local:$v_local, iv:$iv")
    r = ContactMag(pen_pos_delta_local, v_local, k, f1, AF, iv)
    fdMag=0.0
    fd=[0.,0.]
    if r>0.0
        fdMag=k*damper*pen_vel_glo
        if ub[2]>0.0
            fd=[-fdMag*ub[1],-fdMag*ub[2]]
        elseif ub[2]<0.0
            fd=[-fdMag*ub[1],-fdMag*ub[2]]
        else
        end
    end
    # fd=[0.,0.]
    println("ContactMag: ", r)
    println("fd:$fd ")
    return [-r * ub[1]+fd[1], -r * ub[2]+fd[2]]
end

function ContactMag(delta_pos_local, v_local, Eeq, F2, E2, init_vel_local)
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

# function compute_and_assign_result(p0, disp, vel, jm, init_vel_id, clearance, ub, cnf_type, k, f1, AF, b)
#     delta, pen_vel, init_vel, r_ub = calculate_delta_spline(p0, disp, vel, jm, init_vel_id, clearance, ub)
#     if cnf_type == 1
#         result = cnf_hz(ub, delta, pen_vel, init_vel, k, f1, AF)
#     elseif cnf_type == 2
#         result = cnf_hz_d(ub, delta, pen_vel, init_vel, k, b)
#     elseif cnf_type == 3
#         result = cnf_flore(ub, delta, pen_vel, init_vel, k, f1, AF)
#     end
#     return [result[1], result[2], 0]
# end



function myfunhertz!(du, u, p, t)
    ub,f1,AF,M,g,B5,spl,damper = p
    q1,q2, v1,v2 = u
    pos_part=[q1,q2]
    vel=[v1,v2]
    Eeq = 1e5
    init_vel = B5
    clearance = 0.05
    pen_pos_delta_local, pen_vel_glo, r_init_vel, ub = calculate_delta_spline(pos_part, vel, spl, clearance)
    println("t: ", t)
    println("ub: $ub, pen_vel_glo:$pen_vel_glo")
    Fx , Fy =  calculate_resultant_force(ub, pen_pos_delta_local, pen_vel_glo, init_vel, Eeq, f1, AF,damper)
    println("Fx: ", Fx, " Fy: ", Fy)
    println("vx: ", v1, " vy: ", v2)
    du[1] = v1 # 穿透深度的导数是速度
    du[3] = Fx / M  # 加速度
    du[2] = v2 # 穿透深度的导数是速度
    du[4] = Fy / M - g # 加速度
end



# testdir=1
#testdir=-1

# 构造数据点
x = [0.0, 100.0, 200.0, 300.0, 400.0,600,800,1000]
y = [1.0, 1.0, 1.0, 2.0, 1.0,1,1,1]

# 使用 Dierckx.jl 中的 `Spline1D` 函数创建样条曲线
spl = Spline1D(x, y)

ub=[1.0,0.0]
f1=0.97
AF=0.6
M = 20.0 # 示例值
g = 9.81 # 重力加速度
B5=-1.1#*testdir
damper=0.005
p=Any[ub,f1,AF,M,g,B5,spl,damper]
u0=[200.0,1.0,1000.0,10.0] # 初始穿透深度和初始速度，需要根据实际情况调整

# Initial conditions
#u0 = [0.0, -1.0*testdir] # Initial position and velocity

# Time span
tspan = (0.0, 0.8) # From 0 to 10 seconds

# Define the problem
prob = ODEProblem(myfunhertz!, u0, tspan, p)

# Solve the problem
sol = solve(prob,Tsit5(), reltol=1e-6, abstol=1e-8, maxiters=1000,dtmin=1e-10, dtmax=1e-2)

# Plot the results
p0=plot(legend=:topright)
plot!(sol[1,:],sol[2,:], xlabel="Position", ylabel="Position", label="Position", lw=2)
display(p0)
p1=plot(legend=:topright)
plot!(sol.t,sol[1,:], xlabel="Time", ylabel="Position", label="Position", lw=2)
display(p1)
p2=plot(legend=:topright)
plot!(sol.t,sol[2,:], xlabel="Time", ylabel="Position", label="Position", lw=2)
display(p2)
p3=plot(legend=:topright)
plot!(sol.t,sol[3,:], xlabel="Time", ylabel="V", label="V", lw=2)
display(p3)
p4=plot(legend=:topright)
plot!(sol.t,sol[4,:], xlabel="Time", ylabel="V", label="V", lw=2)
display(p4)
#plot(legend=:topright)
#plot!(sol, vars=(1,2), xlabel="Time", ylabel="Velocity", label="Velocity", lw=2)

# Generate a finer set of x values for plotting the spline smoothly
x_fine = range(100, 1000, length=300)

# Evaluate the spline at the finer set of x values
y_fine = [spl(x) for x in x_fine]

# Plot the original data points
# a=plot(x, y, seriestype=:scatter, label="Data Points")

# Plot the spline curve
a=plot(x_fine, y_fine, label="Spline Interpolation")
plot!(sol[1,:],sol[2,:], xlabel="Position", ylabel="Position", label="Position", lw=2)
# Display the plot
display(a)
