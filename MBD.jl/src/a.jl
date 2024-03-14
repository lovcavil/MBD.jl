using DifferentialEquations, Plots
include("./eval/QAC.jl")


const F2 = 0.97 # 示例值
	  E2 = 0.6   
const M = 20.0 # 示例值
const g = 9.81 # 重力加速度
#	  V = 0.4
B5=-0.01
const E=1e9

function myfunhertz!(du, u, p, t)
    q, v = u
    init_vel=B5
    flag_move_dir=checkSign(q,v)
    if q<0
        Fd=calculate_F_flore_Abs_modified(q, v,  checkSign(q,v), E, F2, E2, init_vel)
    else
        Fd=0.0
    end
    du[1] = v # 穿透深度的导数是速度
    du[2] = Fd / M - g # 加速度
end

u0 = [0.0, B5] # 初始穿透深度和初始速度，需要根据实际情况调整
tspan = (0.0,5.0)

prob = ODEProblem(myfunhertz!, u0, tspan)
sol = solve(prob, Tsit5(), reltol=1e-15, abstol=1e-15)
q_vals = sol[1, :] # 穿透深度
v= sol[2, :]
	

q_vals = sol[1, :] # Penetration depth
v = sol[2, :]      # Velocity
init_vel = B5
Fd_vals = map((q, v) -> q < 0 ? calculate_F_flore_Abs_modified(q, v,  checkSign(q,v), E, F2, E2, init_vel) : 0.0, q_vals, v)

# 绘图
p = plot()
plot!(q_vals, Fd_vals, label="Contact Force vs. Penetration Depth", xlabel="Penetration Depth (m)", ylabel="Contact Force (N)", title="Hertz Contact Model")
display(p)
p = plot()
plot!(sol.t, q_vals)
display(p)
p = plot()
plot!(sol.t, v)
display(p)
