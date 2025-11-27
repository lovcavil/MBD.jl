import numpy as np
import matplotlib.pyplot as plt
from numpy import cos, sin, sqrt

# 参数设置
def set_parameters():
    """设置系统参数"""
    params = {
        'm1': 30,  # 质量
        'g': 9.81,  # 重力加速度
        'dt': 0.001,  # 时间步长
        't_end': 5.0,  # 模拟结束时间
        'q0': [1, 1, 0],  # 初始位置
        'q_v0': [0, 0, 0]  # 初始速度
    }
    return params

# 初始化状态变量
def initialize_state(params):
    """初始化状态变量"""
    t = np.arange(0, params['t_end'] + params['dt'], params['dt'])
    n = len(t)
    q = np.zeros((3, n))  # 位置 (x, y, z)
    q_v = np.zeros((3, n))  # 速度 (vx, vy, vz)
    q_ac = np.zeros((3, n))  # 加速度 (ax, ay, az)
    
    # 设置初始条件
    q[:, 0] = params['q0']
    q_v[:, 0] = params['q_v0']
    
    return t, q, q_v, q_ac

# 构建矩阵 A 和向量 B
def build_matrices(q, q_v, params):
    """构建矩阵 A 和向量 B"""
    q1, q2, q3 = q
    qv1, qv2, qv3 = q_v
    m1 = params['m1']
    g = params['g']
    
    A = np.array([
        [m1, 0, 0, q1],
        [0, m1, 0, q2],
        [0, 0, m1, q3],
        [q1, q2, q3, 0]
    ])
    Qa = np.array([0, 0, -m1 * g])  # 外力
    garma = np.array([-qv1**2 - qv2**2 - qv3**2])  # 约束项
    B = np.concatenate((Qa, garma))
    
    return A, B

# 计算导数（速度 q_v 和加速度 q_ac）
def compute_derivatives(q, q_v, params):
    """计算速度 q_v 和加速度 q_ac"""
    A, B = build_matrices(q, q_v, params)
    temp = np.linalg.solve(A, B)
    q_ac = temp[:3]  # 加速度
    return q_v, q_ac

# 更新状态
def update_state(q, q_v, q_ac, i, params):
    """更新状态变量"""
    dt = params['dt']
    
    if i == 0:
        q_v[:, i + 1] = q_v[:, i] + dt * q_ac[:, i]
        q[:, i + 1] = q[:, i] + dt * q_v[:, i]
    else:
        q_v[:, i + 1] = q_v[:, i] + dt * (q_ac[:, i - 1] + q_ac[:, i]) / 2
        q[:, i + 1] = q[:, i] + dt * (q_v[:, i - 1] + q_v[:, i]) / 2

# 动力学求解
def solve_dynamics(params):
    """求解动力学方程"""
    t, q, q_v, q_ac = initialize_state(params)
    n = len(t)
    
    for i in range(n - 1):
        q_v[:, i], q_ac[:, i] = compute_derivatives(q[:, i], q_v[:, i], params)
        update_state(q, q_v, q_ac, i, params)
    
    return t, q, q_v, q_ac

# 绘图
def plot_results(t, q):
    """绘制结果"""
    # 位置随时间变化
    fig, axs = plt.subplots(3, 1, figsize=(6, 18))
    axs[0].plot(t, q[0, :])
    axs[0].set_title("Time vs q1")
    axs[0].set_xlabel("Time")
    axs[0].set_ylabel("q1")
    
    axs[1].plot(t, q[1, :])
    axs[1].set_title("Time vs q2")
    axs[1].set_xlabel("Time")
    axs[1].set_ylabel("q2")
    
    axs[2].plot(t, q[2, :])
    axs[2].set_title("Time vs q3")
    axs[2].set_xlabel("Time")
    axs[2].set_ylabel("q3")
    
    # 3D 轨迹图
    fig_3d = plt.figure(figsize=(8, 6))
    ax_3d = fig_3d.add_subplot(111, projection='3d')
    
    # 颜色映射
    colors = plt.cm.viridis(np.linspace(0, 1, len(t)))
    
    # 绘制轨迹
    for i in range(len(t) - 1):
        ax_3d.plot(q[0, i:i+2], q[1, i:i+2], q[2, i:i+2], color=colors[i])
    
    ax_3d.set_title("3D Trajectory with Color Map")
    ax_3d.set_xlabel("q1")
    ax_3d.set_ylabel("q2")
    ax_3d.set_zlabel("q3")
    ax_3d.set_xlim(-2, 2)
    ax_3d.set_ylim(-2, 2)
    ax_3d.set_zlim(-2, 2)
    
    plt.tight_layout()
    plt.show()

# 主函数
def main():
    """主函数"""
    params = set_parameters()
    t, q, q_v, q_ac = solve_dynamics(params)
    plot_results(t, q)

# 运行程序
if __name__ == "__main__":
    main()