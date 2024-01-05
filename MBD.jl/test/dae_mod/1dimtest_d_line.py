import numpy as np
import matplotlib.pyplot as plt
from numpy import cos,sin,sqrt
# Given parameters
a = 1000
b = 1000
m1 = 30
m2 = 0.19
m3 = 0.70
I1 = 7360
I2 = 880
I3 = 360
g = 9.81
l1 = 50
l2 = 220
w = 5 * np.pi
dt = 0.001
t = np.arange(0, 5.1, dt)
n = len(t)
q = np.zeros((9, n))
q_v = np.zeros((9, n))
q_ac = np.zeros((9, n))
q[:, 0] = [1, 1, 0, 0, 0, 2*np.pi, 270, 0, 0]
q_v[:, 0] = [0, 0, 0, 0, 0, 0, 0, 0, 0]
force = np.zeros((9, n))

for i in range(n):
    q1 = q[0, i]
    q2 = q[1, i]
    q3 = q[2, i]
    qv1 = q_v[0, i]
    qv2 = q_v[1, i]
    qv3 = q_v[2, i]

    # Building matrix A
    # A = np.array([
    #     [m1,  qv1],
    #     [ qv1, 0]
    # ])
    #sph
    A=np.array([[m1,0,0,q1],
                [0,m1,0,q2],
                [0,0,m1,q3],
                [q1,q2,q3,0]])
    #sph+plain 2x-y-1=0
    A=np.array([[m1,0,0,q1,2],
                [0,m1,0,q2,-1],
                [0,0,m1,q3,0],
                [q1,q2,q3,0,0],
                [2,-1,0,0,0]])

    # Calculate Qa
    Qa = np.array([0,0,-m1*g])
    #sph+plain x=0
    garma = np.array(
        [-qv1*qv1-qv2*qv2-qv3*qv3,0])
    #sph+plain y^2/2=0
    # garma = np.array(
    #     [-qv1*qv1-qv2*qv2-qv3*qv3,-qv1*qv1])
    #garma = np.array(
    #    [0])
    #garma=np.array([0])
    
    B = np.concatenate((Qa, garma))
    #B=Qa
    # Solve the system
    #print(A)
    #print(B)
    #print('\n')
    temp = np.linalg.solve(A, B)
    q_ac[0, i] = temp[0]
    q_ac[1, i] = temp[1]
    q_ac[2, i] = temp[2]
    q_ac[3, i] = temp[3]
    q_ac[4, i] = temp[4]
    #q_ac[7, i] = 0
    #q_ac[8, i] = 0
    #force[:, i] = temp[9:]
    
    if i == n-1:
        break
    if i == 0:
        q_v[:, i+1] = q_v[:, i] + dt*q_ac[:, i]
        q[:, i+1] = q[:, i] + dt*q_v[:, i]
    else:
        q_v[:, i+1] = q_v[:, i] + dt*(q_ac[:, i-1] + q_ac[:, i])/2
        q[:, i+1] = q[:, i] + dt*(q_v[:, i-1] + q_v[:, i])/2
        
    #q_v[7, i+1] = 0
    #q_v[8, i+1] = 0
    #q[7, i+1] = 0
    #q[8, i+1] = 0


# Creating subplots
fig, axs = plt.subplots(3, 1, figsize=(6, 18))

# First subplot
axs[0].plot(t, q[0, :])
axs[0].set_title("Time vs q1")
axs[0].set_xlabel("Time")
axs[0].set_ylabel("q_ac[7]")

# Second subplot
axs[1].plot(t, q[1, :])
axs[1].set_title("Time vs q2")
axs[1].set_xlabel("Time")
axs[1].set_ylabel("q_ac[7]")

# Third subplot
axs[2].plot(t, q[2, :])
axs[2].set_title("Time vs q3")
axs[2].set_xlabel("Time")
axs[2].set_ylabel("q_ac[7]")

# Additional 3D plot

fig_3d_colored = plt.figure(figsize=(8, 6))
ax_3d_colored = fig_3d_colored.add_subplot(111, projection='3d')

# Creating a color map based on time
colors = plt.cm.viridis(np.linspace(0, 1, len(t)))

# Plotting each segment of the line with a color corresponding to its time point
for i in range(len(t)-1):
    ax_3d_colored.plot3D(q[0, i:i+2], q[1, i:i+2], q[2, i:i+2], color=colors[i])

ax_3d_colored.set_title("3D Plot with Color Map")
ax_3d_colored.set_xlabel("q1")
ax_3d_colored.set_ylabel("q2")
ax_3d_colored.set_zlabel("q3")
ax_3d_colored.set_xlim(-2, 2)
ax_3d_colored.set_ylim(-2, 2)
ax_3d_colored.set_zlim(-2, 2)
plt.tight_layout()
plt.show()