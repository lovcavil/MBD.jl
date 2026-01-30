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
dt = 0.0001
t = np.arange(0, 10, dt)
n = len(t)
q = np.zeros((9, n))
q_v = np.zeros((9, n))
q_ac = np.zeros((9, n))
q[:, 0] = [0, 1, 0, 0, 0, 2*np.pi, 270, 0, 0]
q_v[:, 0] = [0, 0, 0, 0, 0, 0, 0, 0, 0]
force = np.zeros((9, n))

for i in range(n):
    q1 = q[0, i]
    q2 = q[1, i]
    q3 = q[2, i]
    qv1 = q_v[0, i]
    qv2 = q_v[1, i]
    qv3 = q_v[2, i]

    #sph
    A=np.array([[m1,0,0,q1],
                [0,m1,0,q2],
                [0,0,m1,q3],
                [q1,q2,q3,0]])

    # Calculate Qa
    Qa = np.array([0,0,-m1*g])
    #sph
    garma = np.array(
        [-qv1*qv1-qv2*qv2-qv3*qv3])

    
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

# Plotting
plt.figure()
plt.plot(t, q[0, :])
plt.title("d Time vs q1")
plt.xlabel("Time")
plt.ylabel("q_ac[7]")
plt.figure()
plt.plot(t, q[1, :])
plt.title("d Time vs q2")
plt.xlabel("Time")
plt.ylabel("q_ac[7]")
plt.figure()
plt.plot(t, q[2, :])
plt.title("d Time vs q3")
plt.xlabel("Time")
plt.ylabel("q_ac[7]")
# plt.figure()
# plt.plot(t, q[2, :])
# plt.title("Time vs l1")
# plt.xlabel("Time")
# plt.ylabel("q_ac[7]")
# plt.figure()
# plt.plot(t, q[3, :])
# plt.title("Time vs l2")
# plt.xlabel("Time")
# plt.ylabel("q_ac[7]")
plt.show()
