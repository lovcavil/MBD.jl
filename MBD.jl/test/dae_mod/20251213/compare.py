
import numpy as np
import matplotlib.pyplot as plt

# --- Data from pend.py ---
te=100
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
t = np.arange(0, te, dt)
n = len(t)
q_pend = np.zeros((9, n))
q_v_pend = np.zeros((9, n))
q_ac_pend = np.zeros((9, n))
q_pend[:, 0] = [0, 1, 0, 0, 0, 2*np.pi, 270, 0, 0]
q_v_pend[:, 0] = [0, 0, 0, 0, 0, 0, 0, 0, 0]

for i in range(n-1):
    q1 = q_pend[0, i]
    q2 = q_pend[1, i]
    q3 = q_pend[2, i]
    qv1 = q_v_pend[0, i]
    qv2 = q_v_pend[1, i]
    qv3 = q_v_pend[2, i]

    A=np.array([[m1,0,0,q1],
                [0,m1,0,q2],
                [0,0,m1,q3],
                [q1,q2,q3,0]])
    Qa = np.array([0,0,-m1*g])
    garma = np.array([-qv1*qv1-qv2*qv2-qv3*qv3])
    B = np.concatenate((Qa, garma))
    temp = np.linalg.solve(A, B)
    q_ac_pend[0, i] = temp[0]
    q_ac_pend[1, i] = temp[1]
    q_ac_pend[2, i] = temp[2]
    q_ac_pend[3, i] = temp[3]
    
    if i == 0:
        q_v_pend[:, i+1] = q_v_pend[:, i] + dt*q_ac_pend[:, i]
        q_pend[:, i+1] = q_pend[:, i] + dt*q_v_pend[:, i]
    else:
        q_v_pend[:, i+1] = q_v_pend[:, i] + dt*(q_ac_pend[:, i-1] + q_ac_pend[:, i])/2
        q_pend[:, i+1] = q_pend[:, i] + dt*(q_v_pend[:, i-1] + q_v_pend[:, i])/2

# --- Data from pend_BS.py ---
dt = 0.0001
t_bs = np.arange(0, te, dt)
n_bs = len(t_bs)
q_bs = np.zeros((9, n_bs))
q_v_bs = np.zeros((9, n_bs))
q_ac_bs = np.zeros((9, n_bs))
q_bs[:, 0] = [0, 1, 0, 0, 0, 2*np.pi, 270, 0, 0]
q_v_bs[:, 0] = [0, 0, 0, 0, 0, 0, 0, 0, 0]

alpha = 400
beta = 200

for i in range(n_bs - 1):
    q1 = q_bs[0, i]
    q2 = q_bs[1, i]
    q3 = q_bs[2, i]
    qv1 = q_v_bs[0, i]
    qv2 = q_v_bs[1, i]
    qv3 = q_v_bs[2, i]

    phi = (q1**2 + q2**2 + q3**2 -1) / 2
    dphi = q1 * qv1 + q2 * qv2 + q3 * qv3
    garma = -qv1**2 - qv2**2 - qv3**2
    Baums = np.array([garma-2*alpha * dphi - beta*beta * phi])
    A = np.array([[m1, 0, 0, q1],
                  [0, m1, 0, q2],
                  [0, 0, m1, q3],
                  [q1, q2, q3, 0]])
    Qa = np.array([0, 0, -m1 * g])
    B = np.concatenate((Qa, Baums))
    temp = np.linalg.solve(A, B)
    q_ac_bs[0, i] = temp[0]
    q_ac_bs[1, i] = temp[1]
    q_ac_bs[2, i] = temp[2]
    q_ac_bs[3, i] = temp[3]

    if i == 0:
        q_v_bs[:, i+1] = q_v_bs[:, i] + dt * q_ac_bs[:, i]
        q_bs[:, i+1] = q_bs[:, i] + dt * q_v_bs[:, i]
    else:
        q_v_bs[:, i+1] = q_v_bs[:, i] + dt * (q_ac_bs[:, i-1] + q_ac_bs[:, i]) / 2
        q_bs[:, i+1] = q_bs[:, i] + dt * (q_v_bs[:, i-1] + q_v_bs[:, i]) / 2

# --- Plotting ---
plt.figure(figsize=(12, 12))

# precompute constraint and speed energy terms
pos_sum_sq_pend = q_pend[0, :]**2 + q_pend[1, :]**2 + q_pend[2, :]**2
pos_sum_sq_bs = q_bs[0, :]**2 + q_bs[1, :]**2 + q_bs[2, :]**2
vel_sum_sq_pend = q_v_pend[0, :]**2 + q_v_pend[1, :]**2 + q_v_pend[2, :]**2
vel_sum_sq_bs = q_v_bs[0, :]**2 + q_v_bs[1, :]**2 + q_v_bs[2, :]**2

plt.subplot(5, 1, 1)
plt.plot(t, q_pend[0, :], label='pend.py')
plt.plot(t_bs, q_bs[0, :], label='pend_BS.py', linestyle='--')
plt.title("Comparison of q1")
plt.xlabel("Time")
plt.ylabel("q1")
plt.legend()

plt.subplot(5, 1, 2)
plt.plot(t, q_pend[1, :], label='pend.py')
plt.plot(t_bs, q_bs[1, :], label='pend_BS.py', linestyle='--')
plt.title("Comparison of q2")
plt.xlabel("Time")
plt.ylabel("q2")
plt.legend()

plt.subplot(5, 1, 3)
plt.plot(t, q_pend[2, :], label='pend.py')
plt.plot(t_bs, q_bs[2, :], label='pend_BS.py', linestyle='--')
plt.title("Comparison of q3")
plt.xlabel("Time")
plt.ylabel("q3")
plt.legend()

plt.subplot(5, 1, 4)
plt.plot(t, pos_sum_sq_pend, label='pend.py')
plt.plot(t_bs, pos_sum_sq_bs, label='pend_BS.py', linestyle='--')
plt.title("Comparison of q1^2 + q2^2 + q3^2")
plt.xlabel("Time")
plt.ylabel("Sum of squares (pos)")
plt.legend()

plt.subplot(5, 1, 5)
plt.plot(t, vel_sum_sq_pend, label='pend.py')
plt.plot(t_bs, vel_sum_sq_bs, label='pend_BS.py', linestyle='--')
plt.title("Comparison of qv1^2 + qv2^2 + qv3^2")
plt.xlabel("Time")
plt.ylabel("Sum of squares (vel)")
plt.legend()

plt.tight_layout()
plt.show()
