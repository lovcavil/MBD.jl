import numpy as np
import matplotlib.pyplot as plt
from numpy import cos, sin, sqrt

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

# Baumgarte stabilization parameters
alpha = 50  # Position correction factor
beta = 50   # Velocity correction factor

for i in range(n - 1):
    q1 = q[0, i]
    q2 = q[1, i]
    q3 = q[2, i]
    qv1 = q_v[0, i]
    qv2 = q_v[1, i]
    qv3 = q_v[2, i]

    # Constraint function phi = (q1^2 + q2^2 + q3^2) / 2 = 0
    phi = (q1**2 + q2**2 + q3**2-1) / 2  # Constraint function

    # Time derivative of the constraint (phi_dot)
    dphi = q1 * qv1 + q2 * qv2 + q3 * qv3  # Dot of phi (constraint velocity)

    # Baumgarte stabilization term
    garma = -qv1**2 - qv2**2 - qv3**2  # Original damping term (adjusted)
    Baums = np.array([-qv1**2 - qv2**2 - qv3**2 - alpha * dphi - beta *beta* phi])  # Baumgarte stabilization terms

    # The system's generalized acceleration matrix
    A = np.array([[m1, 0, 0, q1],
                  [0, m1, 0, q2],
                  [0, 0, m1, q3],
                  [q1, q2, q3, 0]])

    # Calculate Qa (gravity term)
    Qa = np.array([0, 0, -m1 * g])  # Gravity force

    # Concatenate the gravity and Baumgarte stabilization forces
    B = np.concatenate((Qa, Baums))  # Total generalized forces

    # Solve for accelerations
    temp = np.linalg.solve(A, B)
    q_ac[0, i] = temp[0]
    q_ac[1, i] = temp[1]
    q_ac[2, i] = temp[2]
    q_ac[3, i] = temp[3]

    # Update velocities and positions using an explicit scheme (RK2)
    if i == 0:
        q_v[:, i+1] = q_v[:, i] + dt * q_ac[:, i]
        q[:, i+1] = q[:, i] + dt * q_v[:, i]
    else:
        q_v[:, i+1] = q_v[:, i] + dt * (q_ac[:, i-1] + q_ac[:, i]) / 2
        q[:, i+1] = q[:, i] + dt * (q_v[:, i-1] + q_v[:, i]) / 2

# Plotting
plt.figure()
plt.plot(t, q[0, :])
plt.title("Time vs q1")
plt.xlabel("Time")
plt.ylabel("q[0]")
plt.figure()
plt.plot(t, q[1, :])
plt.title("Time vs q2")
plt.xlabel("Time")
plt.ylabel("q[1]")
plt.figure()
plt.plot(t, q[2, :])
plt.title("Time vs q3")
plt.xlabel("Time")
plt.ylabel("q[2]")
plt.show()
