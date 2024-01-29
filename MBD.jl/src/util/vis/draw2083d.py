
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import cos,sin,sqrt
data_jl = pd.read_csv(f"MBD.jl/src/util/208jl_solver.csv")
data_jl = pd.read_csv(f"MBD.jl/src/util/208jl_solver_II8.csv")
data_jl = pd.read_csv(f"MBD.jl/src/util/208jl_solver01281.csv")
# Extract 't' column as a numpy array
t = data_jl['t'].to_numpy()

# Create a 3xN numpy array from 'x1', 'y1', 'z1' columns
q = data_jl[['x1', 'y1', 'z1']].T.to_numpy()


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