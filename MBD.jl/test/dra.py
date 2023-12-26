import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the curve function
def curve(t):
    return np.cos(t), t#np.sin(t)

# Define the distance function to handle 2D grid
def distance(x, y, curve_func, t_vals):
    u, v = curve_func(t_vals)
    dist = np.sqrt((x[..., None] - u) ** 2 + (y[..., None] - v) ** 2)
    return np.min(dist, axis=2) 

# Create a meshgrid for x, y values
x = np.linspace(-10, 10, 400)
y = np.linspace(-5, 5, 400)
x, y = np.meshgrid(x, y)

# Calculate the distance for each point on the grid
t_vals = np.linspace(-2 * np.pi, 2 * np.pi, 400)
z = distance(x, y, curve, t_vals)

# Plotting the corrected surface
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Distance to Curve')

plt.title("Surface Representing Distance to the Curve")
plt.show()

# Define the curve function
def curve2(t):
    return 3*t+2, t#np.sin(t)

# Calculate the distance for each point on the grid
t_vals = np.linspace(-2 * np.pi, 2 * np.pi, 400)
z = distance(x, y, curve2, t_vals)

# Plotting the corrected surface
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Distance to Curve')

plt.title("Surface Representing Distance to the Curve")
plt.show()