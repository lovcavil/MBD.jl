import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Set the font to Times New Roman and adjust font sizes
matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['font.size'] = 10  # Base font size

# Load data from the file, skipping the header
data = np.genfromtxt(r'C:\Users\lovca\OneDrive\Articles\10.Working\MBD.jl\MBD.jl\src\problem\GUIDE_M.csv', delimiter=',', skip_header=1)
x = data[:, 0]
y = data[:, 1]
# z = data[:, 2]  # Z is ignored as per your instruction

# Plot XY curve in Figure 1 with equal scaling
plt.figure(1, figsize=(5, 2))
plt.plot(x, y, linestyle='-', color='black', linewidth=1)
plt.xlabel('X', fontsize=12)
plt.ylabel('Y', fontsize=12)
plt.title('XY Curve with Equal Scaling', fontsize=14)
plt.axis('equal')  # Set equal scaling on both axes
plt.gca().invert_xaxis()
plt.grid(True)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

# Save Figure 1 to the script folder with DPI 600
plt.savefig('Figure1_XYCurve.png', dpi=600, bbox_inches='tight')

# Compute first derivatives
dx = np.gradient(x)
dy = np.gradient(y)

# Compute second derivatives
ddx = np.gradient(dx)
ddy = np.gradient(dy)

# Compute curvature
curvature = np.abs(dx * ddy - dy * ddx) / (dx**2 + dy**2)**1.5

# Plot curvature vs X in Figure 2 with X-axis reversed
plt.figure(2, figsize=(5, 2))
plt.plot(x, curvature, linestyle='-', color='black', linewidth=1)
plt.xlabel('X', fontsize=12)
plt.ylabel('Curvature', fontsize=12)
plt.title('Curvature vs X (X-axis Reversed)', fontsize=14)
plt.gca().invert_xaxis()  # Reverse the X-axis
plt.grid(True)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

# Save Figure 2 to the script folder with DPI 600
plt.savefig('Figure2_Curvature.png', dpi=600, bbox_inches='tight')

# Display the plots
plt.show()
