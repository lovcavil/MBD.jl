import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Set the font to Times New Roman and adjust font sizes
matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['font.size'] = 10

# Load data from the file, skipping the header
data = np.genfromtxt(
    r'C:\Users\lovca\OneDrive\Articles\10.Working\MBD.jl\MBD.jl\src\problem\GUIDE_M_A26MID_Report.csv',
    delimiter=',',
    skip_header=1
)

x = data[:, 0]
y = data[:, 1]

# =========================
# Figure 1: XY curve
# =========================
plt.figure(1, figsize=(5, 2))
plt.plot(x, y, linestyle='-', color='black', linewidth=1)
plt.xlabel('X', fontsize=12)
plt.ylabel('Y', fontsize=12)
plt.title('XY Curve with Equal Scaling', fontsize=14)
plt.axis('equal')
plt.gca().invert_xaxis()
plt.grid(True)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.savefig('Figure1_XYCurveA28.png', dpi=600, bbox_inches='tight')

# =========================
# Curvature calculation using arc length
# =========================
# Arc length increment
ds = np.sqrt(np.diff(x)**2 + np.diff(y)**2)

# Construct cumulative arc length
s = np.zeros(len(x))
s[1:] = np.cumsum(ds)

# First derivatives with respect to arc length s
dx_ds = np.gradient(x, s)
dy_ds = np.gradient(y, s)

# Second derivatives with respect to arc length s
d2x_ds2 = np.gradient(dx_ds, s)
d2y_ds2 = np.gradient(dy_ds, s)

# Curvature formula
curvature = np.abs(dx_ds * d2y_ds2 - dy_ds * d2x_ds2) / (dx_ds**2 + dy_ds**2)**1.5

# Replace NaN or inf if any
curvature = np.nan_to_num(curvature, nan=0.0, posinf=0.0, neginf=0.0)

# =========================
# Exclude zero-curvature / near-zero-curvature points
# =========================
tol = 1e-6  # threshold, can adjust
valid_mask = curvature > 0.0025
valid_curvature = curvature[valid_mask]

if len(valid_curvature) > 0:
    avg_curvature = np.mean(valid_curvature)
    max_curvature = np.max(valid_curvature)
    min_curvature = np.min(valid_curvature)
else:
    avg_curvature = 0.0
    max_curvature = 0.0
    min_curvature = 0.0

print(f'Number of total points      : {len(curvature)}')
print(f'Number of valid curve points: {len(valid_curvature)}')
print(f'Average curvature           : {avg_curvature:.8f}')
print(f'Min curvature (valid part)  : {min_curvature:.8f}')
print(f'Max curvature (valid part)  : {max_curvature:.8f}')

# =========================
# Figure 2: Curvature vs X
# =========================
plt.figure(2, figsize=(5, 2))
plt.plot(x, curvature, linestyle='-', color='black', linewidth=1, label='Curvature')
plt.plot(x[valid_mask], curvature[valid_mask], 'r.', markersize=2, label='Valid curvature')
plt.xlabel('X', fontsize=12)
plt.ylabel('Curvature', fontsize=12)
plt.title('Curvature vs X (X-axis Reversed)', fontsize=14)
plt.gca().invert_xaxis()
plt.grid(True)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(fontsize=8)
plt.savefig('Figure2_CurvatureA28.png', dpi=600, bbox_inches='tight')

plt.show()