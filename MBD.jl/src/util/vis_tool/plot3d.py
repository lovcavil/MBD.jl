import json
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import os

script_dir = os.path.dirname(__file__)
# Load CSV data
csv_file = os.path.join(script_dir, 'adams_output.csv')  # Update this with your CSV file path
data = pd.read_csv(csv_file)

# Load column groups from JSON
json_file = os.path.join(script_dir, 'column_groups.json')  # Update this with your JSON file path
with open(json_file, 'r') as file:
    column_groups = json.load(file)

# Create 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Get max value of 't' for color mapping
max_t = data['t'].max()

# Plot each group of columns as a 3D curve
for curve_name, columns in column_groups.items():
    x, y, z = columns
    t = data['t'] / max_t  # Normalize 't' for color mapping

    # Plot each segment of the curve with a color corresponding to the 't' value
    for i in range(len(data) - 1):
        ax.plot(data[x][i:i+2], data[y][i:i+2], data[z][i:i+2], color=cm.viridis(t.iloc[i]))

# Set labels and title
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.set_title('3D Curves Plot')

# Show legend
ax.legend()

# Show plot
plt.show()
