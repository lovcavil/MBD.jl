import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd

# Define the path to your Excel file
file_path = 'E:/New Folder/A16/con4.xlsx'

# Read the Excel file into a DataFrame
df = pd.read_excel(file_path)

# Define sample data
df = df[df["ΔE"] != 0.0].reset_index(drop=True)
df = df[np.abs(df["Avg Time (s)"]) < 115].reset_index(drop=True)
n = num_rows = df.shape[0]  # Total number of points
#y = df["ΔE"].apply(lambda value: np.log10(np.abs(value)))  # y-coordinates with log10(abs) transformation
y = df["ΔE"].apply(lambda value: (np.abs(value)))
x = df["Avg Time (s)"]  # x-coordinates

# # Clip values to ensure boundary points are included
y_clipped = np.clip(y, None, 5)  # Clip y at 5 (10^5)
x_clipped = np.clip(x, None, 119)  # Clip x at 119

# Provided name labels
name_labels = df["Solver"]

# Sort unique names alphabetically and create a colormap
unique_names = sorted(set(name_labels))  # Sort unique labels alphabetically
num_groups = len(unique_names)
cmap = cm.get_cmap('turbo', num_groups)  # Choose a colormap

# Assign a color to each unique name based on alphabetical order
color_map = {name: cmap(i) for i, name in enumerate(unique_names)}

# Plot the scatter plot with colors based on alphabetical order
plt.figure(figsize=(10, 8))
for name in unique_names:
    group_idx = np.array(name_labels) == name  # Logical index for the current group
    plt.scatter(y_clipped[group_idx], y_clipped[group_idx], color=color_map[name], label=name, s=50)

# Add labels and arrows for points at the boundary
for i in range(n):
    xi, yi = y_clipped.iloc[i], y_clipped.iloc[i]
    plt.text(xi, yi, name_labels.iloc[i], fontsize=8, ha='right', va='bottom')
    # Add an arrow indicator if the point was clipped
    if y.iloc[i] > 5:
        plt.annotate('', xy=(xi, 5), xytext=(xi, 5.1), arrowprops=dict(arrowstyle="->", color='red'))
    if x.iloc[i] > 119:
        plt.annotate('', xy=(119, yi), xytext=(119.5, yi), arrowprops=dict(arrowstyle="->", color='red'))

# Title, axis labels, and legend
plt.title('Scatter Plot with Alphabetical Color Mapping and Boundary Indicators (Swapped Axes)')
plt.xlabel('Avg Time (s)')
plt.ylabel('Log10(ΔE)')
plt.legend(title="Groups", bbox_to_anchor=(1.0, 1), loc='upper left')  # Show legend with group names

plt.show()
