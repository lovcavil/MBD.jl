import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd

# Define the path to your Excel file
file_path = 'E:/New Folder/A16/2024-10-25-16-24-09/nihao.xlsx'

# Read the Excel file into a DataFrame
df = pd.read_excel(file_path)

# Filter the DataFrame
df = df[df["abstol"] > 0.0001].reset_index(drop=True)  # Reset index after filtering
df = df[np.abs(df["ΔE"]) < 100000].reset_index(drop=True)
# Define sample data (Ensure these match in length)
n = num_rows = df.shape[0]  # Total number of points
x = df["ΔE"].apply(lambda value: np.log10(np.abs(value)))  # x-coordinates with log10(abs) transformation
y = df["Avg Time (s)"]  # y-coordinates

# Provided name labels
name_labels = df["Solver"]

# Sort unique names alphabetically and create a colormap
unique_names = sorted(set(name_labels))  # Sort unique labels alphabetically
num_groups = len(unique_names)
cmap = cm.get_cmap('viridis', num_groups)  # Choose a colormap
cmap = cm.get_cmap('turbo', num_groups)  # Choose a colormap
# Assign a color to each unique name based on alphabetical order
color_map = {name: cmap(i) for i, name in enumerate(unique_names)}

# Plot the scatter plot with colors based on alphabetical order
plt.figure(figsize=(10, 8))
for name in unique_names:
    group_idx = np.array(name_labels) == name  # Logical index for the current group
    plt.scatter(x[group_idx], y[group_idx], color=color_map[name], label=name, s=50)

# Add labels to each point
for i in range(n):
    plt.text(x.iloc[i], y.iloc[i], name_labels.iloc[i], fontsize=8, ha='right', va='bottom')

# Title, axis labels, and legend
plt.title('Scatter Plot with Alphabetical Color Mapping')
plt.xlabel('Log10(ΔE)')
plt.ylabel('Avg Time (s)')
plt.legend(title="Groups", bbox_to_anchor=(1.0, 1), loc='upper left')  # Show legend with group names

plt.show()
