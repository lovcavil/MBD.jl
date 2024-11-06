import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm

# Define the paths to your Excel files
file_path1 = r'E:/New Folder/A16/tol2.xlsx'
file_path2 = r'E:/New Folder/A16/tol3.xlsx'
file_path3 = r'E:/New Folder/A16/tol4.xlsx'

# Load the DataFrames
df1 = pd.read_excel(file_path1)
df2 = pd.read_excel(file_path2)
df3 = pd.read_excel(file_path3)

# Set 'Method' as the index for each DataFrame to make combining simpler
df1.set_index('Method', inplace=True)
df2.set_index('Method', inplace=True)
df3.set_index('Method', inplace=True)

# Combine DataFrames with df3 as priority, then df2, then df1
combined_df = df3.combine_first(df2).combine_first(df1).reset_index()

# Filter based on specified criteria
filtered_df = combined_df[(combined_df["autodiff"].isin(["Any[false, Val{:forward}]"])) & 
                          (combined_df["linsolve"] == "LUFactorization")]

# Provided name labels for line groups (e.g., "Solver" or any other relevant category)
name_labels = filtered_df["Solver"]

# Sort unique names alphabetically and create a colormap for different lines
unique_names = sorted(set(name_labels))
num_groups = len(unique_names)
cmap = cm.get_cmap('tab10', num_groups)  # Choose a colormap

# Assign a color to each unique name
color_map = {name: cmap(i) for i, name in enumerate(unique_names)}

# Plot 1: Line plot with x-axis on log scale
plt.figure(figsize=(10, 8))
for name in unique_names:
    # Filter and sort data for each group
    group_data = filtered_df[filtered_df["Solver"] == name].sort_values(by="abstol")
    plt.plot(group_data["abstol"], group_data["Avg Time (s)"], label=name, color=color_map[name], marker='o', linestyle='-')

# Set x-axis to log scale
plt.xscale("log")

# Title, axis labels, and legend
plt.title('Line Plot of tol vs. Avg Time (s) for Selected autodiff and linsolve Criteria')
plt.xlabel('Tolerance (log scale)')
plt.ylabel('Average Time (s)')
plt.legend(title="Solver Groups", bbox_to_anchor=(1.0, 1), loc='upper left')  # Show legend with group names

plt.show()

# Plot 2: Log-log plot with both x-axis and y-axis on log scale
plt.figure(figsize=(10, 8))
for name in unique_names:
    # Filter and sort data for each group
    group_data = filtered_df[filtered_df["Solver"] == name].sort_values(by="abstol")
    plt.plot(group_data["abstol"], group_data["Avg Time (s)"], label=name, color=color_map[name], marker='o', linestyle='-')

# Set both x-axis and y-axis to log scale
plt.xscale("log")
plt.yscale("log")

# Title, axis labels, and legend
plt.title('Log-Log Plot of tol vs. Avg Time (s) for Selected autodiff and linsolve Criteria')
plt.xlabel('Tolerance (log scale)')
plt.ylabel('Average Time (s) (log scale)')
plt.legend(title="Solver Groups", bbox_to_anchor=(1.0, 1), loc='upper left')  # Show legend with group names

plt.show()
