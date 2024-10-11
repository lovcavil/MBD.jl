import pandas as pd
import matplotlib.pyplot as plt
import os

# Read the CSV data
df = pd.read_csv('D:/OneDrive/Articles/10.Working/A18.20240920.PerfectJoint/timing/a.csv')

# Group data by 'property' and calculate mean for test1 and test2
mean_data = df.groupby('x').mean().reset_index()

# Plotting
plt.figure(figsize=(10, 5))
textstyle = {'family': 'Times New Roman', 'size': 28}
tick_fontsize = 20
tick_font = 'Times New Roman'

# Plot test1 (mean curve)
plt.plot(mean_data['x'], mean_data['data_matlab'], marker='o', label='Matlab Mean', color='#1f77b4')

# Plot test2 (mean curve)
plt.plot(mean_data['x'], mean_data['data_julia'], marker='o', label='Julia Mean', color='#ff7f0e')

# Scatter plot for original data points
plt.scatter(df['x'], df['data_matlab'], color='#1f77b4', alpha=0.7, label='Matlab Run Cases', marker='x')
plt.scatter(df['x'], df['data_julia'], color='#ff7f0e', alpha=0.7, label='Julia Run Cases', marker='x')

# Adding titles and labels
plt.xlabel('Exponent of tolerance', fontdict=textstyle)
plt.ylabel('CPU runtime', fontdict=textstyle)

# Modify xticks to add "-" before each value
modified_xticks = ['-' + str(x) for x in mean_data['x']]
plt.xticks(mean_data['x'], modified_xticks, fontsize=tick_fontsize, fontfamily=tick_font)
plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)

# Increase legend font size
plt.legend(fontsize=18)

# Show grid and save plot
plt.grid()
output_folder = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\Post\A18"
output_path = os.path.join(output_folder, 'implcomp.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
