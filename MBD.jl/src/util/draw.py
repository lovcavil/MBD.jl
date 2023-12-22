import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("1111/data.csv")
#df = pd.read_csv("demo6/data.csv")
# Plotting
def d(label,pname):
    plt.figure(figsize=(16,10))
    plt.plot(df['t'],df[label], label=label)
    plt.savefig(f"{pname}/{label}.png")  # Saves the plot as a PNG file
    plt.close()  # Close the plot to free up memory
    
    
# Define the possible values for each part
values = ['1', '-1']
values = ['1']
# Use list comprehension to generate all combinations
combinations = [v1 + v2 + v3 + v4 for v1 in values for v2 in values for v3 in values for v4 in values]

# Print all combinations
for combination in combinations:

    for label in df.columns:
        d(label,combination)
pass
pname="1111"
fig=plt.figure(figsize=(16,10))
ax = fig.add_subplot(111, projection='3d')

# Plot the 3D curve
ax.plot(df['x1'], df['y1'], df['z1'])

# Set labels
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

plt.savefig(f"{pname}/3d1.png")  # Saves the plot as a PNG file
plt.close()  # Close the plot to free up memory

# fig=plt.figure(figsize=(16,10))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the 3D curve
# ax.plot(df['x2'], df['y2'], df['z2'])

# # Set labels
# ax.set_xlabel('X axis')
# ax.set_ylabel('Y axis')
# ax.set_zlabel('Z axis')

# plt.savefig(f"{pname}/3d2.png")  # Saves the plot as a PNG file
# plt.close()  # Close the plot to free up memory
# fig=plt.figure(figsize=(16,10))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the 3D curve
# ax.plot(df['x3'], df['y3'], df['z3'])

# # Set labels
# ax.set_xlabel('X axis')
# ax.set_ylabel('Y axis')
# ax.set_zlabel('Z axis')

# plt.savefig(f"{pname}/3d3.png")  # Saves the plot as a PNG file
# plt.close()  # Close the plot to free up memory