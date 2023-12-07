import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("1slider/slider.tab", sep='\t')

# Plotting
def d(label,pname):
    plt.figure(figsize=(16,10))
    plt.plot(df['.MODEL_1.Last_Run.PART_4_XFORM.TIME'],df[label], label=label)
    plt.savefig(f"{pname}/{label}.png")  # Saves the plot as a PNG file
    plt.close()  # Close the plot to free up memory
    
    

pname="1slider"
for label in df.columns:
    d(label,pname)


# fig=plt.figure(figsize=(16,10))
# ax = fig.add_subplot(111, projection='3d')

# # Plot the 3D curve
# ax.plot(df['x1'], df['y1'], df['z1'])

# # Set labels
# ax.set_xlabel('X axis')
# ax.set_ylabel('Y axis')
# ax.set_zlabel('Z axis')

# plt.savefig(f"{pname}/3d1.png")  # Saves the plot as a PNG file
# plt.close()  # Close the plot to free up memory

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