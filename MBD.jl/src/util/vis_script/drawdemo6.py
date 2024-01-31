import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
#df = pd.read_csv("1111/data.csv")
df = pd.read_csv("demo6/data.csv")
dfmt = pd.read_csv("demo6/t.csv")
dfmQ = pd.read_csv("demo6/Q.csv")
# Plotting
def draw_jl_res(label,pname):
    plt.figure(figsize=(16,10))
    plt.plot(df['t'],df[label], label=label)
    plt.savefig(f"{pname}/{label}.png")  # Saves the plot as a PNG file
    plt.close()  # Close the plot to free up memory
    
def draw_m_res(label,pname):
    plt.figure(figsize=(16,10))
    plt.plot(dfmt['t'],dfmQ[label], label=label)
    plt.savefig(f"{pname}/{label}m.png")  # Saves the plot as a PNG file
    plt.close()  # Close the plot to free up memory

for label in df.columns:
    draw_jl_res(label,"demo6")
for label in dfmQ.columns:
    draw_m_res(label,"demo6")    
pass

pname="demo6"
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

fig=plt.figure(figsize=(16,10))
ax = fig.add_subplot(111, projection='3d')

# Plot the 3D curve
ax.plot(df['x2'], df['y2'], df['z2'])

# Set labels
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

plt.savefig(f"{pname}/3d2.png")  # Saves the plot as a PNG file
plt.close()  # Close the plot to free up memory

#_______________________________________________________________________
fig=plt.figure(figsize=(16,10))
plt.plot(df['t'],df['z2'], label='z2-jl')
plt.plot(dfmt['t'],dfmQ['z2'], label='z2-mat')
plt.legend()
plt.savefig(f"{pname}/compare-z2.png")  # Saves the plot as a PNG file
plt.close()

