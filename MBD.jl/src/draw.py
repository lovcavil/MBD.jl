import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("1111/data.csv")

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