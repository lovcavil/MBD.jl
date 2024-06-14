import pandas as pd
import matplotlib.pyplot as plt

# Load your CSV file
data = pd.read_csv('./SavedResult/solvercompare/COMPARE.csv')

# Extract x and y data for the scatter plot
x_data = data.iloc[:, 1]  # Second column for x-axis
y_data = data.iloc[:, 2]  # Third column for y-axis

# Using the first column for labels
labels = data.iloc[:, 0]  # First column for labels

plt.figure(figsize=(10, 6))

# Create a scatter plot
plt.scatter(x_data, y_data)

# Annotate each point with its corresponding label from the first column
for i, label in enumerate(labels):
    plt.annotate(label, (x_data[i], y_data[i]), textcoords="offset points", xytext=(0,10), ha='center')

plt.title('Scatter Plot with Labels from CSV Data')
plt.xlabel(data.columns[1])  # Label x-axis with the name of the second column
plt.ylabel(data.columns[2])  # Label y-axis with the name of the third column
plt.savefig(f'./SavedResult\solvercompare/RSMEvTS.png', dpi=500)
plt.show()