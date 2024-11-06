import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm

# Define the path to your Excel files
file_path1 = r'E:/New Folder/A16/tol2.xlsx'
file_path2 = r'E:/New Folder/A16/tol3.xlsx'
file_path3 = r'E:/New Folder/A16/tol4.xlsx'

# Load the DataFrames
df1 = pd.read_excel(file_path1)
df2 = pd.read_excel(file_path2)
df3 = pd.read_excel(file_path3)

# Set 'Method' as the index for each DataFrame
df1.set_index('Method', inplace=True)
df2.set_index('Method', inplace=True)
df3.set_index('Method', inplace=True)
Solver="QBDF"
Solver="FBDF"
#Solver="TRBDF2"
# Filter each DataFrame for the solver "QBNF" with specified criteria
filtered_df1 = df1[(df1["Solver"] == Solver) & (df1["autodiff"] == "Any[false, Val{:forward}]") & (df1["linsolve"] == "LUFactorization")]
filtered_df2 = df2[(df2["Solver"] == Solver) & (df2["autodiff"] == "Any[false, Val{:forward}]") & (df2["linsolve"] == "LUFactorization")]
filtered_df3 = df3[(df3["Solver"] == Solver) & (df3["autodiff"] == "Any[false, Val{:forward}]") & (df3["linsolve"] == "LUFactorization")]

# Initialize the plot
plt.figure(figsize=(10, 8))

# Define colors for each DataFrame to distinguish them in the plot
colors = ['blue', 'green', 'red']
labels = ['tol2.xlsx', 'tol3.xlsx', 'tol4.xlsx']

# Plot each filtered DataFrame for solver "QBNF" on the same figure
for i, (filtered_df, color, label) in enumerate(zip([filtered_df1, filtered_df2, filtered_df3], colors, labels)):
    filtered_df = filtered_df.sort_values(by="abstol")  # Sort by tolerance if not sorted
    plt.plot(filtered_df["abstol"], filtered_df["Avg Time (s)"], label=label, color=color, marker='o', linestyle='-')

# Set x-axis to log scale
plt.xscale("log")
plt.yscale("log")

# Title, axis labels, and legend
plt.title('Log-Log Plot of Tolerance vs. Avg Time (s) for Solver "QBNF" Across Different Files')
plt.xlabel('Tolerance (log scale)')
plt.ylabel('Average Time (s) (log scale)')
plt.legend(title="Data Sources", bbox_to_anchor=(1.0, 1), loc='upper left')  # Show legend with file labels

plt.show()
