import pandas as pd
import matplotlib.pyplot as plt
import os
from multiprocessing import Pool

# Create a folder named 'plots' if it does not exist
output_folder = 'D:\OneDrive\Articles/10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots/20241028a'
os.makedirs(output_folder, exist_ok=True)

# Read the CSV data into a DataFrame
file_path = './csv/Euler_run1.csv'
file_path = 'E:/New Folder/A16/con1/1.csv.2'

df = pd.read_csv(file_path)

# Function to plot and save each figure
def plot_and_save(column):
    plt.figure()
    plt.plot(df['Time'], df[column])
    plt.title(f'{column} vs Time')
    plt.xlabel('Time')
    plt.ylabel(column)
    plt.grid(True)
    plt.savefig(f'{output_folder}/{column}_vs_Time.png')
    plt.close()

# Get the columns to plot
columns_to_plot = [column for column in df.columns if column != 'Time']

# Use Pool to run the plotting in parallel
if __name__ == '__main__':
    with Pool(processes=os.cpu_count()) as pool:
        pool.map(plot_and_save, columns_to_plot)

    print(f"Plots saved in the folder: {output_folder}")
