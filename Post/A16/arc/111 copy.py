import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
import os

# Scaling dictionary
scaling_dict = {
    'Strain@B1_1.RN_6': 0.339,
    'Strain@B2_1.RN_6': 0.3505,
    'Strain@B2_2.RN_6': 0.5762,
    'Strain@B2_3.RN_6': -0.5072,
    'Strain@B2_4.RN_6': 0.0763,
    'Strain@B3_1.RN_6': 0.8693,
    'Strain@B3_2.RN_6': 0.8576,
    'Strain@B3_3.RN_6': -0.2255,
    'Strain@B3_4.RN_6': 0.2455,
}

# Define paths
path1 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\MR_door (run 29)_out2.csv"
path2 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\request_ACC_door U2_export.csv"

# Load data
df1 = pd.read_csv(path1)
df2 = pd.read_csv(path2)

# Define parameters
selected_columns1 = ['acc@A7_X.RN_6', "t "]
start_time1 = 98.415
end_time1 = 99.930
selected_columns2 = ['Time', '413.asc', '414.asc']  # Include '414.asc'
start_time2 = 0
end_time2 = 1.515
data_amount = 8
smoothing_method = 'Rolling Mean'
window_sizes1 = [1] * len(selected_columns1)
window_sizes2 = [1] * len(selected_columns2)

# Filter data
def filter_data(df, start_time, end_time, selected_columns, tn=""):
    return df[(df[tn] >= start_time) & (df[tn] <= end_time)][selected_columns]

filtered_df1 = filter_data(df1, start_time1, end_time1, selected_columns1, "t ")
filtered_df2 = filter_data(df2, start_time2, end_time2, selected_columns2, "Time")

# Apply scaling
def apply_scaling(df, selected_columns, scaling_dict):
    for col in selected_columns:
        scale = scaling_dict.get(col, 1)
        df[col] /= scale
    return df

# Scale both filtered dataframes
filtered_df1 = apply_scaling(filtered_df1, selected_columns1, scaling_dict)
filtered_df2 = apply_scaling(filtered_df2, selected_columns2, scaling_dict)  # Apply scaling to 414.asc

# Apply smoothing
def apply_smoothing(df, selected_columns, method, window_sizes):
    if method == 'Rolling Mean':
        for col, win in zip(selected_columns, window_sizes):
            df[col] = df[col].rolling(window=win, min_periods=1).mean()
    elif method == 'Fourier Filter':
        for col, win in zip(selected_columns, window_sizes):
            signal = df[col].values
            fft_vals = fft(signal)
            fft_vals[win:] = 0
            df[col] = np.real(ifft(fft_vals))
    return df

filtered_df1 = apply_smoothing(filtered_df1, selected_columns1, smoothing_method, window_sizes1)
filtered_df2 = apply_smoothing(filtered_df2, selected_columns2, smoothing_method, window_sizes2)

# Process Data
def process_data(filtered_df1, filtered_df2, start_time1):
    filtered_df1['t_adj'] = filtered_df1['t '] - start_time1
    filtered_df2['t_adj'] = filtered_df2['Time']  # Adjust if needed
    return filtered_df1, filtered_df2

processed_data1, processed_data2 = process_data(filtered_df1, filtered_df2, start_time1)
processed_data1['acc@A7_X.RN_6'] = processed_data1['acc@A7_X.RN_6'] * 9.8
processed_data2['413.asc'] = processed_data2['413.asc'] / 1000 
processed_data2['414.asc'] = processed_data2['414.asc'] / 1000  # Scale 414.asc similarly

# Plotting
def create_figure(processed_data1, processed_data2, selected_columns1, selected_columns2):
    fig, ax = plt.subplots(figsize=(10, 6))
    textstyle = {'family': 'Times New Roman', 'size': 24}
    tick_fontsize = 20
    tick_font = 'Times New Roman'
    
    # Plot Experiment data
    ax.plot(processed_data1['t_adj'], processed_data1['acc@A7_X.RN_6'], label='Experiment', color='k', linewidth=1, linestyle='--')
    # Plot Simulation data for 413.asc
    ax.plot(processed_data2['t_adj'], processed_data2['413.asc'], label='Simulation 413', color='k', linewidth=1)
    # Plot Simulation data for 414.asc
    ax.plot(processed_data2['t_adj'], processed_data2['414.asc'], label='Simulation 414', color='b', linewidth=1, linestyle='--')  # Adjust color and style as needed
    
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Acceleration in x Direction (m/s²)", fontdict=textstyle)
    ax.legend(fontsize=18)
    plt.xticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)
    
    script_dir = os.path.dirname(__file__)
    output_path = os.path.join(script_dir, 'valiacc2.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.show()

# Generate plot
create_figure(processed_data1, processed_data2, selected_columns1, selected_columns2)
