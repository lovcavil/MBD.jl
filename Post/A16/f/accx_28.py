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

timename1 = "t "
timename2 = "Time"

# Define parameters
selected_columns1 = ['acc@A7_X.RN_6']
start_time1 = 98.415
end_time1 = 99.930
selected_columns2 = df2.columns.tolist()
selected_columns2.remove(timename2)

start_time2 = 0
end_time2 = 1.515

# Filter data
def filter_data(df, start_time, end_time, selected_columns, tn=""):
    return df[(df[tn] >= start_time) & (df[tn] <= end_time)][selected_columns + [tn]]

# Apply scaling
def apply_scaling(df, selected_columns, scaling_dict):
    for col in selected_columns:
        scale = scaling_dict.get(col, 1)
        df[col] /= scale
    return df

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

# Process Data
def process_data(filtered_df, start_time, tn):
    filtered_df['t_adj'] = filtered_df[tn] - start_time
    return filtered_df

def assign_styles_and_colors(selected_columns):
    line_styles = ['-', '--', '-.', ':']
    line_colors = [
        (16/256, 70/256, 128/256),
        (49/256, 124/256, 183/256),
        (109/256, 173/256, 209/256),
        (182/256, 215/256, 232/256),
        (246/256, 178/256, 147/256),
        (220/256, 109/256, 87/256),
        (183/256, 34/256, 48/256)
    ]

    styles_and_colors = {}
    for col in selected_columns:
        if len(col) >= 4:
            first_char = col[1]
            third_char = col[3]
            
            lc_index = int(first_char) - 2 if first_char.isdigit() else 0
            ls_index = int(third_char) - 3 if third_char.isdigit() else 0
            
            styles_and_colors[col] = {
                'color': line_colors[lc_index % len(line_colors)],
                'linestyle': line_styles[ls_index % len(line_styles)]
            }
    return styles_and_colors

# Load and process data
filtered_df1 = filter_data(df1, start_time1, end_time1, selected_columns1, timename1)
filtered_df2 = filter_data(df2, start_time2, end_time2, selected_columns2, timename2)

# Scale both filtered dataframes
filtered_df1 = apply_scaling(filtered_df1, selected_columns1, scaling_dict)
filtered_df2 = apply_scaling(filtered_df2, selected_columns2, scaling_dict)

# Apply smoothing
smoothing_method = 'Rolling Mean'
window_sizes1 = [1] * len(selected_columns1)
window_sizes2 = [1] * len(selected_columns2)

filtered_df1 = apply_smoothing(filtered_df1, selected_columns1, smoothing_method, window_sizes1)
filtered_df2 = apply_smoothing(filtered_df2, selected_columns2, smoothing_method, window_sizes2)

processed_data1 = process_data(filtered_df1, start_time1, timename1)
processed_data2 = process_data(filtered_df2, start_time2, timename2)

# Adjust the processed data
processed_data1['acc@A7_X.RN_6'] *= 9.8
for col in selected_columns2:
    processed_data2[col] /= 1000

# Plotting function to create and save figures based on character grouping
def create_and_save_figures(processed_data1, processed_data2, selected_columns2):
    # Group by first character
    first_grouped = {}
    for filename in selected_columns2:
        first_char = filename[0]  # Get the first character
        if first_char not in first_grouped:
            first_grouped[first_char] = []
        first_grouped[first_char].append(filename)

    # Generate figures for each group by first character
    for first_char, columns in first_grouped.items():
        if columns:  # Only create a figure if there are columns to plot
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(processed_data1['t_adj'], processed_data1['acc@A7_X.RN_6'], label='Experiment', color='k', linewidth=1, linestyle='--')
            
            styles_and_colors = assign_styles_and_colors(columns)
            for col in columns:
                style = styles_and_colors.get(col, {})
                ax.plot(processed_data2['t_adj'], processed_data2[col], label=f'Simulation {col}', 
                        color=style.get('color'), linestyle=style.get('linestyle'), linewidth=1)
            
            ax.set_title(f"Grouped by First Character: {first_char}")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Acceleration in x Direction (m/s²)")
            ax.legend()
            output_path = os.path.join(os.path.dirname(__file__), f'group_by_first_char_{first_char}.png')
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)

    # Group by third character
    third_grouped = {}
    for filename in selected_columns2:
        third_char = filename[2]  # Get the third character
        if third_char not in third_grouped:
            third_grouped[third_char] = []
        third_grouped[third_char].append(filename)

    # Generate figures for each group by third character
    for third_char, columns in third_grouped.items():
        if columns:  # Only create a figure if there are columns to plot
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(processed_data1['t_adj'], processed_data1['acc@A7_X.RN_6'], label='Experiment', color='k', linewidth=1, linestyle='--')
            
            styles_and_colors = assign_styles_and_colors(columns)
            for col in columns:
                style = styles_and_colors.get(col, {})
                ax.plot(processed_data2['t_adj'], processed_data2[col], label=f'Simulation {col}', 
                        color=style.get('color'), linestyle=style.get('linestyle'), linewidth=1)
            
            ax.set_title(f"Grouped by Third Character: {third_char}")
            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Acceleration in x Direction (m/s²)")
            ax.legend()
            output_path = os.path.join(os.path.dirname(__file__), f'group_by_third_char_{third_char}.png')
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)

# Generate plots
create_and_save_figures(processed_data1, processed_data2, selected_columns2)
