import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from scipy.fft import fft, ifft

# Default colors and line styles
default_colors = [
    (16/256, 70/256, 128/256),
    (49/256, 124/256, 183/256),
    (109/256, 173/256, 209/256),
    (182/256, 215/256, 232/256),
    (246/256, 178/256, 147/256),
    (220/256, 109/256, 87/256),
    (183/256, 34/256, 48/256)
]
default_line_styles = ['-', '--', '-.', ':']

def load_data(path):
    return pd.read_csv(path)

def filter_data(df, start_time, end_time):
    return df.loc[(df['Time'] >= start_time) & (df['Time'] <= end_time)]

def apply_scaling(filtered_df, selected_columns, scaling_dict):
    for column in selected_columns:
        scale = scaling_dict.get(column, 1)
        filtered_df.loc[:, column] = filtered_df.loc[:, column] / scale
    return filtered_df

def apply_smoothing(filtered_df, selected_columns, smoothing_method, window_sizes):
    for column, window_size in zip(selected_columns, window_sizes):
        original_data = filtered_df[column].copy()  # Save original data for comparison
        if smoothing_method == 'Rolling Mean':
            filtered_df.loc[:, column] = filtered_df.loc[:, column].rolling(window=window_size, min_periods=1).mean()
        elif smoothing_method == 'Fourier Filter':
            signal = filtered_df[column].values
            fft_vals = fft(signal)
            fft_vals[window_size:] = 0
            filtered_signal = np.real(ifft(fft_vals))
            filtered_df.loc[:, column] = filtered_signal
        
        # Debugging: print original vs smoothed
        print(f"Original vs Smoothed for {column}:")
        print("Original:", original_data.head())
        print("Smoothed:", filtered_df[column].head())
        
    return filtered_df

def process_data(path, start_time, end_time, scaling_dict, smoothing_method, window_sizes):
    df = load_data(path)
    # Filter data
    filtered_df = filter_data(df, start_time, end_time)

    # Get all columns except "Time"
    selected_columns = filtered_df.columns.drop('Time')

    # Reflect data for specific time range (if necessary)
    reflected_df = filtered_df[(filtered_df['Time'] > 0.1) & (filtered_df['Time'] < 0.2)].copy()
    reflected_df['Time'] = 0.1 + (0.1 - reflected_df['Time'])
    filtered_df = filtered_df[(filtered_df['Time'] > 0.1)]
    filtered_df = pd.concat([filtered_df, reflected_df]).sort_values(by='Time').reset_index(drop=True)

    # Apply scaling
    filtered_df = apply_scaling(filtered_df, selected_columns, scaling_dict)

    # Apply smoothing
    filtered_df = apply_smoothing(filtered_df, selected_columns, smoothing_method, window_sizes)

    return filtered_df

def plot_all_columns_with_generated_time(fn, filtered_df, output_folder, colors, line_styles, dpi):
    textstyle = {'family': 'Times New Roman', 'size': 16}
    
    time_values = filtered_df['Time'].values  # Assuming 'Time' is a column in filtered_df
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    num_colors = len(colors)
    num_styles = len(line_styles)
    
    for i, column in enumerate(filtered_df.columns):
        if column != 'Time':
            color_index = (i // num_styles) % num_colors
            style_index = i % num_styles
            ax.plot(time_values, filtered_df[column], label=column,
                    color=colors[color_index], linestyle=line_styles[style_index])
    
    plt.xlabel(r"$\it{t}$ [ $\it{s}$ ]", fontdict=textstyle)
    plt.ylabel('Fc of roller  [N]', fontdict=textstyle)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)

    # Save the plot
    output_path = os.path.join(output_folder, f'{fn}.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Plot saved at: {output_path}")

def run_plotting(fn, data_csv_file, output_folder, scaling_dict, smoothing_method, window_sizes, dpi=300):
    os.makedirs(output_folder, exist_ok=True)
    
    # Process data
    start_time = 0
    end_time = 2
    filtered_df = process_data(data_csv_file, start_time, end_time, scaling_dict, smoothing_method, window_sizes)

    # Plot data
    selected_colors = default_colors
    selected_styles = default_line_styles
    plot_all_columns_with_generated_time(fn, filtered_df, output_folder, selected_colors, selected_styles, dpi)

if __name__ == "__main__":
    matplotlib.rcParams['font.family'] = 'Times New Roman'
    
    output_folder = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\Post\A16\f/"  # Replace with your desired output folder path

    scaling_dict = {
        # Define your scaling factors here if needed
    }
    smoothing_method = 'Rolling Mean'
    smoothing_method = ""
    window_sizes = [12]  # Adjust as necessary

    # File paths for plotting
    file_paths = [
        # r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_LG U4_export.csv",
        # r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_LG U3_export.csv",
        # r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_MG U4_export.csv",
        # r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_MG U3_export.csv",
        # r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_MF U3_export.csv",
        # r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_LF U3_export.csv"，
        r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\request_ACC_door U2_export.csv",
        r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\request_ACC_door U3_export.csv",
    ]
    
    for idx, data_csv_file in enumerate(file_paths):
        run_plotting(f'Plot_{idx + 1}', data_csv_file, output_folder, scaling_dict, smoothing_method, window_sizes)
