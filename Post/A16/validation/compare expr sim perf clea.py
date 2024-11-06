import pandas as pd
import numpy as np
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
import os

def load_data(path):
    df = pd.read_csv(path)
    return df

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

def filter_data1(df1, start_time1, end_time1, selected_columns1, data_amount):
    filtered_df1 = df1.loc[(df1['t '] >= start_time1) & (df1['t '] <= end_time1)]
    selected_columns1 = selected_columns1[:data_amount]
    return filtered_df1, selected_columns1

def filter_data2(df2, start_time2, end_time2, selected_columns2, data_amount):
    filtered_df2 = df2.loc[(df2['Time'] >= start_time2) & (df2['Time'] <= end_time2)]
    selected_columns2 = selected_columns2[:data_amount]
    return filtered_df2, selected_columns2

def apply_scaling(filtered_df1, selected_columns1):
    for column in selected_columns1:
        scale = scaling_dict.get(column, 1)
        filtered_df1.loc[:, column] = filtered_df1.loc[:, column] / scale
        print(scale)
    return filtered_df1

def apply_smoothing(filtered_df1, selected_columns1, smoothing_method, window_sizes1):
    if smoothing_method == 'Rolling Mean':
        for column, window_size in zip(selected_columns1, window_sizes1):
            filtered_df1.loc[:, column] = filtered_df1.loc[:, column].rolling(window=window_size, min_periods=1).mean()
    elif smoothing_method == 'Fourier Filter':
        for column, window_size in zip(selected_columns1, window_sizes1):
            signal = filtered_df1[column].values
            fft_vals = fft(signal)
            fft_vals[window_size:] = 0
            filtered_signal = np.real(ifft(fft_vals))
            filtered_df1.loc[:, column] = filtered_signal
    return filtered_df1

def create_figure1_matplotlib(ax, filtered_df1, selected_columns1, start_time1, label, color='blue', line_style='-'):
    textstyle = {'family': 'Times New Roman', 'size': 32}
    tick_fontsize = 20
    tick_font = 'Times New Roman'
    
    for column in selected_columns1:
        ax.plot(filtered_df1['t '] - start_time1, filtered_df1[column], label=label, color=color, linestyle=line_style)
    
    #ax.set_title('Plot of Selected Columns vs Time', fontsize=20, fontweight='bold')
    ax.set_xlabel('Time (s)', fontdict=textstyle)
    ax.set_ylabel('Normal Contact Force(N)', fontdict=textstyle)
    #ax.grid(True)
    plt.xticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.tight_layout()
    
    return ax

def create_figure2_matplotlib(ax, filtered_df2, selected_columns2, label, color='orange', line_style='--'):
    textstyle = {'family': 'Times New Roman', 'size': 32}
    tick_fontsize = 20
    tick_font = 'Times New Roman'
    
    for column in selected_columns2:
        ax.plot(filtered_df2['Time'], filtered_df2[column], label=label, color=color, linestyle=line_style)
    
    ax.legend(fontsize=20, loc='upper left')
    
    return ax

def create_figure_matplotlib():
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.set_facecolor('white')
    return fig, ax

# Define paths
path1 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\MR_door (run 29)_out2.csv"
path2 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\ic3_exprcontace\GFO_MF U3_export.csv"
path5 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\perfectA18\perfect.csv"

# Load data
df1 = load_data(path1)
df2 = load_data(path2)
df5 = load_data(path5)

# Define parameters
selected_columns1 = ['Strain@B3_2.RN_6']
start_time1 = 98.4
end_time1 = 98.4 + 1.518
selected_columns2 = ['clea3_07-29-2024_1.asc']
selected_columns5 = ['FY']
start_time2 = 0
end_time2 = 2
start_time3 = 0
end_time3 = 1.46
data_amount = 8
smoothing_method = 'Rolling Mean'
window_sizes1 = [1] * 9
window_sizes2 = [12] * 9
window_sizes5 = [2] * 9

# Filter data
filtered_df1, selected_columns1 = filter_data1(df1, start_time1, end_time1, selected_columns1, data_amount)
filtered_df2, selected_columns2 = filter_data2(df2, start_time2, end_time2, selected_columns2, data_amount)
filtered_df5, selected_columns5 = filter_data2(df5, start_time3, end_time3, selected_columns5, data_amount)

# Reflect data for t = 0.1 to t = 0.2 in filtered_df2
reflected_df = filtered_df2[(filtered_df2['Time'] > 0.1) & (filtered_df2['Time'] < 0.2)].copy()
reflected_df['Time'] = 0.1 + (0.1 - reflected_df['Time'])
filtered_df2 = filtered_df2[(filtered_df2['Time'] > 0.1)]
filtered_df2 = pd.concat([filtered_df2, reflected_df]).sort_values(by='Time').reset_index(drop=True)

# Apply scaling
filtered_df1 = apply_scaling(filtered_df1, selected_columns1)

# Apply smoothing
filtered_df1 = apply_smoothing(filtered_df1, selected_columns1, smoothing_method, window_sizes1)
filtered_df2 = apply_smoothing(filtered_df2, selected_columns2, smoothing_method, window_sizes2)
filtered_df5 = apply_smoothing(filtered_df5, selected_columns5, smoothing_method, window_sizes5)

# Create figure
fig, ax = create_figure_matplotlib()

# Plot with solid and dashed lines
ax = create_figure1_matplotlib(ax, filtered_df1, selected_columns1, start_time1, "Experiment", color='k', line_style='--')
ax = create_figure2_matplotlib(ax, filtered_df2, selected_columns2, "Simulation", color='k', line_style='-')

# Show the figure (or save it)
fig.show()
output_folder = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\Post\A16\validation"
output_path = os.path.join(output_folder, 'comp.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')
