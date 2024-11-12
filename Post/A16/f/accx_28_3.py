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

folder_path = os.getenv("ONEDRIVE")

path1 = os.path.join(folder_path,r"Articles\10.Working\MBD.jl\plots\adams2\MR_door (run 29)_out2.csv")
path2 = os.path.join(folder_path,r"Articles\10.Working\MBD.jl\plots\adams2\7\request_ACC_door U2_export.csv")

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
end_time2 = 1.507

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
    line_styles = ['-']
    line_colors = [
        (16/256, 70/256, 128/256),
        (49/256, 124/256, 183/256),
        # (109/256, 173/256, 209/256),
        # (182/256, 215/256, 232/256),
        # (246/256, 178/256, 147/256),
        (220/256, 109/256, 87/256),
        (183/256, 34/256, 48/256)
    ]

    styles_and_colors = {}
    
    for index, col in enumerate(selected_columns):
        # Determine line style and color based on some rule, here we use index for simplicity
        styles_and_colors[col] = {
            'color': line_colors[index % len(line_colors)],
            'linestyle': line_styles[index % len(line_styles)]
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
window_sizes2 = [5] * len(selected_columns2)

#filtered_df1 = apply_smoothing(filtered_df1, selected_columns1, smoothing_method, window_sizes1)
filtered_df1 = filtered_df1
filtered_df2 = apply_smoothing(filtered_df2, selected_columns2, smoothing_method, window_sizes2)

processed_data1 = process_data(filtered_df1, start_time1, timename1)
processed_data2 = process_data(filtered_df2, start_time2, timename2)

# Adjust the processed data
processed_data1['acc@A7_X.RN_6'] *= 9.8
for col in selected_columns2:
    processed_data2[col] /= 1000

def create_and_save_figures(processed_data1, processed_data2, selected_columns2):
    textstyle = {'family': 'Times New Roman', 'size': 24}
    tick_fontsize = 20
    tick_font = 'Times New Roman'
    fig_label = ['Flores', 'Hunt-Crossley', 'Lankarani-Nikravesh', 
                 'Herbert-McWhannell', 'Gonthier', "Hu-Guo", 'Gharib-Hurmuzlu']
    style_label = ['Linear', 'Smooth', 'Modified', 'Classic']

    # Group by third character
    third_grouped = {}
    for filename in selected_columns2:
        third_char = filename[2]  # Get the third character
        if third_char not in third_grouped:
            third_grouped[third_char] = []
        third_grouped[third_char].append(filename)

    # Prepare for subplots: 4 rows, 2 columns
    num_plots = len(third_grouped)
    nrows, ncols = 4, 2
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16, 24))
    axs = axs.flatten()  # Flatten the 2D array of axes to easily iterate over

    # Generate figures for each group by third character
    for index, (third_char, columns) in enumerate(third_grouped.items()):
        ax = axs[index]  # Get the appropriate subplot axis

        if columns:  # Only create a figure if there are columns to plot
            styles_and_colors = assign_styles_and_colors(columns)
            for col_index, col in enumerate(columns):
                style = styles_and_colors.get(col, {})
                na = f"{style_label[col_index]}"
                ax.plot(processed_data2['t_adj'], processed_data2[col], label=f'Simulation: {na} Friction ', 
                        color=style.get('color'), linestyle=style.get('linestyle'), linewidth=1)

                # Annotate peak values in the specified time range
                peak_mask = (processed_data2['t_adj'] >= 1.2) & (processed_data2['t_adj'] <= 1.3)
                if peak_mask.any():
                    peak_value = processed_data2[col][peak_mask].max()
                    peak_time = processed_data2['t_adj'][peak_mask][processed_data2[col][peak_mask].idxmax()]
                    # ax.annotate(f'Peak: {peak_value:.2f}', 
                    #             xy=(peak_time, peak_value), 
                    #             xytext=(peak_time + 0.01, peak_value + 0.5), 
                    #             arrowprops=dict(facecolor='black', arrowstyle='->'),
                    #             fontsize=8, color=style.get('color'))
                    
            ax.plot(processed_data1['t_adj'], processed_data1['acc@A7_X.RN_6'], 
                    label='Experiment', color='k', linewidth=1, linestyle='--')

            # ax.set_title(f"Grouped by Third Character: {third_char}", fontsize=20)
            ax.set_xlabel("t (s)", fontdict=textstyle)
            ax.set_ylabel("Acceleration in x direction (m/s²)", fontdict=textstyle)
            
            # Set tick labels for each axis
            ax.tick_params(axis='both', labelsize=tick_fontsize)  # Set label size
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontname(tick_font)  # Set font family for each tick label

            ax.legend(fontsize=14)
    
    # Remove any unused subplots if necessary
    for remaining_ax in axs[num_plots:]:
        remaining_ax.remove()

    # Save the combined figure
    output_path = os.path.join(os.path.dirname(__file__), 'grouped_subplots.png')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)

# Generate plots with peak annotations
create_and_save_figures(processed_data1, processed_data2, selected_columns2)


