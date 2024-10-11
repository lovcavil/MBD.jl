
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
    
    # Add other columns and their scaling factors here
}

def filter_data1(df1, start_time1, end_time1, selected_columns1, data_amount):
    filtered_df1 = df1.loc[(df1['t '] >= start_time1) & (df1['t '] <= end_time1)]
    selected_columns1 = selected_columns1[:data_amount]
    return filtered_df1, selected_columns1

def filter_data2( df2, start_time2, end_time2,  selected_columns2, data_amount):
    filtered_df2 = df2.loc[(df2['Time'] >= start_time2) & (df2['Time'] <= end_time2)]
    selected_columns2 = selected_columns2[:data_amount]
    return  filtered_df2,  selected_columns2


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
            fft_vals[window_size:] = 0  # Zero out high frequencies
            filtered_signal = np.real(ifft(fft_vals))
            filtered_df1.loc[:, column] = filtered_signal
    return filtered_df1



def create_figure1_matplotlib(ax, filtered_df1, selected_columns1, start_time1, label, color='blue'):
    textstyle = {'family': 'Times New Roman', 'size': 32}
    tick_fontsize = 20
    tick_font = 'Times New Roman'
    
    for column in selected_columns1:
        ax.plot(filtered_df1['t '] - start_time1, filtered_df1[column], label=label, color=color)
    
    ax.set_title('Plot of Selected Columns vs Time', fontsize=20, fontweight='bold')
    ax.set_xlabel('Time (s)', fontdict=textstyle)
    ax.set_ylabel('Value', fontdict=textstyle)
    
    ax.grid(True)
    
    plt.xticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)
    
    plt.tight_layout()
    
    return ax

def create_figure2_matplotlib(ax, filtered_df2, selected_columns2, label, color='orange'):
    textstyle = {'family': 'Times New Roman', 'size': 32}
    tick_fontsize = 20
    tick_font = 'Times New Roman'
    
    for column in selected_columns2:
        ax.plot(filtered_df2['Time'], filtered_df2[column], label=label, color=color)
    
    ax.legend(fontsize=20, loc='upper left')  # Customize legend as needed
    
    return ax


def create_figure_matplotlib():
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.set_facecolor('white')
    return fig, ax
def apply_scaling(filtered_df1, selected_columns1):
    for column in selected_columns1:
        scale = scaling_dict.get(column, 1)
        filtered_df1.loc[:, column] = filtered_df1.loc[:, column] / scale
        print(scale)
    return filtered_df1


# Define paths
path1 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\MR_door (run 29)_out2.csv"
path2 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\ic3_exprcontace\GFO_MF U3_export.csv"
# path3 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\ic4\GFO_MF U3_export.csv"
# path4 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\ic3\GCONMR U6_export.csv"
# idea data
path5 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\ic4\GCONMR U6_export.csv"
path5 = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\perfectA18\perfect.csv"
# Load data
df1= load_data(path1)
df2= load_data(path2)
# df3= load_data(path3)
# df4= load_data(path4)
df5= load_data(path5)

# Define parameters
selected_columns1 = ['Strain@B3_2.RN_6']
start_time1 = 98.4
end_time1 = 98.4+1.518
selected_columns2 = ['clea3.asc']  # Replace with actual column names
selected_columns2 = ['clea3_07-29-2024_1.asc'] 

# selected_columns3 = ['clea4.asc']  # Replace with actual column names
# selected_columns4 = ['idea3.asc']  # Replace with actual column names
selected_columns5 = ['idea4.asc']  # Replace with actual column names
selected_columns5 = ['FY']
start_time2 = 0
end_time2 = 2
start_time3 = 0
end_time3 = 1.46
data_amount = 8
smoothing_method = 'Fourier Filter'
smoothing_method = 'Rolling Mean'
#smoothing_method = ''
window_sizes1 = [1]*9
window_sizes2 = [10]*9
# window_sizes3 = [15]*9
# window_sizes4 = [20]*9
window_sizes5 = [2]*9

# Filter data
filtered_df1, selected_columns1 = filter_data1(df1, start_time1, end_time1,  selected_columns1, data_amount)
filtered_df2, selected_columns2 = filter_data2(df2, start_time2, end_time2,  selected_columns2, data_amount)
# filtered_df3, selected_columns3 = filter_data2(df3, start_time2, end_time2,  selected_columns3, data_amount)
# filtered_df4, selected_columns4 = filter_data2(df4, start_time2, end_time2,  selected_columns4, data_amount)
filtered_df5, selected_columns5 = filter_data2(df5, start_time3, end_time3,  selected_columns5, data_amount)
# Apply scaling
filtered_df1 = apply_scaling(filtered_df1, selected_columns1)

# Apply smoothing
filtered_df1 = apply_smoothing(filtered_df1, selected_columns1,  smoothing_method, window_sizes1)
filtered_df2 = apply_smoothing(filtered_df2, selected_columns2, smoothing_method,  window_sizes2)
# filtered_df3 = apply_smoothing(filtered_df3, selected_columns3, smoothing_method,  window_sizes2)
# filtered_df4 = apply_smoothing(filtered_df4, selected_columns4, smoothing_method,  window_sizes4)
filtered_df5 = apply_smoothing(filtered_df5, selected_columns5, smoothing_method,  window_sizes5)
# Create figure
fig,ax = create_figure_matplotlib()

ax = create_figure1_matplotlib(ax, filtered_df1, selected_columns1, start_time1,"Experiment", color='#264653')


ax = create_figure2_matplotlib(ax, filtered_df5, selected_columns5,"Simulation for Perfect joint", color='#e76f51')
ax = create_figure2_matplotlib(ax, filtered_df2, selected_columns2,"Simulation for Clearance joint", color='#2a9d8f')
#ax = create_figure2_matplotlib(ax, filtered_df3, selected_columns3)
#ax = create_figure2_matplotlib(ax, filtered_df4, selected_columns4)
# Show the figure (or save it)
fig.show()  # Or use fig.write_image("output.png") to save the figure as an image
output_folder=r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\Post\A18"
output_path = os.path.join(output_folder, f'comp.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')