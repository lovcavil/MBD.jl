import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

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

def plot_all_columns_with_generated_time(fn, data_csv_file, output_folder, colors, line_styles, dpi):
    textstyle = {'family': 'Times New Roman', 'size': 16}
    
    df = pd.read_csv(data_csv_file)
    time_values = 0.001 * df.index
    
    window_size=10
    df.iloc[:, 1:] = df.iloc[:, 1:].rolling(window=window_size, min_periods=1).mean()
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    num_colors = len(colors)
    num_styles = len(line_styles)
    
    for i, column in enumerate(df.columns):
        if column != 'Time':
            color_index = (i // num_styles) % num_colors
            style_index = i % num_styles
            ax.plot(time_values, df[column], label=column,
                    color=colors[color_index], linestyle=line_styles[style_index])
    
    plt.xlabel(r"$\it{t}$ [ $\it{s}$ ]", fontdict=textstyle)
    plt.ylabel('Fc of roller  [N]', fontdict=textstyle)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    
    # Create custom legends
    color_legend_label = ['Flores', 'Hunt-Crossley', 'Lankarani-Nikravesh', 'Herbert-McWhannell', 'Gonthier', "Hu-Guo", 'Gharib-Hurmuzlu']
    style_legend_label = ['Linear', 'Smooth', 'Modified', 'Classic']
    color_legend_handles = [plt.Line2D([0], [0], color=color, linestyle='-', label=f'{color_legend_label[i]}') for i, color in enumerate(colors)]
    style_legend_handles = [plt.Line2D([0], [0], color='black', linestyle=style, label=f'{style_legend_label[i]}') for i, style in enumerate(line_styles)]
    
    color_legend = ax.legend(handles=color_legend_handles, title='Contact Model', loc='upper left', bbox_to_anchor=(0.05, 0.45), prop={'family': 'Times New Roman', 'size': 16}, title_fontsize=16)
    style_legend = ax.legend(handles=style_legend_handles, title='Friction Model', loc='upper left', bbox_to_anchor=(.36, 0.3), prop={'family': 'Times New Roman', 'size': 16}, title_fontsize=16)

    # Add the legends to the plot
    ax.add_artist(color_legend)
    ax.add_artist(style_legend)

    # Save the plot
    output_path = os.path.join(output_folder, f'{fn}.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Plot saved at: {output_path}")

def run_plotting(fn, data_csv_file, output_folder, dpi=300):
    os.makedirs(output_folder, exist_ok=True)
    selected_colors = default_colors
    selected_styles = default_line_styles
    plot_all_columns_with_generated_time(fn, data_csv_file, output_folder, selected_colors, selected_styles, dpi)

if __name__ == "__main__":
    matplotlib.rcParams['font.family'] = 'Times New Roman'
    data_csv_file = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_LG U4_export.csv"  # Replace with your actual data CSV file path
    output_folder = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\Post\A16\f/"  # Replace with your desired output folder path
    run_plotting('LG', data_csv_file, output_folder, dpi=300)

    data_csv_file = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_LG U3_export.csv"  # Replace with your actual data CSV file path
    run_plotting('LG_Y', data_csv_file, output_folder, dpi=300)
    
    data_csv_file = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_MG U4_export.csv"  # Replace with your actual data CSV file path
    run_plotting('MG', data_csv_file, output_folder, dpi=300)

    data_csv_file = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_MG U3_export.csv"  # Replace with your actual data CSV file path
    run_plotting('MG_Y', data_csv_file, output_folder, dpi=300)
    
    data_csv_file = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_MF U3_export.csv"  # Replace with your actual data CSV file path
    run_plotting('MF', data_csv_file, output_folder, dpi=300)
    
    data_csv_file = r"D:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\7\GFO_LF U3_export.csv"  # Replace with your actual data CSV file path
    run_plotting('LF', data_csv_file, output_folder, dpi=300)
# 