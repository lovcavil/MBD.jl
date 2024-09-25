import os
import pandas as pd
import matplotlib.pyplot as plt

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
    df = pd.read_csv(data_csv_file)
    time_values = 0.001 * df.index
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    num_colors = len(colors)
    num_styles = len(line_styles)
    
    for i, column in enumerate(df.columns):
        if column != 'Time':
            color_index = (i // num_styles) % num_colors
            style_index = i % num_styles
            ax.plot(time_values, df[column], label=column,
                    color=colors[color_index], linestyle=line_styles[style_index])
    
    #plt.title('All Columns vs Time')
    plt.xlabel('Simulation Time(s)')
    plt.ylabel('Computation Time Spended(s)')
    plt.grid(True)
    
    # Create custom legends
    color_legend_label=['Flores','Hunt-Crossley','Lankarani-Nikravesh','Herbert-McWhannell','Gonthier',"Hu-Guo",'Gharib-Hurmuzlu']
    style_legend_label=['linear', 'smooth', 'modified', 'classic']
    color_legend_handles = [plt.Line2D([0], [0], color=color, linestyle='-', label=f'{color_legend_label[i]}') for i, color in enumerate(colors)]
    style_legend_handles = [plt.Line2D([0], [0], color='black', linestyle=style, label=f'{style_legend_label[i]}') for i, style in enumerate(line_styles)]
    
    color_legend = ax.legend(handles=color_legend_handles, title='Contact Model', loc='upper left', bbox_to_anchor=(0, 1))
    style_legend = ax.legend(handles=style_legend_handles, title='Friction Model', loc='upper left', bbox_to_anchor=(0, 0.7))

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
    data_csv_file = r"C:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\6\CPU FUNCTION_MEA_CPU_export.csv"  # Replace with your actual data CSV file path
    output_folder = './plots'  # Replace with your desired output folder path
    run_plotting('all_columns_plot', data_csv_file, output_folder, dpi=300)
