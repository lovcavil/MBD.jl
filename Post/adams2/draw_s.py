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
    
    plt.figure(figsize=(12, 8))
    
    num_colors = len(colors)
    num_styles = len(line_styles)
    
    for i, column in enumerate(df.columns):
        if column != 'Time':
            color_index = (i // num_styles) % num_colors
            style_index = i % num_styles
            plt.plot(time_values, df[column], label=column,
                     color=colors[color_index], linestyle=line_styles[style_index])
    
    plt.title('All Columns vs Time')
    plt.xlabel('Time')
    plt.ylabel('Values')
    plt.legend()
    plt.grid(True)
    output_path = os.path.join(output_folder, f'{fn}.png')
    plt.savefig(output_path, dpi=dpi)
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
