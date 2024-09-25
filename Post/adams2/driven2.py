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

def plot_each_column_with_generated_time(fn, data_csv_file, output_folder, colors, line_styles, dpi, time_threshold=None, end_time_threshold=None):
    textstyle = {'family': 'Times New Roman', 'size': 16}
    
    df = pd.read_csv(data_csv_file)
    
    # Filter the data based on the time threshold
    if time_threshold is not None:
        df = df[df['t '] >= time_threshold].copy()
    
    if end_time_threshold is not None:
        df = df[df['t '] <= end_time_threshold].copy()
    
    if time_threshold is not None:
        df['t '] -= time_threshold  # Adjust time to start from 0
    
    num_colors = len(colors)
    num_styles = len(line_styles)
    
    for i, column in enumerate(df.columns):
        if column != 't ':
            color_index = (i // num_styles) % num_colors
            style_index = i % num_styles
            
            fig, ax = plt.subplots(figsize=(12, 8))
            
            ax.plot(df['t '], df[column], label=column,
                    color=colors[color_index], linestyle=line_styles[style_index])
    
            plt.xlabel(r"$\it{t}$ [ $\it{s}$ ]", fontdict=textstyle)
            plt.ylabel('Fc of roller  [N]', fontdict=textstyle)
            plt.xticks(fontsize=20)
            plt.yticks(fontsize=20)
            plt.grid(True)
            plt.legend()
    
            # Save the plot
            output_path = os.path.join(output_folder, f'{fn}_{column}.png')
            plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
            plt.close()
            print(f"Plot saved at: {output_path}")

def run_plotting(fn, data_csv_file, output_folder, dpi=300, time_threshold=None, end_time_threshold=None):
    os.makedirs(output_folder, exist_ok=True)
    selected_colors = default_colors
    selected_styles = default_line_styles
    plot_each_column_with_generated_time(fn, data_csv_file, output_folder, selected_colors, selected_styles, dpi, time_threshold, end_time_threshold)

if __name__ == "__main__":
    matplotlib.rcParams['font.family'] = 'Times New Roman'
    data_csv_file = r"C:\OneDrive\Articles\10.Working\[D21][20211009]ContactMechanics\MBD.jl\plots\adams2\MR_door (run 29)_out2.csv"  # Replace with your actual data CSV file path
    output_folder = './plots/adams2/driv/'  # Replace with your desired output folder path
    time_threshold = 98.0  # Replace with your desired start time threshold
    end_time_threshold = 100.0  # Replace with your desired end time threshold
    run_plotting('driv', data_csv_file, output_folder, dpi=300, time_threshold=time_threshold, end_time_threshold=end_time_threshold)
