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
default_line_styles = ['-', '--', '-.', ':', (0, (1, 10)), (0, (5, 10)), (0, (3, 5, 1, 5)),
                       (0, (5, 1)), (0, (3, 10, 1, 10)), (0, (3, 1, 1, 1))]

def plot_and_save(fn, csv_files, columns, output_folder, colors, line_styles, dpi, mode='single'):
    if mode == 'single':
        plt.figure(figsize=(12, 8))
        for i, file in enumerate(csv_files):
            df = pd.read_csv(file)
            for column in columns:
                if column in df.columns:
                    plt.plot(df['Time'], df[column], label=f"{os.path.basename(file)} - {column}",
                             color=colors[i % len(colors)], linestyle=line_styles[i % len(line_styles)])
        plt.title('Columns vs Time')
        plt.xlabel('Time')
        plt.ylabel('Values')
        plt.legend()
        plt.grid(True)
        output_path = os.path.join(output_folder, f'{fn}.png')
        plt.savefig(output_path, dpi=dpi)
        plt.close()

    elif mode == 'subplot':
        fig, axes = plt.subplots(len(columns), 1, figsize=(12, 8 * len(columns)))
        for i, column in enumerate(columns):
            for j, file in enumerate(csv_files):
                df = pd.read_csv(file)
                if column in df.columns:
                    axes[i].plot(df['Time'], df[column], label=f"{os.path.basename(file)} - {column}",
                                 color=colors[j % len(colors)], linestyle=line_styles[j % len(line_styles)])
            axes[i].set_title(f'{column} vs Time')
            axes[i].set_xlabel('Time')
            axes[i].set_ylabel(column)
            axes[i].legend()
            axes[i].grid(True)
        output_path = os.path.join(output_folder, f'{fn}_subplot.png')
        plt.savefig(output_path, dpi=dpi)
        plt.close()

    print(f"Plot saved at: {output_path}")

def run_plotting(fn, csv_files, columns_to_plot, output_folder, dpi=300, mode='single'):
    os.makedirs(output_folder, exist_ok=True)

    selected_colors = default_colors
    selected_styles = default_line_styles

    plot_and_save(fn, csv_files, columns_to_plot, output_folder, selected_colors, selected_styles, dpi, mode)

if __name__ == "__main__":
    csv_files = ["C:/OneDrive/Articles/10.Working/[D21][20211009]ContactMechanics/MBD.jl/csv/2024-06-20_19-32-40_tb_355__contact_g65_Tsit5_run1.csv"]  # Replace with your actual CSV file paths
    output_folder = './plots'  # Replace with your desired output folder path
    
    columns_to_plot = ['ex_12_mf_Fy', 'ex_14_mf_Fdy']  # Replace with your actual column names
    run_plotting('mf_damper', csv_files, columns_to_plot, output_folder, dpi=300, mode='single')
    
    columns_to_plot = ['ex_17_mf_vel_slide', 'ex_14_mf_Fdy']  # Replace with your actual column names
    run_plotting('mf_v', csv_files, columns_to_plot, output_folder, dpi=300, mode='subplot')
    columns_to_plot = ['ex_16_mf_vy', 'ex_14_mf_Fdy']  # Replace with your actual column names
    run_plotting('mf_vd', csv_files, columns_to_plot, output_folder, dpi=300, mode='subplot')
    columns_to_plot = ['ex_21_mf_testflag','ex_16_mf_vy','ex_8_mf_ub_y', 'ex_14_mf_Fdy']  # Replace with your actual column names
    run_plotting('mf_flagd', csv_files, columns_to_plot, output_folder, dpi=300, mode='subplot')
    columns_to_plot = ['ex_7_mf_ub_x','ex_8_mf_ub_y','ex_11_mf_Fx', 'ex_12_mf_Fy']  # Replace with your actual column names
    run_plotting('mf_', csv_files, columns_to_plot, output_folder, dpi=300, mode='subplot')