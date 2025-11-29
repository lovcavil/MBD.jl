import os
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d  # Added import

SCRIPT_DIR = Path(__file__).resolve().parent


def load_test(file):
    """加载测试数据"""
    column_names = [
        "StrainB1_1", "StrainB2_1", "StrainB2_2", "StrainB2_3", "StrainB2_4", "StrainB3_1", "StrainB3_2",
        "StrainB3_3", "StrainB3_4", "accA1_X", "accA1_Y", "accA1_Z", "accA2_X", "accA2_Y", "accA2_Z",
        "accA3_X", "accA3_Y", "accA3_Z", "accA4_X", "accA4_Y", "accA4_Z", "accA5_X", "accA5_Y", "accA5_Z",
        "accA6_X", "accA6_Y", "accA6_Z", "accA7_X", "accA7_Y", "accA7_Z", "canpress"
    ]
    column_types = {name: 'float64' for name in column_names}
    try:
        df = pd.read_csv(
            file,
            delim_whitespace=True,
            skiprows=26,
            names=column_names,
            dtype=column_types,
            engine='python'
        )
        return df
    except Exception as e:
        print(f"加载文件 {file} 时出错: {e}")
        return None


def post_test(file):
    """处理测试数据，添加时间列"""
    test = load_test(file)
    if test is not None:
        delta = 1.9531250e-3
        t = np.arange(delta, delta * (len(test) + 1), delta)
        test['t'] = t
    return test


def normalize_label_data(label_data):
    """归一化标签数据"""
    min_val = min(label_data.values())
    max_val = max(label_data.values())
    normalized_label_data = {
        key: (val - min_val) / (max_val - min_val) for key, val in label_data.items()
    }
    return normalized_label_data, min_val, max_val


def plot_smoothed_pressure(target, label_data):
    """绘制单个数据文件的平滑压力曲线"""
    dforcesme = {}
    dforcesme0 = []
    saved_plot_path = None
    plt.figure(figsize=(10, 6.5), facecolor='white')
    ax = plt.subplot(1, 1, 1)
    ax.grid(True)

    i1, i2 = target
    file_path = os.path.join(
        "C:/Vault/10.Project/A009.BYD_R2/20220710/DATA_PROC/CUT/data_allclose_22919",
        f"MR_door_{i1}{i2}.asc"
    )

    if (i1, i2) in label_data and os.path.isfile(file_path):
        test = post_test(file_path)
        if test is None:
            print(f"File {file_path} could not be loaded.")
            return dforcesme, dforcesme0, saved_plot_path

        x = test['t'].to_numpy()
        y = test['canpress'].to_numpy()
        aaa = y > 10
        V1 = np.argmax(aaa) if np.any(aaa) else None
        M = np.max(y)
        V2 = np.argmax(y)
        V3 = len(y) - np.argmax(aaa[::-1]) - 1 if np.any(aaa) else None

        if V1 is not None and V3 is not None:
            dforcesme[(i1, i2)] = [x[V1], x[V2], x[V3], M]
            dforcesme0.append([x[V1], x[V2], x[V3], M])

            window_length = min(17, len(y) if len(y) % 2 != 0 else len(y) - 1)
            window_length = max(window_length, 5)
            if window_length % 2 == 0:
                window_length += 1

            y_smoothed = savgol_filter(y, window_length=window_length, polyorder=3)

            color = 'tab:blue'
            st = 1
            en = min(V3 + 40, len(x))
            ax.plot(x[st:en], y_smoothed[st:en], color=color, alpha=0.8)
            ax.fill_between(x[st:en], y_smoothed[st:en], color=color, alpha=0.3)
    else:
        print(f"Target {target} not found in label_data or file does not exist: {file_path}")

    textstyle = {'family': 'Times New Roman', 'size': 28}
    tick_fontsize = 20
    tick_font = 'Times New Roman'
    ax.set_xlabel('t (s)', fontdict=textstyle)
    ax.set_ylabel('Driving Forces (N)', fontdict=textstyle)
    plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.xticks(fontsize=tick_fontsize, fontfamily=tick_font)
    output_folder = SCRIPT_DIR / "driving_force_outputs"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_path = output_folder / f'drives_{target[0]}_{target[1]}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    saved_plot_path = output_path
    return dforcesme, dforcesme0, saved_plot_path


def save_smoothed_pressure_to_excel(target, label_data):
    """类似绘图函数，但将平滑后的数据保存为 Excel 文件"""
    i1, i2 = target
    file_path = os.path.join(
        "C:/Vault/10.Project/A009.BYD_R2/20220710/DATA_PROC/CUT/data_allclose_22919",
        f"MR_door_{i1}{i2}.asc"
    )

    if (i1, i2) not in label_data or not os.path.isfile(file_path):
        print(f"Target {target} not found in label_data or file does not exist: {file_path}")
        return None

    test = post_test(file_path)
    if test is None:
        print(f"File {file_path} could not be loaded.")
        return None

    x = test['t'].to_numpy()
    y = test['canpress'].to_numpy()
    aaa = y > 10
    V3 = len(y) - np.argmax(aaa[::-1]) - 1 if np.any(aaa) else None

    if V3 is None:
        print(f"Could not determine end index for target {target}.")
        return None

    window_length = min(17, len(y) if len(y) % 2 != 0 else len(y) - 1)
    window_length = max(window_length, 5)
    if window_length % 2 == 0:
        window_length += 1

    y_smoothed = savgol_filter(y, window_length=window_length, polyorder=3)
    st = 1
    en = min(V3 + 40, len(x))
    export_df = pd.DataFrame({
        "t": x[st:en],
        "canpress": y[st:en],
        "canpress_smoothed": y_smoothed[st:en],
    })

    output_folder = SCRIPT_DIR / "driving_force_data"
    output_folder.mkdir(parents=True, exist_ok=True)
    base_name = f'drives_{target[0]}_{target[1]}'
    output_path = output_folder / f"{base_name}.xlsx"

    try:
        export_df.to_excel(output_path, index=False)
    except PermissionError:
        # Fall back to a timestamped name if the base file is locked
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_path = output_folder / f"{base_name}_{timestamp}.xlsx"
        export_df.to_excel(output_path, index=False)
        print(f"Base Excel locked; saved to {output_path.resolve()} instead (generated by {Path(__file__).name})")
    else:
        print(f"Data saved to {output_path.resolve()} (generated by {Path(__file__).name})")

    return output_path


def _format_float_list(values):
    """Format an iterable of numbers as a comma-separated string for ADM output."""
    return ', '.join(f"{v:.6g}" for v in values)


def write_spline_adm_from_excel(excel_file, spline_name, adams_id, model_prefix=".MODEL_CLOSE_ST2", comments=""):
    """
    Read the first two columns (x, y) from an Excel file and write an Adams spline definition.
    Returns the output file path, or None on failure.
    """
    try:
        df = pd.read_excel(excel_file)
    except Exception as exc:
        print(f"Failed to read Excel file {excel_file}: {exc}")
        return None

    if df.shape[1] < 2:
        print(f"Excel file {excel_file} must have at least two columns (x, y).")
        return None

    x_vals = df.iloc[:, 0].to_numpy()
    y_vals = df.iloc[:, 1].to_numpy()

    output_folder = SCRIPT_DIR / "driving_force_splines"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_path = output_folder / f"{spline_name}.adm"

    content = (
        f"      data_element create spline &\n"
        f"         spline={model_prefix}.{spline_name} &\n"
        f"         adams_id={adams_id} &\n"
        f"         x={_format_float_list(x_vals)} &\n"
        f"         y={_format_float_list(y_vals)} &\n"
        f"         linear_extrapolate=no &\n"
        f"         x_units=no_units &\n"
        f"         y_units=no_units &\n"
        f"         comments=\"{comments}\"\n"
    )

    output_path.write_text(content, encoding="utf-8")
    print(f"ADM spline written to {output_path.resolve()} (generated by {Path(__file__).name})")
    return output_path


def main():
    """主函数"""
    # 标签数据
    label_data = {
        (24, 1): 0.86, (24, 2): 0.80, (24, 3): 0.82, (24, 4): 1.00, (24, 5): 1.04, (24, 6): 1.20, (24, 7): 1.21,
        (25, 1): 1.21, (25, 2): 1.17, (25, 3): 1.06, (25, 4): 0.94, (25, 5): 0.89, (25, 6): 0.68, (25, 7): 0.95, (25, 8): 0.73,
        (26, 1): 0.80, (26, 2): 0.87, (26, 3): 0.87, (26, 4): 0.92, (26, 5): 0.99, (26, 6): 1.17, (26, 7): 0.80,
        (27, 1): 0.69, (27, 2): 0.74, (27, 3): 0.92, (27, 4): 1.05, (27, 5): 1.00,
        (28, 1): 1.16, (28, 2): 0.82, (28, 3): 0.90, (28, 4): 1.04,
        (29, 1): 0.56, (29, 2): 0.67, (29, 3): 0.80, (29, 4): 0.78, (29, 5): 0.77, (29, 6): 0.95, (29, 7): 0.95, (29, 8): 1.03,
        (30, 1): 0.95, (30, 2): 1.14, (30, 3): 0.73, (30, 4): 0.92, (30, 5): 0.97,
    }

    # 选择要处理的目标 (i1, i2)
    target = (29, 4)

    dforcesme, dforcesme0, saved_plot_path = plot_smoothed_pressure(target, label_data)
    saved_excel_path = save_smoothed_pressure_to_excel(target, label_data)
    saved_adm_path = None
    if saved_excel_path:
        spline_name = f"SPLINE_{target[0]}_{target[1]}"
        adams_id = int(f"{target[0]}{target[1]}")
        saved_adm_path = write_spline_adm_from_excel(
            saved_excel_path,
            spline_name=spline_name,
            adams_id=adams_id,
            model_prefix=".MODEL_CLOSE_ST2",
            comments=""
        )

    if not dforcesme:
        print("Warning: No data files were found or processed for the target. Halting script.")
        return

    if saved_plot_path:
        print(f"Plot saved to {saved_plot_path.resolve()} (generated by {Path(__file__).name})")
    if saved_excel_path:
        print(f"Excel saved to {saved_excel_path.resolve()} (generated by {Path(__file__).name})")
    if saved_adm_path:
        print(f"ADM saved to {saved_adm_path.resolve()} (generated by {Path(__file__).name})")


if __name__ == "__main__":
    main()
