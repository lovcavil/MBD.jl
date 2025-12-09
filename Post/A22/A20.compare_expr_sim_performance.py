from pathlib import Path

missing_deps: list[str] = []

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError as error:
    plt = None
    missing_deps.append("matplotlib")
else:
    _MATPLOTLIB_ERROR = None

try:
    import numpy as np
except ModuleNotFoundError:
    np = None  # type: ignore[assignment]
    missing_deps.append("numpy")

try:
    import pandas as pd
except ModuleNotFoundError:
    pd = None  # type: ignore[assignment]
    missing_deps.append("pandas")

try:
    from scipy.fft import fft, ifft
except ModuleNotFoundError:
    fft = None  # type: ignore[assignment]
    ifft = None  # type: ignore[assignment]
    missing_deps.append("scipy")

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
PLOTS_DIR = REPO_ROOT / "plots" / "adams2"
EXPR_DIR = SCRIPT_DIR / "DATA" / "expr"
SIM_DIR = SCRIPT_DIR / "DATA" / "sim"

scaling_dict = {
    "Strain@B1_1.RN_6": 0.339,
    "Strain@B2_1.RN_6": 0.3505,
    "Strain@B2_2.RN_6": 0.5762,
    "Strain@B2_3.RN_6": -0.5072,
    "Strain@B2_4.RN_6": 0.0763,
    "Strain@B3_1.RN_6": 0.8693,
    "Strain@B3_2.RN_6": 0.8576,
    "Strain@B3_3.RN_6": -0.2255,
    "Strain@B3_4.RN_6": 0.2455,
}


def load_data(path):
    return pd.read_csv(path, low_memory=False)


def filter_data(df, time_col, start_time, end_time, columns, max_columns):
    time_series = pd.to_numeric(df[time_col], errors="coerce")
    mask = (time_series >= start_time) & (time_series <= end_time)
    filtered = df.loc[mask].copy()
    filtered[time_col] = time_series[mask]
    return filtered, columns[:max_columns]


def apply_scaling(filtered_df, columns):
    for column in columns:
        filtered_df.loc[:, column] = pd.to_numeric(filtered_df[column], errors="coerce")
        scale = scaling_dict.get(column, 1)
        filtered_df.loc[:, column] = filtered_df.loc[:, column] / scale
    return filtered_df


def apply_smoothing(filtered_df, columns, method, window_sizes):
    if method == "Rolling Mean":
        for column, window_size in zip(columns, window_sizes):
            filtered_df.loc[:, column] = filtered_df.loc[:, column].rolling(window=window_size, min_periods=1).mean()
    elif method == "Fourier Filter":
        for column, window_size in zip(columns, window_sizes):
            signal = filtered_df[column].values
            fft_vals = fft(signal)
            fft_vals[window_size:] = 0
            filtered_df.loc[:, column] = np.real(ifft(fft_vals))
    return filtered_df


def create_axes():
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.set_facecolor("white")
    return fig, ax


def plot_experiment(ax, df, columns, start_time, label, color="blue", linestyle="-"):
    textstyle = {"family": "Times New Roman", "size": 32}
    tick_font = "Times New Roman"
    tick_fontsize = 20
    for column in columns:
        ax.plot(df["t"] - start_time, df[column], label=label, color=color, linestyle=linestyle)
    ax.set_xlabel("Time (s)", fontdict=textstyle)
    ax.set_ylabel("Normal Contact Force(N)", fontdict=textstyle)
    plt.xticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.tight_layout()
    return ax


def plot_simulation(ax, df, columns, label, color="orange", linestyle="--"):
    textstyle = {"family": "Times New Roman", "size": 32}
    tick_font = "Times New Roman"
    tick_fontsize = 20
    for column in columns:
        ax.plot(df["Time"], df[column], label=label, color=color, linestyle=linestyle)
    ax.legend(fontsize=20, loc="upper left")
    plt.xticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)
    return ax


SIM_TIME_COL = "Time"


def build_plot_paths():
    return (
        EXPR_DIR / "MR_door (run 24)_out.csv",
        SIM_DIR / "request_CONTACT_MID_F_KALKER_SS_clean.csv",
        PLOTS_DIR / "perfectA18" / "perfect.csv",
    )


def _check_dependencies():
    if missing_deps:
        print("Missing dependencies:", ", ".join(sorted(set(missing_deps))))
        print("Install them before running this script.")
        return False
    return True


def main():
    if not _check_dependencies():
        return
    if plt is None:
        print("matplotlib is not installed; cannot generate the plot.")
        return
    if plt is None:
        print("matplotlib is not installed in this environment; install it to run the comparison.")
        return
    path1, path2, path5 = build_plot_paths()

    df1 = load_data(path1)
    df2 = load_data(path2)
    df5 = load_data(path5)

    # Parameters
    selected_columns1 = ["Strain@B3_2.RN_1"]#"Strain@B3_2.RN_1" Strain@B3_2.RN_6
    selected_columns2 = ["U3"]
    selected_columns5 = ["FY"]
    start_time1 = 60.8
    end_time1 = start_time1 + 1.35
    start_time2 = 0
    end_time2 = 2
    start_time3 = 0
    end_time3 = 1.46
    data_amount = 8
    smoothing_method = "Rolling Mean"
    window_sizes1 = [1] * 9
    window_sizes2 = [1] * 9
    window_sizes5 = [2] * 9

    filtered_df1, selected_columns1 = filter_data(df1, "t", start_time1, end_time1, selected_columns1, data_amount)
    filtered_df2, selected_columns2 = filter_data(df2, SIM_TIME_COL, start_time2, end_time2, selected_columns2, data_amount)
    filtered_df5, selected_columns5 = filter_data(df5, "Time", start_time3, end_time3, selected_columns5, data_amount)

    # reflected_df = filtered_df2[(filtered_df2[SIM_TIME_COL] > 0.1) & (filtered_df2[SIM_TIME_COL] < 0.2)].copy()
    # reflected_df[SIM_TIME_COL] = 0.1 + (0.1 - reflected_df[SIM_TIME_COL])
    # filtered_df2 = filtered_df2[(filtered_df2[SIM_TIME_COL] > 0.1)]
    # filtered_df2 = pd.concat([filtered_df2, reflected_df]).sort_values(by=SIM_TIME_COL).reset_index(drop=True)

    filtered_df1 = apply_scaling(filtered_df1, selected_columns1)

    filtered_df1 = apply_smoothing(filtered_df1, selected_columns1, smoothing_method, window_sizes1)
    filtered_df2 = apply_smoothing(filtered_df2, selected_columns2, smoothing_method, window_sizes2)
    filtered_df5 = apply_smoothing(filtered_df5, selected_columns5, smoothing_method, window_sizes5)

    fig, ax = create_axes()
    ax = plot_experiment(ax, filtered_df1, selected_columns1, start_time1, "Experiment", color="r", linestyle="--")
    ax = plot_simulation(
        ax,
        filtered_df2,
        selected_columns2,
        "CONTACT_MID_F_KALKER / Y",
        color="k",
        linestyle="-",
    )

    output_folder = SCRIPT_DIR / "validation_outputs"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_path = output_folder / "A20.compare_expr_sim_performance.png"
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved to {output_path.resolve()} (generated by {Path(__file__).name})")


if __name__ == "__main__":
    main()
