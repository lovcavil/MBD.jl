#!/usr/bin/env python3
"""Compare experiment vs simulation performance with interactive folder selection."""

from __future__ import annotations

import json
from pathlib import Path
from typing import List, Tuple, Dict, Any

missing_deps: list[str] = []

try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    plt = None  # type: ignore[assignment]
    missing_deps.append("matplotlib")

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
DATA_DIR = SCRIPT_DIR / "DATA"

# Default parameters (used if not in config)
DEFAULT_SIM_FILE = "req_KALKER_MID_F_clean.csv"
DEFAULT_SELECTED_COLUMNS_SIM = ["U3", "U2"]
DEFAULT_START_TIME_SIM = 0
DEFAULT_END_TIME_SIM = 2
DEFAULT_DATA_AMOUNT = 8
DEFAULT_SMOOTHING_METHOD = "Rolling Mean"
DEFAULT_WINDOW_SIZES_EXPR = [1] * 9
DEFAULT_WINDOW_SIZES_SIM = [1]

SIM_TIME_COL = "TIME"


def _check_dependencies() -> bool:
    """Check if required dependencies are installed."""
    if missing_deps:
        print("Missing dependencies:", ", ".join(sorted(set(missing_deps))))
        print("Install them before running this script.")
        return False
    return True


def _load_config(sim_folder: Path) -> Dict[str, Any] | None:
    """Load A20_config.json from sim folder."""
    config_path = sim_folder / "A20_config.json"
    if not config_path.exists():
        return None

    try:
        with config_path.open("r") as f:
            config = json.load(f)
        return config
    except json.JSONDecodeError:
        print(f"  Warning: {sim_folder.name} has invalid JSON format")
        return None


def _discover_configured_sim_folders() -> List[Tuple[Path, Dict[str, Any]]]:
    """Find all sim* folders with valid A20_config.json."""
    configured = []

    for sim_folder in sorted(DATA_DIR.glob("sim*")):
        if not sim_folder.is_dir():
            continue

        config = _load_config(sim_folder)
        if config is not None:
            configured.append((sim_folder, config))
        else:
            print(f"  Skipping {sim_folder.name}: no A20_config.json found")

    return configured


def _display_menu_and_select(folders: List[Tuple[Path, Dict[str, Any]]]) -> Tuple[Path, Dict[str, Any]]:
    """Display numbered menu and get user selection."""
    print("\nAvailable sim folders:")
    for i, (folder, _) in enumerate(folders, 1):
        print(f"  [{i}] {folder.name}")

    while True:
        try:
            choice = input(f"\nSelect folder (1-{len(folders)}): ").strip()
            index = int(choice) - 1
            if 0 <= index < len(folders):
                return folders[index]
            print(f"Please enter a number between 1 and {len(folders)}")
        except ValueError:
            print("Please enter a valid number")
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            raise SystemExit(0)


def load_data(path: Path):
    """Load CSV data."""
    return pd.read_csv(path, low_memory=False)


def filter_data(df, time_col: str, start_time: float, end_time: float, columns: List[str], max_columns: int):
    """Filter data by time range and limit columns."""
    time_series = pd.to_numeric(df[time_col], errors="coerce")
    mask = (time_series >= start_time) & (time_series <= end_time)
    filtered = df.loc[mask].copy()
    filtered[time_col] = time_series[mask]
    return filtered, columns[:max_columns]


def calculate_magnitude(filtered_df, columns: List[str], magnitude_column: str = "magnitude"):
    """Calculate sqrt(U3^2 + U2^2) for given columns and replace with single magnitude column."""
    if len(columns) >= 2 and "U3" in columns and "U2" in columns:
        # Convert to numeric
        u3 = pd.to_numeric(filtered_df["U3"], errors="coerce").values
        u2 = pd.to_numeric(filtered_df["U2"], errors="coerce").values
        # Calculate magnitude
        magnitude = np.sqrt(u3**2 + u2**2)
        # Add magnitude column
        filtered_df.loc[:, magnitude_column] = magnitude
        return [magnitude_column]
    return columns


def apply_scaling(filtered_df, columns: List[str], scaling_dict: Dict[str, float]):
    """Apply scaling factors to columns."""
    for column in columns:
        filtered_df.loc[:, column] = pd.to_numeric(filtered_df[column], errors="coerce")
        scale = scaling_dict.get(column, 1)
        filtered_df.loc[:, column] = filtered_df.loc[:, column] / scale
    return filtered_df


def apply_smoothing(filtered_df, columns: List[str], method: str, window_sizes: List[int]):
    """Apply smoothing to data."""
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
    """Create matplotlib figure and axes."""
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.set_facecolor("white")
    return fig, ax


def plot_experiment(ax, df, columns: List[str], start_time: float, label: str, color: str = "blue", linestyle: str = "-"):
    """Plot experiment data."""
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


def plot_simulation(ax, df, columns: List[str], label: str, color: str = "orange", linestyle: str = "--"):
    """Plot simulation data."""
    textstyle = {"family": "Times New Roman", "size": 32}
    tick_font = "Times New Roman"
    tick_fontsize = 20
    for column in columns:
        ax.plot(df["TIME"], df[column], label=label, color=color, linestyle=linestyle)
    ax.legend(fontsize=20, loc="upper left")
    plt.xticks(fontsize=tick_fontsize, fontfamily=tick_font)
    plt.yticks(fontsize=tick_fontsize, fontfamily=tick_font)
    return ax


def main():
    """Main execution function."""
    if not _check_dependencies():
        return
    if plt is None:
        print("matplotlib is not installed; cannot generate the plot.")
        return

    # Discover configured sim folders
    print("Scanning for configured sim folders...")
    configured_folders = _discover_configured_sim_folders()

    if not configured_folders:
        print("\nNo sim folders with A20_config.json found.")
        print("Create A20_config.json in sim folders to use this script.")
        return

    # Display menu and get selection
    sim_folder, config = _display_menu_and_select(configured_folders)

    # Extract parameters from config
    scaling_dict = config.get("scaling_dict", {})
    expr_file = config.get("expr_file", "MR_door (run 30)_out.csv")
    selected_columns_expr = config.get("selected_columns_expr", ["Strain@B3_2.RN_7"])
    start_time_expr = config.get("start_time_expr", 62.1)

    # Use default parameters for simulation
    sim_file = config.get("sim_file", DEFAULT_SIM_FILE)
    selected_columns_sim = config.get("selected_columns_sim", DEFAULT_SELECTED_COLUMNS_SIM)

    # Build paths
    path_expr = EXPR_DIR / expr_file
    path_sim = sim_folder / sim_file

    # Load data
    print(f"\nLoading data from:")
    print(f"  Experiment: {path_expr}")
    print(f"  Simulation: {path_sim}")

    df_expr = load_data(path_expr)
    df_sim = load_data(path_sim)

    # Time parameters
    end_time_expr = start_time_expr + 2.0
    start_time_sim = DEFAULT_START_TIME_SIM
    end_time_sim = DEFAULT_END_TIME_SIM

    # Filter data
    filtered_df_expr, selected_columns_expr = filter_data(
        df_expr, "t", start_time_expr, end_time_expr, selected_columns_expr, DEFAULT_DATA_AMOUNT
    )
    filtered_df_sim, selected_columns_sim = filter_data(
        df_sim, SIM_TIME_COL, start_time_sim, end_time_sim, selected_columns_sim, DEFAULT_DATA_AMOUNT
    )

    # Calculate magnitude for simulation data (sqrt(U3^2 + U2^2))
    selected_columns_sim = calculate_magnitude(filtered_df_sim, selected_columns_sim, magnitude_column="magnitude")

    # Apply scaling
    filtered_df_expr = apply_scaling(filtered_df_expr, selected_columns_expr, scaling_dict)

    # Apply smoothing
    filtered_df_expr = apply_smoothing(filtered_df_expr, selected_columns_expr, DEFAULT_SMOOTHING_METHOD, DEFAULT_WINDOW_SIZES_EXPR)
    filtered_df_sim = apply_smoothing(filtered_df_sim, selected_columns_sim, DEFAULT_SMOOTHING_METHOD, DEFAULT_WINDOW_SIZES_SIM)

    # Create plot
    fig, ax = create_axes()
    ax = plot_experiment(ax, filtered_df_expr, selected_columns_expr, start_time_expr, "Experiment", color="r", linestyle="--")
    ax = plot_simulation(
        ax,
        filtered_df_sim,
        selected_columns_sim,
        "CONTACT_MID_F_KALKER / Y",
        color="k",
        linestyle="-",
    )

    # Generate output filename based on folder name
    folder_name = sim_folder.name
    # Remove "sim" prefix if present
    if folder_name.startswith("sim"):
        folder_id = folder_name[3:]  # Remove "sim" prefix
    else:
        folder_id = folder_name

    output_folder = SCRIPT_DIR / "validation_outputs"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_path = output_folder / f"A20.compare_expr_sim_performance_{folder_id}.png"

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"\n✓ Plot saved to {output_path.resolve()}")


if __name__ == "__main__":
    main()
