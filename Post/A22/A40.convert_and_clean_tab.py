#!/usr/bin/env python3
"""Convert Adams *.tab exports directly into cleaned CSVs under DATA/sim."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, List

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "DATA"


def _iter_tab_files(base_dir: Path) -> Iterable[Path]:
    yield from sorted(base_dir.glob("*.tab"))


def _build_header_from_row(row: List[str]) -> List[str]:
    sanitized: List[str] = []
    for cell in row:
        trimmed = cell.strip()
        sanitized.append(trimmed.split(".")[-1] if trimmed else "")
    return sanitized


def _convert_single(tab_path: Path, dest_dir: Path) -> Path:
    dest_path = dest_dir / f"{tab_path.stem}_clean.csv"
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    header: List[str] | None = None
    data_rows: List[List[str]] = []

    with tab_path.open("r", encoding="utf-8", errors="replace", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        for row in reader:
            if not any(cell.strip() for cell in row):
                continue
            cells = [cell.strip() for cell in row]
            first_cell = cells[0].lstrip("\ufeff") if cells else ""
            if header is None and first_cell.startswith("."):
                header = _build_header_from_row(cells)
                continue
            if header is not None:
                data_rows.append(cells)

    if header is None:
        raise ValueError(f"No header line found in {tab_path.name}")

    with dest_path.open("w", newline="", encoding="utf-8") as writer_file:
        writer = csv.writer(writer_file)
        writer.writerow(header)
        writer.writerows(data_rows)
    return dest_path


def _discover_adams_folders(base_dir: Path) -> List[Path]:
    """Find all Adams* folders in the base directory."""
    return sorted([f for f in base_dir.glob("Adams*") if f.is_dir()])


def _get_sim_folder(adams_folder: Path, base_dir: Path) -> Path:
    """Generate sim folder name from Adams folder name."""
    folder_name = adams_folder.name
    sim_name = folder_name.replace("Adams", "sim")
    return base_dir / sim_name


def main() -> None:
    """Convert all .tab files from all Adams* folders to corresponding sim folders."""
    adams_folders = _discover_adams_folders(DATA_DIR)

    if not adams_folders:
        print("No Adams* folders found in DATA directory.")
        return

    total_converted = 0

    for adams_folder in adams_folders:
        sim_folder = _get_sim_folder(adams_folder, DATA_DIR)

        # Create sim folder if it doesn't exist
        sim_folder.mkdir(parents=True, exist_ok=True)

        # Convert all .tab files in this folder
        converted = []
        for tab_file in _iter_tab_files(adams_folder):
            converted_path = _convert_single(tab_file, sim_folder)
            converted.append(converted_path)

        if converted:
            print(f"Processing {adams_folder.name} → {sim_folder.name}")
            for path in converted:
                print(f"  Created: {path.relative_to(DATA_DIR)}")
            print(f"  ✓ {len(converted)} file{'s' if len(converted) > 1 else ''} converted\n")
            total_converted += len(converted)

    if total_converted > 0:
        print(f"Total: {total_converted} files converted across {len(adams_folders)} folder{'s' if len(adams_folders) > 1 else ''}")
    else:
        print("No .tab files found in any Adams folders.")


if __name__ == "__main__":
    main()
