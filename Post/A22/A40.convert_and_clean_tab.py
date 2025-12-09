#!/usr/bin/env python3
"""Convert Adams *.tab exports directly into cleaned CSVs under DATA/sim."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable, List

SCRIPT_DIR = Path(__file__).resolve().parent
ADAMS_DIR = SCRIPT_DIR / "DATA" / "Adams"
SIM_DIR = SCRIPT_DIR / "DATA" / "sim"


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


def main() -> None:
    args = parse_args()
    dest_dir = args.output_dir.resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)
    converted = []
    for tab_file in _iter_tab_files(args.source_dir):
        converted_path = _convert_single(tab_file, dest_dir)
        converted.append(converted_path)

    if converted:
        print("Created cleaned CSVs:")
        for path in converted:
            print(f" - {path}")
    else:
        print("No .tab files found.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert Adams tab exports into cleaned CSVs.")
    parser.add_argument(
        "--source-dir",
        type=Path,
        default=ADAMS_DIR,
        help="Directory containing Adams *.tab exports",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SIM_DIR,
        help="Destination directory for cleaned CSVs",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
