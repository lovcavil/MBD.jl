#!/usr/bin/env python3
"""Convert an Adams-generated CSV into a plain time-series CSV."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, List, Optional


def _next_nonempty_line(lines: Iterable[str]) -> Optional[str]:
    """Return the next non-blank line in the iterator (without newline)."""
    for raw_line in lines:
        line = raw_line.rstrip("\n").rstrip("\r")
        if line.strip():
            return line
    return None


def _build_header(titles: Optional[List[str]], column_count: int) -> List[str]:
    """Build the output header, dropping the index column."""
    if titles:
        header = [cell.strip() for cell in titles[1:]]
        if len(header) < column_count:
            header.extend("" for _ in range(column_count - len(header)))
        else:
            header = header[:column_count]
    else:
        header = ["" for _ in range(column_count)]

    if not header:
        return header

    if not header[0]:
        header[0] = "t"

    for idx, value in enumerate(header):
        if not value:
            header[idx] = f"col{idx + 1}"
    return header


def convert(input_path: Path, output_path: Path) -> None:
    """Read the Adams CSV, drop the step column, and write the clean CSV."""
    titles: Optional[List[str]] = None
    data_mode = False
    header_written = False

    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", newline="") as reader, output_path.open("w", newline="") as writer_file:
        writer = csv.writer(writer_file)
        lines = iter(reader)
        for raw_line in lines:
            line = raw_line.rstrip("\n").rstrip("\r")
            stripped = line.strip()
            if not stripped:
                continue

            if stripped.startswith("#"):
                if stripped.startswith("#TITLES"):
                    next_line = _next_nonempty_line(lines)
                    if next_line is None:
                        raise ValueError("Missing title row after #TITLES")
                    titles = next(csv.reader([next_line]))
                elif stripped.startswith("#DATA"):
                    data_mode = True
                continue

            if not data_mode:
                continue

            row = next(csv.reader([line]))
            if len(row) <= 1:
                continue

            if not header_written:
                header = _build_header(titles, len(row) - 1)
                writer.writerow(header)
                header_written = True

            writer.writerow(row[1:])

        if not header_written:
            raise ValueError("No data rows found after #DATA")




INPUT_PATH =  Path("Post/A22/DATA/csv/MR_door (run 24)_out.csv")
OUTPUT_PATH = Path("Post/A22/DATA/expr//MR_door (run 24)_out.csv")


def main() -> None:
    convert(INPUT_PATH, OUTPUT_PATH)


if __name__ == "__main__":
    main()
