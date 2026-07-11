"""
Standalone utility: bulk-convert Agilent binary .ch files to CSV.

NOT part of the main pipeline -- MOCCA2 reads .D folders natively, so you only
need this if you want CSV exports for Excel or another tool. Reads every .D
folder in a directory with the `rainbow` library (in parallel), exports each
data channel to CSV, then re-encodes as UTF-16 for Excel compatibility.

Usage:
    python tools/agilent_to_csv.py /path/to/folder/with/.D/subfolders

The dataset path is a CLI argument -- unlike the original, nothing is hardcoded.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from pathlib import Path

import rainbow as rb


def convert_utf8_to_utf16(input_file: str, output_file: str) -> None:
    """Re-encode a UTF-8 CSV as UTF-16 (Excel-friendly)."""
    try:
        with open(input_file, "r", encoding="utf-8") as f:
            content = f.read()
        with open(output_file, "w", encoding="utf-16") as f:
            f.write(content)
        print(f"  OK  {output_file}")
    except FileNotFoundError:
        print(f"  ERR file not found: {input_file}")
    except Exception as exc:
        print(f"  ERR {exc}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Convert Agilent .ch files in .D folders to CSV.")
    parser.add_argument("dataset", type=Path, help="Directory containing .D folders")
    args = parser.parse_args(argv)

    dataset = str(args.dataset)
    dirpaths = [os.path.join(dataset, name) for name in os.listdir(dataset) if name.endswith(".D")]
    if not dirpaths:
        print(f"No .D folders found in {dataset}")
        return 1

    with mp.Pool() as pool:
        datadirs = pool.map(rb.read, dirpaths)

    for injection in datadirs:
        for datafile in injection.datafiles:
            out_filename = os.path.join(dataset, injection.name, datafile.name[:-2] + "CSV")
            injection.export_csv(datafile.name, out_filename, labels=None, delim=",")
            convert_utf8_to_utf16(out_filename, out_filename)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
