#!/usr/bin/env python3

import argparse
import os
import subprocess


def main():
    parser = argparse.ArgumentParser(description="Run log_to_excel and make_charts pipeline")
    parser.add_argument(
        "-b", "--base_dir",
        required=True,
        help="Base directory (e.g. update_tests/log13_250slot_rms_evm_normalization)"
    )
    parser.add_argument(
        "-c", "--charts_config",
        default="charts.yml",
        help="Path to charts.yml (default: charts.yml)"
    )

    args = parser.parse_args()

    base = os.path.abspath(args.base_dir)

    # ---- Construct paths ----
    log_dir = os.path.join(base, "log")
    excel_dir = base
    excel_file = os.path.join(base, "runlog_results.xlsx")
    img_dir = os.path.join(base, "img")

    print("\n=== Running log_to_excel_v2 ===")
    subprocess.run([
        "python3", "-m", "log_to_excel_v2",
        # "python3", "-m", "log_to_excel_v3",
        "-d", log_dir,
        "-ex", excel_dir
    ], check=True)

    print("\n=== Running make_charts ===")
    subprocess.run([
        "python3", "-m", "make_charts",
        "-e", excel_file,
        "-c", args.charts_config,
        "-o", img_dir
    ], check=True)

    print("\n✅ Pipeline complete")
    print(f"Excel: {excel_file}")
    print(f"Charts: {img_dir}")


if __name__ == "__main__":
    main()