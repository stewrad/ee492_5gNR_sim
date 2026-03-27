#!/usr/bin/env python3
"""
Parse 5G NR HARQ runlog_*.txt files and append metrics to an Excel file.

Requires:
  pip install pandas openpyxl

Usage examples:
  # Logs in current folder, excel_path is a directory -> creates ./runlog_results.xlsx
  python log_to_excel.py --log_dir . --excel_path .

  # Logs in a folder, output file named explicitly
  python log_to_excel.py --log_dir tests/AustinS/logs --excel_path results.xlsx

  # Avoid duplicates based on (Source File, Timestamp)
  python log_to_excel.py --log_dir . --excel_path . --dedupe
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from typing import Dict, Any, Optional, List

import pandas as pd


# ---------- Regex helpers ----------

def _re_float(pattern: str, text: str, flags=0) -> Optional[float]:
    m = re.search(pattern, text, flags)
    return float(m.group(1)) if m else None


def _re_int(pattern: str, text: str, flags=0) -> Optional[int]:
    m = re.search(pattern, text, flags)
    return int(m.group(1)) if m else None


def _re_str(pattern: str, text: str, flags=0) -> Optional[str]:
    m = re.search(pattern, text, flags)
    return m.group(1).strip() if m else None


def _parse_sample_rate_to_hz(sample_rate_str: str) -> Optional[float]:
    """
    Accepts strings like:
      "15.36 MHz"
      "30.72 Msps"
      "15360 kHz"
    Returns Hz as float.
    """
    if not sample_rate_str:
        return None
    s = sample_rate_str.strip()

    m = re.match(r"^\s*([0-9]*\.?[0-9]+)\s*([a-zA-Z/]+)\s*$", s)
    if not m:
        return None

    val = float(m.group(1))
    unit = m.group(2).lower()

    if unit in ("hz",):
        mult = 1.0
    elif unit in ("khz",):
        mult = 1e3
    elif unit in ("mhz",):
        mult = 1e6
    elif unit in ("ghz",):
        mult = 1e9
    elif unit in ("sps",):
        mult = 1.0
    elif unit in ("ksps",):
        mult = 1e3
    elif unit in ("msps",):
        mult = 1e6
    elif unit in ("gsps",):
        mult = 1e9
    else:
        return None

    return val * mult


# ---------- Main parser ----------

BASE_COLUMNS = [
    "SNRdB",
    "Modulation",
    "HARQ Type",
    "NHARQProcesses",
    "rvSeq",
    "nTxAnts",
    "nRxAnts",
    "DelayProfile",
    "NumLayers",  # optional in logs
    "Code Rate",
    "Total Transmissions",
    "Initial Transmissions",
    "Retransmissions",
    "Average Transmissions per TB",
    "Retransmission Rate (%)",
    "Successful Blocks",
    "Failed Blocks (after max retx)",
    "Attempt BLER % (pre-HARQ)",
    "Final BLER % (post-HARQ)",
    "Block Error Rate (BLER %)",
    "Total Bits Transmitted",
    "Total Bits Received",
    "Throughput Efficiency (%)",
    "Effective Code Rate",
    "Average Throughput (Mbps)",
    "Average Bits per Slot (bits)",
    "Spectral Efficiency (bits/s/Hz)",
    "RMS Delay Spread (ns)",
    "Mean Excess Delay (ns)",
    "Maximum Channel Delay (us)",      # stored in microseconds (us)
    "Average Channel Gain (dB)",
    "Channel Gain Std Dev (dB)",
    "Average Condition Number",
    # Layer SINR columns added dynamically: "Layer 1 - Mean SINR", ...
    "Coherence Bandwidth (50% corr) (kHz)",        # kHz numeric
    "Coherence Time (50% corr) (ms)",             # ms numeric
]

TRACE_COLUMNS = ["Source File", "Timestamp"]


def parse_runlog(text: str, source_file: str) -> Dict[str, Any]:
    out: Dict[str, Any] = {}

    # --- Header / run parameters ---
    out["Timestamp"] = _re_str(r"Timestamp:\s*(.+)", text)

    out["SNRdB"] = _re_float(r"SNRdB:\s*([0-9]*\.?[0-9]+)", text)
    out["Modulation"] = (
        _re_str(r"Modulation:\s*\"([^\"]+)\"", text)
        or _re_str(r"Modulation:\s*([A-Za-z0-9]+)", text)
    )
    out["NHARQProcesses"] = _re_int(r"NHARQProcesses:\s*([0-9]+)", text)
    out["rvSeq"] = _re_str(r"rvSeq:\s*(\[[^\]]+\])", text)
    out["nTxAnts"] = _re_int(r"nTxAnts:\s*([0-9]+)", text)
    out["nRxAnts"] = _re_int(r"nRxAnts:\s*([0-9]+)", text)
    out["DelayProfile"] = (
        _re_str(r"DelayProfile:\s*\"([^\"]+)\"", text)
        or _re_str(r"Delay Profile:\s*([A-Za-z0-9\-]+)", text)
    )
    out["NumLayers"] = (
        _re_int(r"NumLayers:\s*([0-9]+)", text)
        or _re_int(r"Number of Layers:\s*([0-9]+)", text)
    )

    # --- Results section ---
    out["HARQ Type"] = _re_str(r"HARQ Type:\s*(.+)", text)

    out["Code Rate"] = _re_float(r"Code Rate:\s*([0-9]*\.?[0-9]+)", text)
    out["Total Transmissions"] = _re_int(r"Total Transmissions:\s*([0-9]+)", text)
    out["Initial Transmissions"] = _re_int(r"Initial Transmissions:\s*([0-9]+)", text)
    out["Retransmissions"] = _re_int(r"Retransmissions:\s*([0-9]+)", text)
    out["Average Transmissions per TB"] = _re_float(
        r"Average Transmissions per TB:\s*([0-9]*\.?[0-9]+)", text
    )

    rr = _re_float(r"Retransmission Rate:\s*([0-9]*\.?[0-9]+)\s*%", text)
    out["Retransmission Rate (%)"] = rr

    out["Successful Blocks"] = _re_int(r"Successful Blocks:\s*([0-9]+)\s*/\s*[0-9]+", text)
    out["Failed Blocks (after max retx)"] = _re_int(r"Failed Blocks \(after max retx\):\s*([0-9]+)", text)
    out["Block Error Rate (BLER %)"] = _re_float(r"Block Error Rate \(BLER\):\s*([0-9]*\.?[0-9]+)", text)

    # Attempt BLER (pre-HARQ) and Final BLER (post-HARQ) from summary section
    out["Attempt BLER % (pre-HARQ)"] = _re_float(
        r"Attempt BLER \(pre-HARQ\):\s*[0-9]*\.?[0-9]+\s*\(([0-9]*\.?[0-9]+)%\)", text
    )
    out["Final BLER % (post-HARQ)"] = _re_float(
        r"Final BLER \(post-HARQ\):\s*[0-9]*\.?[0-9]+\s*\(([0-9]*\.?[0-9]+)%\)", text
    )

    out["Total Bits Transmitted"] = _re_int(r"Total Bits Transmitted:\s*([0-9]+)", text)
    out["Total Bits Received"] = _re_int(r"Total Bits Received \(success\):\s*([0-9]+)", text)

    out["Throughput Efficiency (%)"] = _re_float(r"Throughput Efficiency:\s*([0-9]*\.?[0-9]+)\s*%", text)
    out["Effective Code Rate"] = _re_float(r"Effective Code Rate:\s*([0-9]*\.?[0-9]+)", text)

    out["Average Throughput (Mbps)"] = _re_float(r"Average Throughput:\s*([0-9]*\.?[0-9]+)\s*Mbps", text)
    out["Average Bits per Slot (bits)"] = _re_int(r"Average Bits per Slot:\s*([0-9]+)\s*bits", text)

    out["Spectral Efficiency (bits/s/Hz)"] = _re_float(r"Spectral Efficiency:\s*([0-9]*\.?[0-9]+)\s*bits/s/Hz", text)

    # --- Channel details ---
    out["RMS Delay Spread (ns)"] = _re_float(r"RMS Delay Spread:\s*([0-9]*\.?[0-9]+)\s*ns", text)
    out["Mean Excess Delay (ns)"] = _re_float(r"Mean Excess Delay:\s*([0-9]*\.?[0-9]+)\s*ns", text)

    max_delay_samples = _re_int(r"Maximum Channel Delay:\s*([0-9]+)\s*samples", text)
    max_delay_us = _re_float(
        r"Maximum Channel Delay:\s*[0-9]+\s*samples\s*\(\s*([0-9]*\.?[0-9]+)\s*μs\s*\)", text
    )

    if max_delay_us is None and max_delay_samples is not None:
        sr_str = _re_str(r"Sample Rate:\s*([0-9]*\.?[0-9]+\s*[A-Za-z]+)", text)
        sr_hz = _parse_sample_rate_to_hz(sr_str) if sr_str else None
        if sr_hz:
            max_delay_us = (max_delay_samples / sr_hz) * 1e6

    out["Maximum Channel Delay (us)"] = max_delay_us

    out["Average Channel Gain (dB)"] = _re_float(r"Average Channel Gain:\s*([\-0-9]*\.?[0-9]+)\s*dB", text)
    out["Channel Gain Std Dev (dB)"] = _re_float(r"Channel Gain Std Dev:\s*([0-9]*\.?[0-9]+)\s*dB", text)
    out["Average Condition Number"] = _re_float(r"Average Condition Number:\s*([0-9]*\.?[0-9]+)", text)

    out["Coherence Bandwidth (50% corr) (kHz)"] = _re_float(r"Coherence Bandwidth.*?:\s*([0-9]*\.?[0-9]+)\s*kHz", text)
    out["Coherence Time (50% corr) (ms)"] = _re_float(r"Coherence Time.*?:\s*([0-9]*\.?[0-9]+)\s*ms", text)

    # --- Layer SINR: dynamic columns ---
    for m in re.finditer(r"Layer\s+([0-9]+)\s*-\s*Mean SINR:\s*([\-0-9]*\.?[0-9]+)\s*dB", text):
        layer_idx = int(m.group(1))
        mean_sinr = float(m.group(2))
        out[f"Layer {layer_idx} - Mean SINR"] = mean_sinr

    # --- Per-layer EVM metrics: average across all slots ---
    # Log line format (one per layer per slot):
    #   Slot  N | Layer L | RMS EVM = X% | Peak EVM = Y% | Avg EVM = Z dB | Peak EVM = W dB | Avg MER = V dB
    evm_pattern = re.compile(
        r"Slot\s+\d+\s*\|\s*Layer\s+(\d+)\s*\|"
        r"\s*RMS EVM\s*=\s*([\-0-9]*\.?[0-9]+)%"
        r"\s*\|\s*Peak EVM\s*=\s*([\-0-9]*\.?[0-9]+)%"
        r"\s*\|\s*Avg EVM\s*=\s*([\-0-9]*\.?[0-9]+)\s*dB"
        r"\s*\|\s*Peak EVM\s*=\s*([\-0-9]*\.?[0-9]+)\s*dB"
        r"\s*\|\s*Avg MER\s*=\s*([\-0-9]*\.?[0-9]+)\s*dB"
    )

    # Accumulate per-layer lists
    evm_accum: Dict[int, Dict[str, List[float]]] = {}
    for m in evm_pattern.finditer(text):
        lyr = int(m.group(1))
        if lyr not in evm_accum:
            evm_accum[lyr] = {
                "rms_pct": [], "peak_pct": [],
                "avg_db": [], "peak_db": [], "mer_db": []
            }
        evm_accum[lyr]["rms_pct"].append(float(m.group(2)))
        evm_accum[lyr]["peak_pct"].append(float(m.group(3)))
        evm_accum[lyr]["avg_db"].append(float(m.group(4)))
        evm_accum[lyr]["peak_db"].append(float(m.group(5)))
        evm_accum[lyr]["mer_db"].append(float(m.group(6)))

    import statistics
    for lyr, vals in evm_accum.items():
        prefix = f"L{lyr}"
        out[f"{prefix} RMS EVM (%)"]   = statistics.mean(vals["rms_pct"])
        out[f"{prefix} Peak EVM (%)"]  = statistics.mean(vals["peak_pct"])
        out[f"{prefix} Avg EVM (dB)"]  = statistics.mean(vals["avg_db"])
        out[f"{prefix} Peak EVM (dB)"] = statistics.mean(vals["peak_db"])
        out[f"{prefix} Avg MER (dB)"]  = statistics.mean(vals["mer_db"])

    out["Source File"] = os.path.basename(source_file)
    return out


def ensure_columns(df: pd.DataFrame, columns: List[str]) -> pd.DataFrame:
    for c in columns:
        if c not in df.columns:
            df[c] = pd.NA
    return df


def normalize_excel_path(excel_path: str, excel_name: str) -> str:
    """
    If excel_path is a directory (e.g. "."), create a default filename inside it.
    Also ensure .xlsx extension.
    """
    p = os.path.expanduser(excel_path)

    # If they passed a directory, create default filename inside it
    if os.path.isdir(p):
        p = os.path.join(p, excel_name)

    # If they passed something without extension, add .xlsx
    if not p.lower().endswith(".xlsx"):
        p = p + ".xlsx"

    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--log_dir", required=True, help="Directory containing runlog_*.txt files")
    ap.add_argument("-ex", "--excel_path", required=True, help="Path to output .xlsx OR a directory like '.'")
    ap.add_argument("-fname", "--excel_name", default="runlog_results.xlsx", help="Excel file name")
    ap.add_argument("-sn", "--sheet", default="Results", help="Worksheet name")
    ap.add_argument("--pattern", default="*.txt", help="Glob pattern inside log_dir")
    ap.add_argument(
        "--dedupe",
        action="store_true",
        help="If set, avoid appending rows whose (Source File, Timestamp) already exist in Excel.",
    )
    args = ap.parse_args()

    # Normalize paths so "." works reliably
    args.log_dir = os.path.abspath(os.path.expanduser(args.log_dir))
    args.excel_path = normalize_excel_path(args.excel_path, args.excel_name)

    files = sorted(glob.glob(os.path.join(args.log_dir, args.pattern)))
    if not files:
        raise SystemExit(f"No files matched '{args.pattern}' in {args.log_dir}")

    rows: List[Dict[str, Any]] = []
    for fp in files:
        with open(fp, "r", encoding="utf-8", errors="replace") as f:
            txt = f.read()
        rows.append(parse_runlog(txt, fp))

    new_df = pd.DataFrame(rows)

    # Determine all columns (base + trace + any Layer SINR columns found)
    layer_sinr_cols = sorted(
        [c for c in new_df.columns if c.startswith("Layer ") and c.endswith(" - Mean SINR")],
        key=lambda x: int(re.search(r"Layer\s+([0-9]+)", x).group(1)),
    )

    # Per-layer EVM columns: group by layer number, then metric order within each layer
    evm_metric_order = ["RMS EVM (%)", "Peak EVM (%)", "Avg EVM (dB)", "Peak EVM (dB)", "Avg MER (dB)"]
    layer_nums = sorted(set(
        int(re.search(r"L(\d+)", c).group(1))
        for c in new_df.columns
        if re.match(r"L\d+ (RMS|Peak|Avg)", c)
    ))
    layer_evm_cols = [
        f"L{lyr} {metric}"
        for lyr in layer_nums
        for metric in evm_metric_order
        if f"L{lyr} {metric}" in new_df.columns
    ]

    all_cols = TRACE_COLUMNS + BASE_COLUMNS + layer_sinr_cols + layer_evm_cols

    new_df = ensure_columns(new_df, all_cols)
    new_df = new_df[all_cols]  # order columns

    # Append to existing Excel or create new
    if os.path.isfile(args.excel_path):
        existing = pd.read_excel(args.excel_path, sheet_name=args.sheet, engine="openpyxl")
        existing = ensure_columns(existing, all_cols)

        if args.dedupe and {"Source File", "Timestamp"}.issubset(existing.columns):
            key_existing = set(zip(existing["Source File"].astype(str), existing["Timestamp"].astype(str)))
            mask_keep = []
            for _, r in new_df.iterrows():
                k = (str(r["Source File"]), str(r["Timestamp"]))
                mask_keep.append(k not in key_existing)
            new_df = new_df[pd.Series(mask_keep).values]

        combined = pd.concat([existing, new_df], ignore_index=True)
    else:
        combined = new_df

    # Ensure output directory exists (in case they pass something like ./out/results.xlsx)
    out_dir = os.path.dirname(os.path.abspath(args.excel_path))
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Write back (replace sheet)
    with pd.ExcelWriter(args.excel_path, engine="openpyxl", mode="w") as xw:
        combined.to_excel(xw, sheet_name=args.sheet, index=False)

    print(
        f"Appended {len(new_df)} rows from {len(files)} files into: {args.excel_path} "
        f"(sheet '{args.sheet}')"
    )


if __name__ == "__main__":
    main()