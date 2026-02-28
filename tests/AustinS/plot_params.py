import argparse
import pandas as pd
import matplotlib.pyplot as plt


# Columns that typically define a "test configuration" (controls)
DEFAULT_CONTROL_COLS = [
    "DelayProfile",
    "HARQ Type",
    "NHARQProcesses",
    "rvSeq",
    "nTxAnts",
    "nRxAnts",
    "NumLayers",
    "Code Rate",
]

def parse_kv_list(kv_list):
    """
    Parse repeated --fix key=value arguments into a dict.
    Values are kept as strings initially, then cast if possible.
    """
    fixes = {}
    if not kv_list:
        return fixes

    for item in kv_list:
        if "=" not in item:
            raise ValueError(f"Bad --fix '{item}'. Use key=value.")
        k, v = item.split("=", 1)
        k = k.strip()
        v = v.strip()

        # Try to cast numeric values
        if v.lower() in ("true", "false"):
            v_cast = v.lower() == "true"
        else:
            try:
                # int first if clean
                if v.isdigit() or (v.startswith("-") and v[1:].isdigit()):
                    v_cast = int(v)
                else:
                    v_cast = float(v)
            except ValueError:
                v_cast = v  # leave as string

        fixes[k] = v_cast
    return fixes


def most_common_value(series: pd.Series):
    s = series.dropna()
    if s.empty:
        return None
    # mode() can return multiple values; take the first
    return s.mode().iloc[0]


def main():
    parser = argparse.ArgumentParser(
        description="Plot all modulation types vs one varied parameter (e.g., SNRdB) for a chosen metric."
    )
    parser.add_argument("--file", default="runlog_results.xlsx", help="Path to the Excel file")
    parser.add_argument("--x", default="SNRdB", help="Varied parameter for x-axis (e.g., SNRdB)")
    parser.add_argument("--y", default="Throughput Efficiency (%)", help="Metric for y-axis")
    parser.add_argument(
        "--controls",
        nargs="*",
        default=DEFAULT_CONTROL_COLS,
        help="Columns treated as controls (held constant). Others may be auto-selected."
    )
    parser.add_argument(
        "--fix",
        action="append",
        default=[],
        help="Force a control value: --fix DelayProfile=TDL-C (can repeat)"
    )
    parser.add_argument(
        "--auto",
        action="store_true",
        help="If a control isn't fixed, auto-select the most common value for a clean slice."
    )
    parser.add_argument(
        "--title",
        default=None,
        help="Custom plot title (optional)"
    )

    args = parser.parse_args()

    df = pd.read_excel(args.file)

    # Basic validation
    for col in [args.x, args.y, "Modulation"]:
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found. Available columns:\n{list(df.columns)}")

    fixes = parse_kv_list(args.fix)

    # Apply forced fixes
    subset = df.copy()
    for k, v in fixes.items():
        if k not in subset.columns:
            raise KeyError(f"Fix column '{k}' not found in the file.")
        subset = subset[subset[k] == v]

    # If auto mode: choose the most common value for any control column not already fixed
    auto_chosen = {}
    if args.auto:
        for c in args.controls:
            if c in subset.columns and c not in fixes:
                mc = most_common_value(subset[c])
                if mc is not None:
                    subset = subset[subset[c] == mc]
                    auto_chosen[c] = mc

    if subset.empty:
        print("No rows match your fixed/auto control selection.")
        print("Try removing some --fix values, or use --auto, or check column names/values.")
        return

    # Plot: one curve per modulation
    plt.figure()
    mods = sorted(subset["Modulation"].dropna().unique().tolist())

    any_plotted = False
    for mod in mods:
        mod_df = subset[subset["Modulation"] == mod].copy()

        # Group by x in case multiple runs exist per x value
        # (e.g., multiple logs with same SNR). Use mean by default.
        mod_df = (
            mod_df.groupby(args.x, as_index=False)[args.y]
            .mean()
            .sort_values(by=args.x)
        )

        if mod_df.empty:
            continue

        plt.plot(mod_df[args.x], mod_df[args.y], marker="o", label=mod)
        any_plotted = True

    if not any_plotted:
        print("No modulation curves were plotted (after filtering).")
        return

    # Build title
    if args.title:
        title = args.title
    else:
        # Summarize control settings used
        control_bits = []
        for k, v in fixes.items():
            control_bits.append(f"{k}={v}")
        for k, v in auto_chosen.items():
            control_bits.append(f"{k}={v} (auto)")
        control_str = ", ".join(control_bits) if control_bits else "no controls fixed"
        title = f"{args.y} vs {args.x}\n({control_str})"

    plt.xlabel(args.x)
    plt.ylabel(args.y)
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
