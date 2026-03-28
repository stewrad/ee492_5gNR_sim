import argparse
import os
from dataclasses import dataclass
from typing import Optional, List, Dict, Any, Union

import pandas as pd
import matplotlib.pyplot as plt

try:
    import yaml  # pip install pyyaml
except ImportError:
    yaml = None


@dataclass
class ChartSpec:
    name: str
    kind: str  # "line", "scatter", "bar"
    x: str
    y: Union[str, List[str]]  # can be a single column name or a list of column names
    groupby: Optional[str] = None
    title: Optional[str] = None
    xlabel: Optional[str] = None
    ylabel: Optional[str] = None
    filter_expr: Optional[str] = None  # pandas query string
    ylabels: Optional[List[str]] = None  # labels corresponding to y list (e.g., ["Pre-HARQ","Post-HARQ"])


def load_config(path: str) -> Dict[str, Any]:
    ext = os.path.splitext(path)[1].lower()
    if ext in [".yml", ".yaml"]:
        if yaml is None:
            raise RuntimeError("pyyaml not installed. Run: pip install pyyaml")
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    elif ext == ".json":
        import json
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    else:
        raise ValueError("Config must be .yml/.yaml or .json")


def to_specs(cfg: Dict[str, Any]) -> List[ChartSpec]:
    specs: List[ChartSpec] = []
    for c in cfg.get("charts", []):
        specs.append(
            ChartSpec(
                name=c["name"],
                kind=str(c.get("kind", "line")).lower(),
                x=c["x"],
                y=c["y"],
                groupby=c.get("groupby"),
                title=c.get("title"),
                xlabel=c.get("xlabel"),
                ylabel=c.get("ylabel"),
                filter_expr=c.get("filter"),
                ylabels=c.get("ylabels"),
            )
        )
    return specs


def _ensure_dataframe(df_or_dict: Any) -> pd.DataFrame:
    """
    pd.read_excel returns:
      - DataFrame normally
      - dict[str, DataFrame] if sheet_name=None is used
    Normalize to a single DataFrame (first sheet if dict).
    """
    if isinstance(df_or_dict, dict):
        return list(df_or_dict.values())[0]
    return df_or_dict


def _coerce_numeric(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    """
    Convert specified columns to numeric when possible.
    Non-numeric values become NaN.
    """
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _safe_filename(name: str) -> str:
    return "".join(ch if ch.isalnum() or ch in (" ", "_", "-", ".") else "_" for ch in name).strip()


def plot_chart(df: pd.DataFrame, spec: ChartSpec, out_dir: str) -> str:
    d = df.copy()

    # Apply filter (pandas query)
    if spec.filter_expr:
        try:
            d = d.query(spec.filter_expr)
        except Exception as e:
            raise RuntimeError(
                f"Filter query failed for chart '{spec.name}'.\n"
                f"Filter: {spec.filter_expr}\n"
                f"Error: {e}"
            )

    if d.empty:
        raise RuntimeError(
            f"No rows matched for chart '{spec.name}'.\n"
            f"Check your filter / column names.\n"
            f"Filter: {spec.filter_expr}"
        )

    # y can be a single column name or list of column names
    y_cols: List[str] = spec.y if isinstance(spec.y, list) else [spec.y]
    y_names: List[str] = spec.ylabels if spec.ylabels else y_cols

    if len(y_names) != len(y_cols):
        raise RuntimeError(
            f"Chart '{spec.name}': 'ylabels' must be the same length as 'y' when y is a list.\n"
            f"y: {y_cols}\n"
            f"ylabels: {y_names}"
        )

    # Validate columns exist
    required_cols = [spec.x] + y_cols
    for col in required_cols:
        if col not in d.columns:
            raise KeyError(
                f"Column '{col}' not found for chart '{spec.name}'. "
                f"Available columns: {list(d.columns)}"
            )

    # Coerce x and y columns to numeric and drop NaNs
    d = _coerce_numeric(d, [spec.x] + y_cols)
    d = d.dropna(subset=[spec.x] + y_cols)

    if d.empty:
        raise RuntimeError(
            f"After converting to numeric, no valid rows remain for chart '{spec.name}'.\n"
            f"This usually means the x/y columns contain non-numeric text or blanks.\n"
            f"x: {spec.x}, y: {y_cols}"
        )

    fig, ax = plt.subplots()

    # Plot helpers
    def _plot_series(xvals, yvals, label: str):
        if spec.kind == "line":
            ax.plot(xvals, yvals, label=label)
        elif spec.kind == "scatter":
            ax.scatter(xvals, yvals, label=label)
        elif spec.kind == "bar":
            ax.bar(xvals, yvals, label=label)
        else:
            raise ValueError(f"Unknown chart kind: {spec.kind} (use line/scatter/bar)")

    if spec.groupby:
        if spec.groupby not in d.columns:
            raise KeyError(
                f"Groupby column '{spec.groupby}' not found for chart '{spec.name}'. "
                f"Available columns: {list(d.columns)}"
            )

        # Special ordering for Modulation
        if spec.groupby == "Modulation":
            mod_order = ["QPSK", "16QAM", "64QAM", "256QAM", "1024QAM"]
            present = list(d[spec.groupby].astype(str).unique())

            # Plot modulations in preferred order first, then any leftovers
            ordered_mods = [m for m in mod_order if m in present] + sorted([m for m in present if m not in mod_order])

            for mod in ordered_mods:
                gdf = d[d[spec.groupby].astype(str) == mod].sort_values(spec.x)

                # For each requested y column, plot as its own series with suffix label
                for ycol, ylab in zip(y_cols, y_names):
                    series_label = f"{mod} {ylab}" if spec.ylabels else f"{mod} {ycol}"
                    _plot_series(gdf[spec.x], gdf[ycol], series_label)

        else:
            # Default groupby behavior (sorted by x within each group)
            for gval, gdf in d.groupby(spec.groupby):
                gdf = gdf.sort_values(spec.x)
                for ycol, ylab in zip(y_cols, y_names):
                    series_label = f"{gval} {ylab}" if spec.ylabels else f"{gval} {ycol}"
                    _plot_series(gdf[spec.x], gdf[ycol], str(series_label))

        ax.legend(title=spec.groupby, fontsize=8, title_fontsize=8, loc="lower right", bbox_to_anchor=(1, 0.10))

    else:
        # No groupby: just sort by x and plot one or multiple y series
        d = d.sort_values(spec.x)
        for ycol, ylab in zip(y_cols, y_names):
            series_label = ylab if spec.ylabels else ycol
            _plot_series(d[spec.x], d[ycol], str(series_label))
        if len(y_cols) > 1:
            ax.legend()

    ax.set_title(spec.title or spec.name)
    ax.set_xlabel(spec.xlabel or spec.x)
    ax.set_ylabel(spec.ylabel or spec.y if isinstance(spec.y, str) else (spec.ylabel or "Value"))
    ax.grid(True, alpha=0.3)

    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{_safe_filename(spec.name)}.png")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def main():
    p = argparse.ArgumentParser(description="Generate charts from Excel using a config file.")
    p.add_argument("--excel", required=True, help="Path to input .xlsx")
    p.add_argument("--sheet", default=None, help="Sheet name (default: first sheet)")
    p.add_argument("--config", required=True, help="Path to config .yml/.yaml or .json")
    p.add_argument("--out", default="charts_out", help="Output directory for PNGs")
    args = p.parse_args()

    # Read Excel
    if args.sheet:
        df = pd.read_excel(args.excel, sheet_name=args.sheet)
    else:
        df = pd.read_excel(args.excel)

    df = _ensure_dataframe(df)

    cfg = load_config(args.config)
    specs = to_specs(cfg)

    made: List[str] = []
    for spec in specs:
        made.append(plot_chart(df, spec, args.out))

    print("Created charts:")
    for m in made:
        print(" -", m)


if __name__ == "__main__":
    main()