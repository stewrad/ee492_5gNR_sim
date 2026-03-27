import argparse
import os
from dataclasses import dataclass
from typing import Optional, List, Dict, Any

import pandas as pd
import matplotlib.pyplot as plt

try:
    import yaml  # pip install pyyaml
except ImportError:
    yaml = None


@dataclass
class ChartSpec:
    name: str
    kind: str                 # "line", "scatter", "bar"
    x: str
    y: str
    groupby: Optional[str] = None
    title: Optional[str] = None
    xlabel: Optional[str] = None
    ylabel: Optional[str] = None
    filter_expr: Optional[str] = None  # pandas query string


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
    specs = []
    for c in cfg.get("charts", []):
        specs.append(
            ChartSpec(
                name=c["name"],
                kind=c.get("kind", "line"),
                x=c["x"],
                y=c["y"],
                groupby=c.get("groupby"),
                title=c.get("title"),
                xlabel=c.get("xlabel"),
                ylabel=c.get("ylabel"),
                filter_expr=c.get("filter"),
            )
        )
    return specs


def _ensure_dataframe(df_or_dict: Any) -> pd.DataFrame:
    """
    pd.read_excel can return:
      - DataFrame (normal case)
      - dict[str, DataFrame] if sheet_name=None was used somewhere
    This normalizes to a single DataFrame (first sheet if dict).
    """
    if isinstance(df_or_dict, dict):
        return list(df_or_dict.values())[0]
    return df_or_dict


def _coerce_numeric(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    """
    Try to convert listed columns to numeric.
    Non-numeric values become NaN (then we can drop them for plotting).
    """
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


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

    # Make sure x/y are in the dataframe
    for required in [spec.x, spec.y]:
        if required not in d.columns:
            raise KeyError(
                f"Column '{required}' not found for chart '{spec.name}'. "
                f"Available columns: {list(d.columns)}"
            )

    d = _coerce_numeric(d, [spec.x, spec.y])
    d = d.dropna(subset=[spec.x, spec.y])

    fig, ax = plt.subplots()

    if spec.groupby:
        if spec.groupby not in d.columns:
            raise KeyError(
                f"Groupby column '{spec.groupby}' not found for chart '{spec.name}'. "
                f"Available columns: {list(d.columns)}"
            )

        # Special-case ordering for Modulation
        if spec.groupby == "Modulation":
            mod_order = ["QPSK", "16QAM", "64QAM", "256QAM", "1024QAM"]
            present = set(d[spec.groupby].astype(str).unique())

            # Plot in requested order first
            for mod in mod_order:
                if mod in present:
                    gdf = d[d[spec.groupby].astype(str) == mod].sort_values(spec.x)
                    if spec.kind == "line":
                        ax.plot(gdf[spec.x], gdf[spec.y], label=mod)
                    elif spec.kind == "scatter":
                        ax.scatter(gdf[spec.x], gdf[spec.y], label=mod)
                    elif spec.kind == "bar":
                        ax.bar(gdf[spec.x], gdf[spec.y], label=mod)
                    else:
                        raise ValueError(f"Unknown chart kind: {spec.kind}")

            # Plot any remaining modulation labels not in mod_order (if any)
            remaining = sorted(present - set(mod_order))
            for mod in remaining:
                gdf = d[d[spec.groupby].astype(str) == mod].sort_values(spec.x)
                if spec.kind == "line":
                    ax.plot(gdf[spec.x], gdf[spec.y], label=mod)
                elif spec.kind == "scatter":
                    ax.scatter(gdf[spec.x], gdf[spec.y], label=mod)
                elif spec.kind == "bar":
                    ax.bar(gdf[spec.x], gdf[spec.y], label=mod)
                else:
                    raise ValueError(f"Unknown chart kind: {spec.kind}")

        else:
            # Default groupby behavior (sorted by x within each group)
            for gval, gdf in d.groupby(spec.groupby):
                gdf = gdf.sort_values(spec.x)
                if spec.kind == "line":
                    ax.plot(gdf[spec.x], gdf[spec.y], label=str(gval))
                elif spec.kind == "scatter":
                    ax.scatter(gdf[spec.x], gdf[spec.y], label=str(gval))
                elif spec.kind == "bar":
                    ax.bar(gdf[spec.x], gdf[spec.y], label=str(gval))
                else:
                    raise ValueError(f"Unknown chart kind: {spec.kind}")

        ax.legend(title=spec.groupby)

    else:
        # No groupby: just sort by x and plot one series
        d = d.sort_values(spec.x)
        if spec.kind == "line":
            ax.plot(d[spec.x], d[spec.y])
        elif spec.kind == "scatter":
            ax.scatter(d[spec.x], d[spec.y])
        elif spec.kind == "bar":
            ax.bar(d[spec.x], d[spec.y])
        else:
            raise ValueError(f"Unknown chart kind: {spec.kind}")

    ax.set_title(spec.title or spec.name)
    ax.set_xlabel(spec.xlabel or spec.x)
    ax.set_ylabel(spec.ylabel or spec.y)
    ax.grid(True, alpha=0.3)

    os.makedirs(out_dir, exist_ok=True)

    # Safe-ish filename
    safe_name = "".join(ch if ch.isalnum() or ch in (" ", "_", "-", ".") else "_" for ch in spec.name).strip()
    out_path = os.path.join(out_dir, f"{safe_name}.png")

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

    # Read excel
    if args.sheet:
        df = pd.read_excel(args.excel, sheet_name=args.sheet)
    else:
        df = pd.read_excel(args.excel)

    df = _ensure_dataframe(df)

    cfg = load_config(args.config)
    specs = to_specs(cfg)

    made = []
    for spec in specs:
        made.append(plot_chart(df, spec, args.out))

    print("Created charts:")
    for m in made:
        print(" -", m)


if __name__ == "__main__":
    main()