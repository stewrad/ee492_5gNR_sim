import argparse
import os
from dataclasses import dataclass
from typing import Optional, List, Dict, Any

import pandas as pd
import matplotlib.pyplot as plt

try:
    import yaml
except ImportError:
    yaml = None


@dataclass
class ChartSpec:
    name: str
    kind: str
    x: str
    y: str
    groupby: Optional[str] = None
    title: Optional[str] = None
    xlabel: Optional[str] = None
    ylabel: Optional[str] = None
    filter_expr: Optional[str] = None
    yscale: Optional[str] = None


def load_config(path: str) -> Dict[str, Any]:
    if yaml is None:
        raise RuntimeError("pyyaml not installed. Run: pip install pyyaml")

    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


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
                yscale=c.get("yscale"),
            )
        )
    return specs


def coerce_numeric(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def safe_filename(name: str) -> str:
    return "".join(
        ch if ch.isalnum() or ch in (" ", "_", "-", ".") else "_"
        for ch in name
    ).strip()


def plot_chart(df: pd.DataFrame, spec: ChartSpec, out_dir: str) -> str:
    d = df.copy()

    if spec.filter_expr:
        d = d.query(spec.filter_expr)

    required_cols = [spec.x, spec.y]
    if spec.groupby:
        required_cols.append(spec.groupby)

    for col in required_cols:
        if col not in d.columns:
            raise KeyError(
                f"Column '{col}' not found for chart '{spec.name}'. "
                f"Available columns: {list(d.columns)}"
            )

    d = coerce_numeric(d, [spec.x, spec.y])
    d = d.dropna(subset=[spec.x, spec.y])

    if d.empty:
        raise RuntimeError(f"No valid rows found for chart '{spec.name}'.")

    fig, ax = plt.subplots(figsize=(7, 4.5))

    if spec.groupby:
        # Average duplicate points, e.g. multiple HARQ values at same SNR
        grouped = (
            d.groupby([spec.groupby, spec.x], as_index=False)[spec.y]
            .mean()
            .sort_values([spec.groupby, spec.x])
        )

        for gval, gdf in grouped.groupby(spec.groupby):
            gdf = gdf.sort_values(spec.x)

            if spec.kind == "line":
                ax.plot(gdf[spec.x], gdf[spec.y], marker="o", label=str(gval))
            elif spec.kind == "scatter":
                ax.scatter(gdf[spec.x], gdf[spec.y], label=str(gval))
            elif spec.kind == "bar":
                ax.bar(gdf[spec.x], gdf[spec.y], label=str(gval))
            else:
                raise ValueError(f"Unknown chart kind: {spec.kind}")

        ax.legend(title=spec.groupby)

    else:
        # Average duplicate x values
        grouped = (
            d.groupby(spec.x, as_index=False)[spec.y]
            .mean()
            .sort_values(spec.x)
        )

        if spec.kind == "line":
            ax.plot(grouped[spec.x], grouped[spec.y], marker="o")
        elif spec.kind == "scatter":
            ax.scatter(grouped[spec.x], grouped[spec.y])
        elif spec.kind == "bar":
            ax.bar(grouped[spec.x], grouped[spec.y])
        else:
            raise ValueError(f"Unknown chart kind: {spec.kind}")

    if spec.yscale:
        ax.set_yscale(spec.yscale)

    ax.set_title(spec.title or spec.name)
    ax.set_xlabel(spec.xlabel or spec.x)
    ax.set_ylabel(spec.ylabel or spec.y)
    ax.grid(True, alpha=0.3)

    os.makedirs(out_dir, exist_ok=True)

    out_path = os.path.join(out_dir, f"{safe_filename(spec.name)}.png")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    return out_path


def main():
    parser = argparse.ArgumentParser(description="Generate charts from Excel using charts.yml.")
    parser.add_argument("--excel", required=True, help="Input Excel file")
    parser.add_argument("--sheet", default=None, help="Sheet name, optional")
    parser.add_argument("--config", required=True, help="charts.yml file")
    parser.add_argument("--out", default="Figures", help="Output folder")
    args = parser.parse_args()

    if args.sheet:
        df = pd.read_excel(args.excel, sheet_name=args.sheet)
    else:
        df = pd.read_excel(args.excel)

    cfg = load_config(args.config)
    specs = to_specs(cfg)

    made = []
    for spec in specs:
        made.append(plot_chart(df, spec, args.out))

    print("Created charts:")
    for path in made:
        print(" -", path)


if __name__ == "__main__":
    main()