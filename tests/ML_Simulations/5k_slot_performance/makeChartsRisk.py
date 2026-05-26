import argparse
import os
from dataclasses import dataclass
from typing import Optional, List, Dict, Any

import pandas as pd
import matplotlib.pyplot as plt
import yaml


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


def load_config(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def safe_filename(name: str) -> str:
    return "".join(
        ch if ch.isalnum() or ch in (" ", "_", "-", ".") else "_"
        for ch in name
    ).strip()


def make_risk_plot_data(df: pd.DataFrame) -> pd.DataFrame:
    # Clean numeric columns
    df["SNRdB"] = pd.to_numeric(df["SNRdB"], errors="coerce")
    df["NHARQProcesses"] = pd.to_numeric(df["NHARQProcesses"], errors="coerce")

    # Make sure text columns are clean
    df["Modulation"] = df["Modulation"].astype(str).str.strip()

    # -------------------------
    # Modulation-specific risk
    # -------------------------
    mod_df = df[
        df["Modulation"].isin(["QPSK", "16QAM", "64QAM"])
    ].copy()

    mod_df["Method"] = mod_df["Modulation"]
    mod_df["Risk Probability"] = pd.to_numeric(
        mod_df["Predicted Risk Probability"],
        errors="coerce"
    )

    mod_df = mod_df[["SNRdB", "NHARQProcesses", "Method", "Risk Probability"]]

    # -------------------------
    # ML recommended risk
    # -------------------------
    ml_df = df[
        df["Is Recommended"] == 1
    ].copy()

    ml_df["Method"] = "ML"
    ml_df["Risk Probability"] = pd.to_numeric(
        ml_df["Recommended Risk Probability"],
        errors="coerce"
    )

    ml_df = ml_df[["SNRdB", "NHARQProcesses", "Method", "Risk Probability"]]

    # Combine
    plot_df = pd.concat([mod_df, ml_df], ignore_index=True)

    # Average across NHARQProcesses for each SNR/method
    plot_df = (
        plot_df
        .dropna(subset=["SNRdB", "Risk Probability"])
        .groupby(["Method", "SNRdB"], as_index=False)["Risk Probability"]
        .mean()
        .sort_values(["Method", "SNRdB"])
    )

    return plot_df


def plot_chart(plot_df: pd.DataFrame, spec: ChartSpec, out_dir: str) -> str:
    fig, ax = plt.subplots(figsize=(7, 4.5))

    for method, gdf in plot_df.groupby(spec.groupby):
        gdf = gdf.sort_values(spec.x)
        ax.plot(
            gdf[spec.x],
            gdf[spec.y],
            marker="o",
            label=str(method)
        )

    ax.set_title(spec.title or spec.name)
    ax.set_xlabel(spec.xlabel or spec.x)
    ax.set_ylabel(spec.ylabel or spec.y)
    ax.grid(True, alpha=0.3)
    ax.legend(title=spec.groupby)

    os.makedirs(out_dir, exist_ok=True)

    out_path = os.path.join(out_dir, f"{safe_filename(spec.name)}.png")
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)

    return out_path


def main():
    parser = argparse.ArgumentParser(description="Generate risk probability chart.")
    parser.add_argument("--excel", required=True, help="Input Excel file")
    parser.add_argument("--sheet", default="All_Original_Rows", help="Sheet name")
    parser.add_argument("--config", required=True, help="Risk chart YAML file")
    parser.add_argument("--out", default="Figures", help="Output folder")
    args = parser.parse_args()

    df = pd.read_excel(args.excel, sheet_name=args.sheet)

    cfg = load_config(args.config)
    chart = cfg["charts"][0]

    spec = ChartSpec(
        name=chart["name"],
        kind=chart.get("kind", "line"),
        x=chart["x"],
        y=chart["y"],
        groupby=chart.get("groupby"),
        title=chart.get("title"),
        xlabel=chart.get("xlabel"),
        ylabel=chart.get("ylabel"),
    )

    plot_df = make_risk_plot_data(df)
    made = plot_chart(plot_df, spec, args.out)

    print("Created chart:")
    print(" -", made)


if __name__ == "__main__":
    main()