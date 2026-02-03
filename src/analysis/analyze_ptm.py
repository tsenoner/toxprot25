#!/usr/bin/env python3
"""
PTM Analysis for ToxProt Dataset

Analyzes and visualizes post-translational modifications (PTMs) from the
PTM Summary column across ToxProt datasets at multiple time points.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .colors import YEAR_COLORS

# Default comparison years (10-year intervals)
DEFAULT_YEARS = [2005, 2015, 2025]

# Top 10 most common PTM types, ordered by frequency
PTM_TYPES = [
    "Disulfide bond",
    "Amidation",
    "Glycosylation",
    "Hydroxylation",
    "Pyrrolidone carboxylic acid",
    "Gamma-carboxyglutamic acid",
    "D-amino acid",
    "Bromination",
    "Sulfation",
    "Lipidation",
]

# Display labels (pluralized)
PTM_LABELS = {ptm: ptm + "s" if not ptm.endswith("acid") else ptm for ptm in PTM_TYPES}
PTM_LABELS["Disulfide bond"] = "Disulfide bonds"


def parse_ptm_summary(ptm_str: str) -> dict[str, int]:
    """Parse PTM Summary string into a dictionary of counts."""
    if pd.isna(ptm_str) or ptm_str == "":
        return {}
    ptm_dict = {}
    for part in ptm_str.split("; "):
        if ":" in part:
            ptm_type, count = part.rsplit(":", 1)
            ptm_dict[ptm_type] = int(count)
        else:
            ptm_dict[part] = 1
    return ptm_dict


def load_data(file_path: Path) -> pd.DataFrame:
    """Load CSV data and parse PTM information."""
    df = pd.read_csv(file_path)
    df["PTM_dict"] = df["PTM Summary"].apply(parse_ptm_summary)

    for ptm_type in PTM_TYPES:
        df[f"has_{ptm_type}"] = df["PTM_dict"].apply(lambda x, pt=ptm_type: 1 if pt in x else 0)
        df[f"count_{ptm_type}"] = df["PTM_dict"].apply(lambda x, pt=ptm_type: x.get(pt, 0))

    return df


def create_combined_figure(datasets: dict[int, pd.DataFrame], output_dir: Path, top_n: int = 10):
    """Create combined figure with PTM type comparison (A) and count distributions (B).

    Args:
        datasets: Dictionary mapping year to DataFrame
        output_dir: Directory to save output files
        top_n: Number of top PTM types to show (default 10)
    """
    years = sorted(datasets.keys())
    fig = plt.figure(figsize=(18, 16))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.4], hspace=0.35)

    # ========== PLOT A: PTM Type Frequency ==========
    ax_a = fig.add_subplot(gs[0])

    ptm_data = []
    for ptm_type in PTM_TYPES[:top_n]:
        row_data = {"PTM Type": PTM_LABELS[ptm_type]}
        for year in years:
            df = datasets[year]
            count = df[f"has_{ptm_type}"].sum()
            row_data[f"Count {year}"] = count
        ptm_data.append(row_data)

    df_ptm = pd.DataFrame(ptm_data).sort_values(f"Count {years[-1]}", ascending=False)
    y_pos = np.arange(len(df_ptm))

    # Plot stacked bars (latest year at back, earliest at front)
    for year in reversed(years):
        ax_a.barh(
            y_pos,
            df_ptm[f"Count {year}"],
            color=YEAR_COLORS[year],
            edgecolor="black",
            label=str(year),
        )

    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(df_ptm["PTM Type"])
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Number of Proteins", fontsize=16)
    years_str = ", ".join(str(y) for y in years)
    ax_a.set_title(f"PTM Type Frequency ({years_str})", fontsize=16, fontweight="bold", pad=10)
    ax_a.tick_params(axis="both", labelsize=14)
    ax_a.text(-0.1, 1.02, "A", transform=ax_a.transAxes, fontsize=28, fontweight="bold")
    ax_a.set_axisbelow(True)
    ax_a.grid(axis="x", ls="--", lw=2, alpha=0.5)

    # Annotations
    for row_idx, (_, row) in enumerate(df_ptm.iterrows()):
        counts_str = " → ".join(f"{int(row[f'Count {y}'])}" for y in years)
        text_x = max(row[f"Count {y}"] for y in years) + ax_a.get_xlim()[1] * 0.01
        ax_a.text(text_x, y_pos[row_idx], counts_str, va="center", fontsize=12, fontweight="bold")

    ax_a.set_xlim(ax_a.get_xlim()[0], ax_a.get_xlim()[1] * 1.12)

    # ========== PLOT B: PTM Count Distributions (2x3 grid) ==========
    gs_b = gs[1].subgridspec(2, 3, wspace=0.15, hspace=0.35)
    ptm_types_to_plot = PTM_TYPES[:6]
    num_bins = 10

    # Collect histogram data
    all_hist_data = []
    for ptm_type in ptm_types_to_plot:
        year_hists = {}
        for year in years:
            df = datasets[year]
            counts = df[df[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]
            hist = np.array(
                [(counts == i).sum() for i in range(1, num_bins)] + [(counts >= num_bins).sum()]
            )
            year_hists[year] = hist
        all_hist_data.append(year_hists)

    axes_b = []
    for idx, ptm_type in enumerate(ptm_types_to_plot):
        row, col = idx // 3, idx % 3
        ax = fig.add_subplot(gs_b[row, col])
        axes_b.append(ax)

        year_hists = all_hist_data[idx]
        x_pos = np.arange(1, 11)

        for year in reversed(years):
            ax.bar(
                x_pos,
                year_hists[year],
                width=0.8,
                color=YEAR_COLORS[year],
                edgecolor="black",
                label=str(year) if idx == 0 else None,
            )

        ax.set_title(PTM_LABELS[ptm_type], fontsize=14, fontweight="bold", pad=8)
        ax.set_xlim(0.5, 10.5)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([str(i) for i in range(1, 10)] + ["10+"], fontsize=11)
        ax.tick_params(axis="both", labelsize=12)

        if row == 1:
            ax.set_xlabel("PTMs per Protein", fontsize=14)
        if col == 0:
            ax.set_ylabel("Number of Proteins", fontsize=14)
        if idx == 0:
            ax.text(-0.35, 1.08, "B", transform=ax.transAxes, fontsize=28, fontweight="bold")

        ax.set_axisbelow(True)
        ax.grid(axis="y", ls="--", lw=1, alpha=0.5)

    # Set y-axis limits
    for ax, year_hists in zip(axes_b, all_hist_data, strict=True):
        max_y = max(h.max() for h in year_hists.values())
        ax.set_ylim(0, max_y * 1.05)

    # Common legend
    handles, labels = ax_a.get_legend_handles_labels()
    fig.legend(
        handles[::-1],
        labels[::-1],
        loc="center",
        fontsize=22,
        ncol=len(years),
        frameon=False,
        bbox_to_anchor=(0.5, 0.54),
        columnspacing=1.0,
        handletextpad=0.4,
    )

    plt.savefig(output_dir / "ptm_overview.png", dpi=300, bbox_inches="tight")
    plt.close()


def main():
    """Main analysis pipeline."""
    data_dir = Path("data/processed/toxprot")
    output_dir = Path("figures/ptm")
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = {}
    for year in DEFAULT_YEARS:
        datasets[year] = load_data(data_dir / f"toxprot_{year}.csv")

    print(f"Loaded {sum(len(df) for df in datasets.values()):,} entries")
    create_combined_figure(datasets, output_dir)
    print(f"Saved to {output_dir}")


if __name__ == "__main__":
    main()
