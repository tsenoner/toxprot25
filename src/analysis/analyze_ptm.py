#!/usr/bin/env python3
"""
PTM Analysis for ToxProt Dataset

Analyzes and visualizes post-translational modifications (PTMs) from the
PTM Summary column across ToxProt datasets.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

# Years available in the dataset
YEARS = [str(y) for y in range(2005, 2026)]

# PTM types - top 10 most common, ordered by frequency
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

# Labels map to display names (pluralized where appropriate)
PTM_LABELS = {ptm: ptm + "s" if not ptm.endswith("acid") else ptm for ptm in PTM_TYPES}
PTM_LABELS["Disulfide bond"] = "Disulfide bonds"

# Colors (using tab10 colormap for consistency)
COLOR_OLD = plt.cm.tab10(0)
COLOR_NEW = plt.cm.tab10(1)


def parse_ptm_summary(ptm_str):
    """Parse PTM Summary string into a dictionary of counts."""
    if pd.isna(ptm_str) or ptm_str == "":
        return {}
    ptm_dict = {}
    for part in ptm_str.split("; "):
        if ":" in part:
            ptm_type, count = part.rsplit(":", 1)  # rsplit for safety
            ptm_dict[ptm_type] = int(count)
        else:
            ptm_dict[part] = 1  # No count means 1
    return ptm_dict


def load_data(file_path):
    """Load CSV data and parse PTM information."""
    df = pd.read_csv(file_path)
    df["PTM_dict"] = df["PTM Summary"].apply(parse_ptm_summary)

    for ptm_type in PTM_TYPES:
        df[f"has_{ptm_type}"] = df["PTM_dict"].apply(lambda x, pt=ptm_type: 1 if pt in x else 0)
        df[f"count_{ptm_type}"] = df["PTM_dict"].apply(lambda x, pt=ptm_type: x.get(pt, 0))

    df["total_ptm_count"] = df["PTM_dict"].apply(lambda x: sum(x.values()))
    df["has_any_ptm"] = df["PTM_dict"].apply(lambda x: len(x) > 0)

    return df


def create_combined_figure(df_old, df_new, output_dir, year_old, year_new, top_n=5):
    """Create combined figure with Plot A (type comparison) and Plot B (distributions)."""
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.4)

    # ========== PLOT A: PTM Type Frequency ==========
    ax_a = fig.add_subplot(gs[0])

    # Prepare data for top_n PTM types
    ptm_data = []
    for ptm_type in PTM_TYPES[:top_n]:
        count_old = df_old[f"has_{ptm_type}"].sum()
        count_new = df_new[f"has_{ptm_type}"].sum()
        pct_old = (count_old / len(df_old)) * 100
        pct_new = (count_new / len(df_new)) * 100
        ptm_data.append(
            {
                "PTM Type": PTM_LABELS[ptm_type],
                f"Count {year_old}": count_old,
                f"Count {year_new}": count_new,
                f"Pct {year_old}": pct_old,
                f"Pct {year_new}": pct_new,
            }
        )

    df_ptm = pd.DataFrame(ptm_data).sort_values(f"Count {year_new}", ascending=False)

    # Stacked bars
    values_base = np.minimum(df_ptm[f"Count {year_old}"], df_ptm[f"Count {year_new}"])
    values_increase = np.maximum(0, df_ptm[f"Count {year_new}"] - df_ptm[f"Count {year_old}"])
    y_pos = np.arange(len(df_ptm))

    ax_a.barh(y_pos, values_base, color=COLOR_OLD, edgecolor="black")
    ax_a.barh(y_pos, values_increase, left=values_base, color=COLOR_NEW, edgecolor="black")

    # Formatting
    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(df_ptm["PTM Type"])
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Number of Proteins", fontsize=14)
    ax_a.set_title(
        f"PTM Type Frequency ({year_old} vs {year_new})", fontsize=14, fontweight="bold", pad=10
    )
    ax_a.tick_params(axis="both", labelsize=12)

    # Add panel label A
    ax_a.text(
        -0.1,
        0.98,
        "A",
        transform=ax_a.transAxes,
        fontsize=28,
        fontweight="bold",
        va="bottom",
        ha="left",
    )
    ax_a.set_axisbelow(True)
    ax_a.grid(axis="x", ls="--", lw=2, alpha=0.5)

    # Annotations
    for i, row in df_ptm.iterrows():
        idx = list(df_ptm.index).index(i)
        text_x = row[f"Count {year_new}"] + ax_a.get_xlim()[1] * 0.01
        annotation = (
            f"{int(row[f'Count {year_old}'])} ({row[f'Pct {year_old}']:.1f}%) → "
            f"{int(row[f'Count {year_new}'])} ({row[f'Pct {year_new}']:.1f}%)"
        )
        ax_a.text(
            text_x,
            y_pos[idx],
            annotation,
            va="center",
            ha="left",
            fontsize=10,
            fontweight="bold",
        )

    ax_a.set_xlim(ax_a.get_xlim()[0], ax_a.get_xlim()[1] * 1.175)

    # Add legend to Plot A
    legend_elements = [
        Patch(facecolor=COLOR_OLD, label=str(year_old), edgecolor="black"),
        Patch(facecolor=COLOR_NEW, label=str(year_new), edgecolor="black"),
    ]
    ax_a.legend(handles=legend_elements, loc="lower right", fontsize=18)

    # ========== PLOT B: PTM Count Distributions ==========
    gs_b = gs[1].subgridspec(1, 3, wspace=0.15)
    # Top 3 PTM types for distribution plots
    ptm_types_to_plot = PTM_TYPES[:3]

    # First, collect all data to determine shared y-axis range
    all_hist_data = []
    for ptm_type in ptm_types_to_plot:
        counts_old = df_old[df_old[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]
        counts_new = df_new[df_new[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]

        num_bins = 10
        hist_old = np.array(
            [(counts_old == i).sum() for i in range(1, num_bins)] + [(counts_old >= num_bins).sum()]
        )
        hist_new = np.array(
            [(counts_new == i).sum() for i in range(1, num_bins)] + [(counts_new >= num_bins).sum()]
        )
        all_hist_data.append((hist_old, hist_new))

    # Create subplots with shared y-axis
    axes_b = []
    for idx, ptm_type in enumerate(ptm_types_to_plot):
        ax = fig.add_subplot(gs_b[idx])
        axes_b.append(ax)

        hist_old, hist_new = all_hist_data[idx]

        # Plot overlapping bars
        x_pos = np.arange(1, 11)
        ax.bar(
            x_pos,
            hist_new,
            width=1.0,
            color=COLOR_NEW,
            edgecolor="black",
            align="edge",
        )
        ax.bar(
            x_pos,
            hist_old,
            width=1.0,
            color=COLOR_OLD,
            edgecolor="black",
            align="edge",
        )

        # Formatting
        ax.set_title(PTM_LABELS[ptm_type], fontsize=14, fontweight="bold", pad=10)
        ax.set_xlabel("Number of PTMs", fontsize=14)
        ax.set_xlim(1, 11)
        ax.set_xticks([i + 0.5 for i in range(1, 11)])
        ax.set_xticklabels([str(i) for i in range(1, 10)] + ["10+"])
        ax.tick_params(axis="both", labelsize=12)

        if idx == 0:
            ax.set_ylabel("Number of Proteins", fontsize=14)
            # Add panel label B to the leftmost subplot
            ax.text(
                -0.33,
                1.02,
                "B",
                transform=ax.transAxes,
                fontsize=28,
                fontweight="bold",
                va="bottom",
                ha="left",
            )

        ax.set_axisbelow(True)
        ax.grid(axis="y", ls="--", lw=2, alpha=0.5)

    # Set shared y-axis limits
    max_y = max(max(hist_old.max(), hist_new.max()) for hist_old, hist_new in all_hist_data)
    for ax in axes_b:
        ax.set_ylim(0, max_y * 1.05)

    # Save
    plt.savefig(output_dir / "ptm_overview.png", dpi=300, bbox_inches="tight")
    plt.savefig(output_dir / "ptm_overview.pdf", bbox_inches="tight")
    plt.close()

    # Save data
    df_ptm.to_csv(output_dir / "ptm_type_comparison_data.csv", index=False)

    # Distribution data
    dist_data = []
    for ptm_type in ptm_types_to_plot:
        counts_old = df_old[df_old[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]
        counts_new = df_new[df_new[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]

        for i in range(1, 10):
            dist_data.append(
                {
                    "PTM_Type": PTM_LABELS[ptm_type],
                    "PTM_Count": str(i),
                    f"Proteins_{year_old}": (counts_old == i).sum(),
                    f"Proteins_{year_new}": (counts_new == i).sum(),
                }
            )
        dist_data.append(
            {
                "PTM_Type": PTM_LABELS[ptm_type],
                "PTM_Count": "10+",
                f"Proteins_{year_old}": (counts_old >= 10).sum(),
                f"Proteins_{year_new}": (counts_new >= 10).sum(),
            }
        )

    pd.DataFrame(dist_data).to_csv(output_dir / "ptm_count_distribution_data.csv", index=False)


def generate_summary_table(df_old, df_new, output_dir, year_old, year_new):
    """Generate summary statistics table."""
    table_data = [
        ["Metric", str(year_old), str(year_new), "Change"],
        [
            "Total Proteins",
            f"{len(df_old):,}",
            f"{len(df_new):,}",
            f"+{len(df_new) - len(df_old):,}",
        ],
        [
            "Proteins with PTMs",
            f"{df_old['has_any_ptm'].sum():,}",
            f"{df_new['has_any_ptm'].sum():,}",
            f"+{df_new['has_any_ptm'].sum() - df_old['has_any_ptm'].sum():,}",
        ],
        [
            "PTM Coverage (%)",
            f"{(df_old['has_any_ptm'].sum() / len(df_old) * 100):.2f}%",
            f"{(df_new['has_any_ptm'].sum() / len(df_new) * 100):.2f}%",
            f"{(df_new['has_any_ptm'].sum() / len(df_new) - df_old['has_any_ptm'].sum() / len(df_old)) * 100:+.2f}%",
        ],
        [
            "Mean PTMs/Protein",
            f"{df_old[df_old['has_any_ptm']]['total_ptm_count'].mean():.2f}",
            f"{df_new[df_new['has_any_ptm']]['total_ptm_count'].mean():.2f}",
            f"{df_new[df_new['has_any_ptm']]['total_ptm_count'].mean() - df_old[df_old['has_any_ptm']]['total_ptm_count'].mean():+.2f}",
        ],
        ["", "", "", ""],
        ["PTM Type", str(year_old), str(year_new), "Change"],
    ]

    for ptm_type in PTM_TYPES:
        count_old = df_old[f"has_{ptm_type}"].sum()
        count_new = df_new[f"has_{ptm_type}"].sum()
        pct_old = (count_old / len(df_old)) * 100
        pct_new = (count_new / len(df_new)) * 100
        table_data.append(
            [
                PTM_LABELS[ptm_type],
                f"{count_old:,} ({pct_old:.1f}%)",
                f"{count_new:,} ({pct_new:.1f}%)",
                f"+{count_new - count_old:,}",
            ]
        )

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.axis("tight")
    ax.axis("off")

    table = ax.table(
        cellText=table_data,
        cellLoc="left",
        loc="center",
        colWidths=[0.4, 0.2, 0.2, 0.2],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 2.5)

    for i in [0, 6]:
        for j in range(4):
            table[(i, j)].set_facecolor("#34495e")
            table[(i, j)].set_text_props(weight="bold", color="white")

    for j in range(4):
        table[(5, j)].set_facecolor("#ecf0f1")
        table[(5, j)].set_alpha(0.3)

    for i in range(1, len(table_data)):
        if i not in [5, 6]:
            for j in range(4):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor("#f8f9fa")

    plt.title(
        f"PTM Summary Statistics ({year_old} vs {year_new})",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    plt.savefig(output_dir / "ptm_summary_table.png", dpi=300, bbox_inches="tight")
    plt.close()

    pd.DataFrame(table_data[1:], columns=table_data[0]).to_csv(
        output_dir / "ptm_summary_statistics.csv", index=False
    )


def main():
    """Main analysis pipeline."""
    print("PTM Analysis")
    print("=" * 40)

    # Setup paths
    data_dir = Path("data/processed/toxprot")
    output_dir = Path("figures/ptm")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Default years for comparison
    year_old = 2017
    year_new = 2025

    # Load data
    print("Loading data...")
    df_old = load_data(data_dir / f"toxprot_{year_old}.csv")
    df_new = load_data(data_dir / f"toxprot_{year_new}.csv")

    # Generate outputs
    print("Generating figures...")
    create_combined_figure(df_old, df_new, output_dir, year_old, year_new)
    generate_summary_table(df_old, df_new, output_dir, year_old, year_new)

    print(f"Complete. Output: {output_dir.absolute()}")
    print("=" * 40)


if __name__ == "__main__":
    main()
