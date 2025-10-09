#!/usr/bin/env python3
"""
PTM Analysis for ToxProt Dataset

Analyzes and visualizes post-translational modifications (PTMs) from the
PTM Summary column across 2017 and 2025 ToxProt datasets.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

# PTM definitions
PTM_TYPES = ["DISULFID", "CARBOHYD", "MOD_RES", "CROSSLNK", "LIPID"]
PTM_LABELS = {
    "DISULFID": "Disulfide bonds",
    "CARBOHYD": "Glycosylation",
    "MOD_RES": "Modified residues",
    "CROSSLNK": "Cross-links",
    "LIPID": "Lipidation",
}

# Colors
COLOR_2017 = plt.cm.tab10(0)
COLOR_2025 = plt.cm.tab10(1)


def parse_ptm_summary(ptm_str):
    """Parse PTM Summary string into a dictionary of counts"""
    if pd.isna(ptm_str) or ptm_str == "":
        return {}
    ptm_dict = {}
    for part in ptm_str.split("; "):
        if ":" in part:
            ptm_type, count = part.split(":")
            ptm_dict[ptm_type] = int(count)
    return ptm_dict


def load_data(file_path):
    """Load TSV data and parse PTM information"""
    df = pd.read_csv(file_path, sep="\t")
    df["PTM_dict"] = df["PTM Summary"].apply(parse_ptm_summary)

    for ptm_type in PTM_TYPES:
        df[f"has_{ptm_type}"] = df["PTM_dict"].apply(
            lambda x: 1 if ptm_type in x else 0
        )
        df[f"count_{ptm_type}"] = df["PTM_dict"].apply(lambda x: x.get(ptm_type, 0))

    df["total_ptm_count"] = df["PTM_dict"].apply(lambda x: sum(x.values()))
    df["has_any_ptm"] = df["PTM_dict"].apply(lambda x: len(x) > 0)

    return df


def create_combined_figure(df_2017, df_2025, output_dir):
    """Create combined figure with Plot A (type comparison) and Plot B (distributions)"""

    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1], hspace=0.4)

    # ========== PLOT A: PTM Type Frequency ==========
    ax_a = fig.add_subplot(gs[0])

    # Prepare data
    ptm_data = []
    for ptm_type in PTM_TYPES:
        count_2017 = df_2017[f"has_{ptm_type}"].sum()
        count_2025 = df_2025[f"has_{ptm_type}"].sum()
        pct_2017 = (count_2017 / len(df_2017)) * 100
        pct_2025 = (count_2025 / len(df_2025)) * 100
        ptm_data.append(
            {
                "PTM Type": PTM_LABELS[ptm_type],
                "Count 2017": count_2017,
                "Count 2025": count_2025,
                "Pct 2017": pct_2017,
                "Pct 2025": pct_2025,
            }
        )

    df_ptm = pd.DataFrame(ptm_data).sort_values("Count 2025", ascending=False)

    # Stacked bars
    values_base = np.minimum(df_ptm["Count 2017"], df_ptm["Count 2025"])
    values_increase = np.maximum(0, df_ptm["Count 2025"] - df_ptm["Count 2017"])
    y_pos = np.arange(len(df_ptm))

    ax_a.barh(y_pos, values_base, color=COLOR_2017, edgecolor="black")
    ax_a.barh(
        y_pos, values_increase, left=values_base, color=COLOR_2025, edgecolor="black"
    )

    # Formatting
    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(df_ptm["PTM Type"])
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Number of Proteins", fontsize=14)
    ax_a.set_title(
        "PTM Type Frequency (2017 vs 2025)", fontsize=14, fontweight="bold", pad=10
    )
    # Increase the fontsize of both tick labels
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
        text_x = row["Count 2025"] + ax_a.get_xlim()[1] * 0.01
        annotation = f"{int(row['Count 2017'])} ({row['Pct 2017']:.1f}%) → {int(row['Count 2025'])} ({row['Pct 2025']:.1f}%)"
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
        Patch(facecolor=COLOR_2017, label="2017", edgecolor="black"),
        Patch(facecolor=COLOR_2025, label="2025", edgecolor="black"),
    ]
    ax_a.legend(handles=legend_elements, loc="lower right", fontsize=18)

    # ========== PLOT B: PTM Count Distributions ==========
    gs_b = gs[1].subgridspec(1, 3, wspace=0.15)
    ptm_types_to_plot = ["DISULFID", "MOD_RES", "CARBOHYD"]

    # First, collect all data to determine shared y-axis range
    all_hist_data = []
    for ptm_type in ptm_types_to_plot:
        counts_2017 = df_2017[df_2017[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]
        counts_2025 = df_2025[df_2025[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]

        num_bins = 10
        hist_2017 = np.array(
            [(counts_2017 == i).sum() for i in range(1, num_bins)]
            + [(counts_2017 >= num_bins).sum()]
        )
        hist_2025 = np.array(
            [(counts_2025 == i).sum() for i in range(1, num_bins)]
            + [(counts_2025 >= num_bins).sum()]
        )
        all_hist_data.append((hist_2017, hist_2025))

    # Create subplots with shared y-axis
    axes_b = []
    for idx, ptm_type in enumerate(ptm_types_to_plot):
        ax = fig.add_subplot(gs_b[idx])
        axes_b.append(ax)

        hist_2017, hist_2025 = all_hist_data[idx]

        # Plot overlapping bars
        x_pos = np.arange(1, 11)
        ax.bar(
            x_pos,
            hist_2025,
            width=1.0,
            color=COLOR_2025,
            edgecolor="black",
            align="edge",
        )
        ax.bar(
            x_pos,
            hist_2017,
            width=1.0,
            color=COLOR_2017,
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
    max_y = max(
        max(hist_2017.max(), hist_2025.max()) for hist_2017, hist_2025 in all_hist_data
    )
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
        counts_2017 = df_2017[df_2017[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]
        counts_2025 = df_2025[df_2025[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]

        for i in range(1, 10):
            dist_data.append(
                {
                    "PTM_Type": PTM_LABELS[ptm_type],
                    "PTM_Count": str(i),
                    "Proteins_2017": (counts_2017 == i).sum(),
                    "Proteins_2025": (counts_2025 == i).sum(),
                }
            )
        dist_data.append(
            {
                "PTM_Type": PTM_LABELS[ptm_type],
                "PTM_Count": "10+",
                "Proteins_2017": (counts_2017 >= 10).sum(),
                "Proteins_2025": (counts_2025 >= 10).sum(),
            }
        )

    pd.DataFrame(dist_data).to_csv(
        output_dir / "ptm_count_distribution_data.csv", index=False
    )


def generate_summary_table(df_2017, df_2025, output_dir):
    """Generate summary statistics table"""
    table_data = [
        ["Metric", "2017", "2025", "Change"],
        [
            "Total Proteins",
            f"{len(df_2017):,}",
            f"{len(df_2025):,}",
            f"+{len(df_2025) - len(df_2017):,}",
        ],
        [
            "Proteins with PTMs",
            f"{df_2017['has_any_ptm'].sum():,}",
            f"{df_2025['has_any_ptm'].sum():,}",
            f"+{df_2025['has_any_ptm'].sum() - df_2017['has_any_ptm'].sum():,}",
        ],
        [
            "PTM Coverage (%)",
            f"{(df_2017['has_any_ptm'].sum() / len(df_2017) * 100):.1f}%",
            f"{(df_2025['has_any_ptm'].sum() / len(df_2025) * 100):.1f}%",
            f"{(df_2025['has_any_ptm'].sum() / len(df_2025) - df_2017['has_any_ptm'].sum() / len(df_2017)) * 100:+.1f}%",
        ],
        [
            "Mean PTMs/Protein",
            f"{df_2017[df_2017['has_any_ptm']]['total_ptm_count'].mean():.2f}",
            f"{df_2025[df_2025['has_any_ptm']]['total_ptm_count'].mean():.2f}",
            f"{df_2025[df_2025['has_any_ptm']]['total_ptm_count'].mean() - df_2017[df_2017['has_any_ptm']]['total_ptm_count'].mean():+.2f}",
        ],
        ["", "", "", ""],
        ["PTM Type", "2017", "2025", "Change"],
    ]

    for ptm_type in PTM_TYPES:
        count_2017 = df_2017[f"has_{ptm_type}"].sum()
        count_2025 = df_2025[f"has_{ptm_type}"].sum()
        pct_2017 = (count_2017 / len(df_2017)) * 100
        pct_2025 = (count_2025 / len(df_2025)) * 100
        table_data.append(
            [
                PTM_LABELS[ptm_type],
                f"{count_2017:,} ({pct_2017:.1f}%)",
                f"{count_2025:,} ({pct_2025:.1f}%)",
                f"+{count_2025 - count_2017:,}",
            ]
        )

    fig, ax = plt.subplots(figsize=(12, 8))
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

    plt.title("PTM Summary Statistics", fontsize=16, fontweight="bold", pad=20)
    plt.savefig(output_dir / "ptm_summary_table.png", dpi=300, bbox_inches="tight")
    plt.close()

    pd.DataFrame(table_data[1:], columns=table_data[0]).to_csv(
        output_dir / "ptm_summary_statistics.csv", index=False
    )


def main():
    """Main analysis pipeline"""
    print("PTM Analysis")
    print("=" * 40)

    # Setup paths
    data_dir = Path("data/interim")
    output_dir = Path("figures/ptm")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("Loading data...")
    df_2017 = load_data(data_dir / "toxprot_2017.tsv")
    df_2025 = load_data(data_dir / "toxprot_2025.tsv")

    # Generate outputs
    print("Generating figures...")
    create_combined_figure(df_2017, df_2025, output_dir)
    generate_summary_table(df_2017, df_2025, output_dir)

    print(f"✓ Complete. Output: {output_dir.absolute()}")
    print("=" * 40)


if __name__ == "__main__":
    main()
