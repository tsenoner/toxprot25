#!/usr/bin/env python3
"""
PTM Analysis for ToxProt Dataset

Analyzes and visualizes post-translational modifications (PTMs) from the
PTM Summary column across ToxProt datasets at three time points.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Years available in the dataset
YEARS = [str(y) for y in range(2005, 2026)]

# Default comparison years (10-year intervals)
DEFAULT_YEARS = [2005, 2015, 2025]

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

# Colors for three time points
COLORS = {
    0: plt.cm.tab10(0),  # Blue for earliest
    1: plt.cm.tab10(2),  # Green for middle
    2: plt.cm.tab10(1),  # Orange for latest
}


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


def create_combined_figure(datasets: dict, output_dir: Path, top_n: int = 10):
    """Create combined figure with Plot A (type comparison) and Plot B (distributions).

    Args:
        datasets: Dictionary mapping year (int) to DataFrame
        output_dir: Directory to save output files
        top_n: Number of top PTM types to show in bar chart (default 10 = all)
    """
    years = sorted(datasets.keys())
    fig = plt.figure(figsize=(18, 16))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.4], hspace=0.25)

    # ========== PLOT A: PTM Type Frequency (stacked bars) ==========
    ax_a = fig.add_subplot(gs[0])

    # Prepare data for top_n PTM types
    ptm_data = []
    for ptm_type in PTM_TYPES[:top_n]:
        row_data = {"PTM Type": PTM_LABELS[ptm_type]}
        for year in years:
            df = datasets[year]
            count = df[f"has_{ptm_type}"].sum()
            pct = (count / len(df)) * 100
            row_data[f"Count {year}"] = count
            row_data[f"Pct {year}"] = pct
        ptm_data.append(row_data)

    df_ptm = pd.DataFrame(ptm_data).sort_values(f"Count {years[-1]}", ascending=False)

    # Stacked bars: plot in reverse order (latest year first/back, earliest year last/front)
    y_pos = np.arange(len(df_ptm))

    # Plot years from latest to earliest so earlier years appear on top
    for i, year in enumerate(reversed(years)):
        color_idx = len(years) - 1 - i  # Reverse color index
        ax_a.barh(
            y_pos,
            df_ptm[f"Count {year}"],
            color=COLORS[color_idx],
            edgecolor="black",
            label=str(year),
        )

    # Formatting
    ax_a.set_yticks(y_pos)
    ax_a.set_yticklabels(df_ptm["PTM Type"])
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Number of Proteins", fontsize=14)
    years_str = ", ".join(str(y) for y in years)
    ax_a.set_title(f"PTM Type Frequency ({years_str})", fontsize=14, fontweight="bold", pad=10)
    ax_a.tick_params(axis="both", labelsize=12)

    # Add panel label A
    ax_a.text(
        -0.1,
        1.02,
        "A",
        transform=ax_a.transAxes,
        fontsize=28,
        fontweight="bold",
        va="bottom",
        ha="left",
    )
    ax_a.set_axisbelow(True)
    ax_a.grid(axis="x", ls="--", lw=2, alpha=0.5)

    # Annotations - show counts for each year
    for row_idx, (_, row) in enumerate(df_ptm.iterrows()):
        counts_str = " → ".join(f"{int(row[f'Count {y}'])}" for y in years)
        text_x = max(row[f"Count {y}"] for y in years) + ax_a.get_xlim()[1] * 0.01
        ax_a.text(
            text_x,
            y_pos[row_idx],
            counts_str,
            va="center",
            ha="left",
            fontsize=10,
            fontweight="bold",
        )

    ax_a.set_xlim(ax_a.get_xlim()[0], ax_a.get_xlim()[1] * 1.12)
    # Reverse legend order so it matches visual (earliest on top)
    handles, labels = ax_a.get_legend_handles_labels()
    ax_a.legend(handles[::-1], labels[::-1], loc="lower right", fontsize=14)

    # ========== PLOT B: PTM Count Distributions (2x3 grid) ==========
    gs_b = gs[1].subgridspec(2, 3, wspace=0.15, hspace=0.35)
    # Top 6 PTM types for distribution plots
    ptm_types_to_plot = PTM_TYPES[:6]

    # First, collect all data to determine shared y-axis range
    all_hist_data = []
    num_bins = 10
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

    # Create subplots with shared y-axis
    axes_b = []
    for idx, ptm_type in enumerate(ptm_types_to_plot):
        row, col = idx // 3, idx % 3
        ax = fig.add_subplot(gs_b[row, col])
        axes_b.append(ax)

        year_hists = all_hist_data[idx]

        # Plot stacked bars: latest year first (back), earliest year last (front)
        x_pos = np.arange(1, 11)
        for i, year in enumerate(reversed(years)):
            color_idx = len(years) - 1 - i
            ax.bar(
                x_pos,
                year_hists[year],
                width=0.8,
                color=COLORS[color_idx],
                edgecolor="black",
                label=str(year) if idx == 0 else None,
            )

        # Formatting
        ax.set_title(PTM_LABELS[ptm_type], fontsize=12, fontweight="bold", pad=8)
        ax.set_xlim(0.5, 10.5)
        ax.set_xticks(x_pos)
        ax.set_xticklabels([str(i) for i in range(1, 10)] + ["10+"], fontsize=9)
        ax.tick_params(axis="both", labelsize=10)

        # Only show x-label on bottom row
        if row == 1:
            ax.set_xlabel("Number of PTMs", fontsize=12)

        # Only show y-label on left column
        if col == 0:
            ax.set_ylabel("Number of Proteins", fontsize=12)

        if idx == 0:
            # Add panel label B to the first subplot
            ax.text(
                -0.35,
                1.08,
                "B",
                transform=ax.transAxes,
                fontsize=28,
                fontweight="bold",
                va="bottom",
                ha="left",
            )
            # Reverse legend order so it matches visual (earliest on top)
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], loc="upper right", fontsize=9)

        ax.set_axisbelow(True)
        ax.grid(axis="y", ls="--", lw=1, alpha=0.5)

    # Set shared y-axis limits
    max_y = max(max(h.max() for h in year_hists.values()) for year_hists in all_hist_data)
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
        for i in range(1, 10):
            row_data = {"PTM_Type": PTM_LABELS[ptm_type], "PTM_Count": str(i)}
            for year in years:
                df = datasets[year]
                counts = df[df[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]
                row_data[f"Proteins_{year}"] = (counts == i).sum()
            dist_data.append(row_data)
        # 10+ bin
        row_data = {"PTM_Type": PTM_LABELS[ptm_type], "PTM_Count": "10+"}
        for year in years:
            df = datasets[year]
            counts = df[df[f"count_{ptm_type}"] > 0][f"count_{ptm_type}"]
            row_data[f"Proteins_{year}"] = (counts >= 10).sum()
        dist_data.append(row_data)

    pd.DataFrame(dist_data).to_csv(output_dir / "ptm_count_distribution_data.csv", index=False)


def generate_summary_table(datasets: dict, output_dir: Path):
    """Generate summary statistics table.

    Args:
        datasets: Dictionary mapping year (int) to DataFrame
        output_dir: Directory to save output files
    """
    years = sorted(datasets.keys())
    num_cols = len(years) + 2  # Metric + years + Total Change

    # Header row
    header = ["Metric"] + [str(y) for y in years] + ["Total Change"]
    table_data = [header]

    # Helper to compute change
    def fmt_change(new_val, old_val):
        diff = new_val - old_val
        return f"+{diff:,}" if diff >= 0 else f"{diff:,}"

    # Total proteins
    row = ["Total Proteins"]
    for year in years:
        row.append(f"{len(datasets[year]):,}")
    row.append(fmt_change(len(datasets[years[-1]]), len(datasets[years[0]])))
    table_data.append(row)

    # Proteins with PTMs
    row = ["Proteins with PTMs"]
    for year in years:
        row.append(f"{datasets[year]['has_any_ptm'].sum():,}")
    row.append(
        fmt_change(
            datasets[years[-1]]["has_any_ptm"].sum(), datasets[years[0]]["has_any_ptm"].sum()
        )
    )
    table_data.append(row)

    # PTM Coverage
    row = ["PTM Coverage (%)"]
    for year in years:
        df = datasets[year]
        pct = df["has_any_ptm"].sum() / len(df) * 100
        row.append(f"{pct:.2f}%")
    pct_first = datasets[years[0]]["has_any_ptm"].sum() / len(datasets[years[0]]) * 100
    pct_last = datasets[years[-1]]["has_any_ptm"].sum() / len(datasets[years[-1]]) * 100
    diff = pct_last - pct_first
    row.append(f"{diff:+.2f}%")
    table_data.append(row)

    # Mean PTMs/Protein
    row = ["Mean PTMs/Protein"]
    for year in years:
        df = datasets[year]
        mean_ptm = df[df["has_any_ptm"]]["total_ptm_count"].mean()
        row.append(f"{mean_ptm:.2f}")
    mean_first = datasets[years[0]][datasets[years[0]]["has_any_ptm"]]["total_ptm_count"].mean()
    mean_last = datasets[years[-1]][datasets[years[-1]]["has_any_ptm"]]["total_ptm_count"].mean()
    diff = mean_last - mean_first
    row.append(f"{diff:+.2f}")
    table_data.append(row)

    # Separator row
    table_data.append([""] * num_cols)

    # PTM Type header
    ptm_header = ["PTM Type"] + [str(y) for y in years] + ["Total Change"]
    table_data.append(ptm_header)

    # PTM types
    for ptm_type in PTM_TYPES:
        row = [PTM_LABELS[ptm_type]]
        for year in years:
            df = datasets[year]
            count = df[f"has_{ptm_type}"].sum()
            pct = (count / len(df)) * 100
            row.append(f"{count:,} ({pct:.1f}%)")
        count_first = datasets[years[0]][f"has_{ptm_type}"].sum()
        count_last = datasets[years[-1]][f"has_{ptm_type}"].sum()
        row.append(fmt_change(count_last, count_first))
        table_data.append(row)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.axis("tight")
    ax.axis("off")

    col_widths = [0.35] + [0.18] * len(years) + [0.15]
    table = ax.table(
        cellText=table_data,
        cellLoc="left",
        loc="center",
        colWidths=col_widths,
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2.5)

    # Style header rows
    header_rows = [0, 6]
    for i in header_rows:
        for j in range(num_cols):
            table[(i, j)].set_facecolor("#34495e")
            table[(i, j)].set_text_props(weight="bold", color="white")

    # Style separator row
    for j in range(num_cols):
        table[(5, j)].set_facecolor("#ecf0f1")
        table[(5, j)].set_alpha(0.3)

    # Alternating row colors
    for i in range(1, len(table_data)):
        if i not in [5, 6]:
            for j in range(num_cols):
                if i % 2 == 0:
                    table[(i, j)].set_facecolor("#f8f9fa")

    years_str = ", ".join(str(y) for y in years)
    plt.title(
        f"PTM Summary Statistics ({years_str})",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )
    plt.savefig(output_dir / "ptm_summary_table.png", dpi=300, bbox_inches="tight")
    plt.close()

    pd.DataFrame(table_data[1:], columns=table_data[0]).to_csv(
        output_dir / "ptm_summary_statistics.csv", index=False
    )


def plot_ptm_trends(datasets: dict, output_dir: Path, top_n: int = 10):
    """Plot PTM type trends over all years.

    Args:
        datasets: Dictionary mapping year (int) to DataFrame
        output_dir: Directory to save output files
        top_n: Number of top PTM types to show
    """
    years = sorted(datasets.keys())

    # Collect data for each PTM type across years
    trend_data = {ptm: [] for ptm in PTM_TYPES[:top_n]}
    for year in years:
        df = datasets[year]
        for ptm_type in PTM_TYPES[:top_n]:
            count = df[f"has_{ptm_type}"].sum()
            trend_data[ptm_type].append(count)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))

    # Plot each PTM type
    colors = plt.cm.tab10.colors
    for i, ptm_type in enumerate(PTM_TYPES[:top_n]):
        ax.plot(
            years,
            trend_data[ptm_type],
            marker="o",
            linewidth=2,
            markersize=6,
            color=colors[i],
            label=PTM_LABELS[ptm_type],
        )

    ax.set_xlabel("Year", fontsize=14)
    ax.set_ylabel("Number of Proteins", fontsize=14)
    ax.set_title("PTM Type Trends (2005-2025)", fontsize=14, fontweight="bold")
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xticks(years)
    ax.set_xticklabels([str(y) for y in years], rotation=45, ha="right")
    ax.legend(loc="upper left", fontsize=10)
    ax.set_axisbelow(True)
    ax.grid(True, ls="--", alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_dir / "ptm_trends.png", dpi=300, bbox_inches="tight")
    plt.savefig(output_dir / "ptm_trends.pdf", bbox_inches="tight")
    plt.close()

    # Save trend data
    trend_df = pd.DataFrame({"Year": years})
    for ptm_type in PTM_TYPES[:top_n]:
        trend_df[PTM_LABELS[ptm_type]] = trend_data[ptm_type]
    trend_df.to_csv(output_dir / "ptm_trends_data.csv", index=False)


def main():
    """Main analysis pipeline."""
    print("PTM Analysis")
    print("=" * 40)

    # Setup paths
    data_dir = Path("data/processed/toxprot")
    output_dir = Path("figures/ptm")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data for default years
    print("Loading data...")
    datasets = {}
    for year in DEFAULT_YEARS:
        datasets[year] = load_data(data_dir / f"toxprot_{year}.csv")
        print(f"  {year}: {len(datasets[year]):,} entries")

    # Generate outputs
    print("Generating figures...")
    create_combined_figure(datasets, output_dir)
    generate_summary_table(datasets, output_dir)

    print(f"Complete. Output: {output_dir.absolute()}")
    print("=" * 40)


if __name__ == "__main__":
    main()
