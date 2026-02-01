#!/usr/bin/env python3
"""
Analyze and visualize protein family distributions in ToxProt datasets.

This script generates:
1. Protein families stacked bar chart (figures/protein_families/)
2. Summary statistics table (figures/)
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Constants
FAMILY_NAME_MAP_2017_TO_2025 = {
    "Huwentoxin-1 family": "Neurotoxin 10 (Hwtx-1) family",
    "Spider toxin CSTX superfamily": "Neurotoxin 19 (CSTX) family",
}

# Font sizes for consistent styling
TITLE_FONTSIZE = 20
AXIS_LABEL_FONTSIZE = 16
TICK_LABEL_FONTSIZE = 12
LEGEND_FONTSIZE = 14


def get_top_families_per_dataset(df, column, top_n=15):
    """
    Get top N families plus 'Other' and 'NaN' categories.

    Args:
        df: DataFrame containing protein data
        column: Column name to analyze
        top_n: Number of top families to include

    Returns:
        Series with top families, 'Other', and 'NaN' counts
    """
    counts = df[column].value_counts()
    top_counts = counts.nlargest(top_n)
    other_count = counts.iloc[top_n:].sum()
    nan_count = df[column].isna().sum()

    result_series = top_counts.copy()
    if other_count > 0:
        result_series["Other"] = other_count
    if nan_count > 0:
        result_series["NaN"] = nan_count
    return result_series


def _create_family_order(top_17, top_25):
    """Create ordered list of families merging both years."""
    top_17_sorted = top_17.sort_values(ascending=False)
    top_25_sorted = top_25.sort_values(ascending=False)

    ordered_families = []
    seen_canonical = set()

    # Add 2017 families (mapped to 2025 names)
    for fam_name in top_17_sorted.index:
        canonical = FAMILY_NAME_MAP_2017_TO_2025.get(fam_name, fam_name)
        if canonical not in seen_canonical:
            ordered_families.append(canonical)
            seen_canonical.add(canonical)

    # Add any new 2025 families
    for fam_name in top_25_sorted.index:
        if fam_name not in seen_canonical:
            ordered_families.append(fam_name)
            seen_canonical.add(fam_name)

    return ordered_families


def _create_color_map(families):
    """Create color mapping for protein families."""
    cmap = plt.colormaps["tab20"]
    color_map = {
        "Other": cmap(14),
        "NaN": cmap(15),
    }

    # Colors for regular families
    palette = [c for i, c in enumerate(cmap.colors) if i not in [14, 15]]
    regular_families = [f for f in families if f not in ("Other", "NaN")]

    for i, family_name in enumerate(regular_families):
        color_map[family_name] = palette[i % len(palette)]

    return color_map


def plot_stacked_bar_protein_families(df_2017, df_2025, output_path, top_n=15):
    """
    Create 100% stacked bar chart of top protein families.

    Args:
        df_2017: DataFrame with 2017 data
        df_2025: DataFrame with 2025 data
        output_path: Path to save the figure
        top_n: Number of top families to display
    """
    # Get top families for each year
    top_17 = get_top_families_per_dataset(df_2017, "Protein families", top_n=top_n)
    top_25 = get_top_families_per_dataset(df_2025, "Protein families", top_n=top_n)

    # Create unified family order
    all_families = _create_family_order(top_17, top_25)

    # Build data frame with counts and percentages
    df_counts = pd.DataFrame(index=all_families)
    top_17_renamed = top_17.rename(index=FAMILY_NAME_MAP_2017_TO_2025)
    df_counts["2017"] = top_17_renamed.reindex(all_families).fillna(0)
    df_counts["2025"] = top_25.reindex(all_families).fillna(0)

    df_percentages = df_counts.div(df_counts.sum(axis=0), axis=1) * 100

    # Create color mapping
    color_map = _create_color_map(all_families)

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(13, 9))

    # Helper function to create stacked bars
    def create_stacked_bars(year_label, year_col):
        bars_info = []
        bottom = 0
        for family in all_families:
            pct = df_percentages.at[family, year_col]
            cnt = df_counts.at[family, year_col]
            if cnt > 0:
                bar = ax.bar(
                    year_label,
                    pct,
                    bottom=bottom,
                    color=color_map.get(family, "gray"),
                )[0]
                bars_info.append(
                    {
                        "bar": bar,
                        "family": family,
                        "count": cnt,
                        "percentage": pct,
                    }
                )
                bottom += pct
        return bars_info

    # Create bars for both years
    bars_2017 = create_stacked_bars("ToxProt 2017", "2017")
    bars_2025 = create_stacked_bars("ToxProt 2025", "2025")

    # Add labels to bars (only if percentage > 1%)
    for bar_info in bars_2017 + bars_2025:
        if bar_info["percentage"] > 1.0:
            bar = bar_info["bar"]
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_y() + bar.get_height() / 2,
                f"{bar_info['family']}, {int(bar_info['count'])} ({bar_info['percentage']:.1f}%)",
                ha="center",
                va="center",
                fontsize=9,
            )

    # Styling
    ax.set_ylabel("Percentage", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_title(
        f"Top {top_n} Protein Families Distribution",
        fontsize=TITLE_FONTSIZE,
    )
    ax.set_ylim(0, 100)
    ax.tick_params(axis="x", labelsize=TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="y", labelsize=TICK_LABEL_FONTSIZE)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def _calculate_statistics(df):
    """Calculate summary statistics for a dataset."""
    total_count = len(df)

    # Protein families
    unique_families = len(df["Protein families"].dropna().unique())
    na_count = df["Protein families"].isna().sum()
    na_percentage = (na_count / total_count) * 100 if total_count > 0 else 0

    # Fragments
    df_copy = df.copy()
    df_copy["Fragment"] = df_copy["Fragment"].astype(str)
    fragment_count = df_copy["Fragment"].str.contains("fragment", case=False, na=False).sum()
    fragment_percentage = (fragment_count / total_count) * 100 if total_count > 0 else 0

    # PTM annotations
    ptm_count = df["PTM Summary"].notna().sum()
    ptm_percentage = (ptm_count / total_count) * 100 if total_count > 0 else 0

    # Toxic dose annotations
    toxic_count = df["Toxic dose"].notna().sum()
    toxic_percentage = (toxic_count / total_count) * 100 if total_count > 0 else 0

    # Taxonomy
    species_count = len(df["Species"].unique())
    order_count = len(df["Order"].unique())

    return {
        "total": total_count,
        "unique_families": unique_families,
        "na_count": na_count,
        "na_pct": na_percentage,
        "fragment_count": fragment_count,
        "fragment_pct": fragment_percentage,
        "ptm_count": ptm_count,
        "ptm_pct": ptm_percentage,
        "toxic_count": toxic_count,
        "toxic_pct": toxic_percentage,
        "species_count": species_count,
        "order_count": order_count,
    }


def generate_summary_table(datasets: dict[str, pd.DataFrame], output_path: Path):
    """
    Generate summary statistics table comparing all datasets.

    Args:
        datasets: Dictionary mapping year to DataFrame
        output_path: Path to save the figure
    """
    # Calculate statistics for all years
    all_stats = {year: _calculate_statistics(df) for year, df in datasets.items()}
    years = list(datasets.keys())

    # Create table data
    cell_text = [
        [all_stats[y]["total"] for y in years],
        [all_stats[y]["unique_families"] for y in years],
        [f"{all_stats[y]['na_count']} ({all_stats[y]['na_pct']:.1f}%)" for y in years],
        [f"{all_stats[y]['fragment_count']} ({all_stats[y]['fragment_pct']:.1f}%)" for y in years],
        [f"{all_stats[y]['ptm_count']} ({all_stats[y]['ptm_pct']:.1f}%)" for y in years],
        [f"{all_stats[y]['toxic_count']} ({all_stats[y]['toxic_pct']:.1f}%)" for y in years],
        [all_stats[y]["species_count"] for y in years],
        [all_stats[y]["order_count"] for y in years],
    ]

    row_headers = [
        "Total entries",
        "Unique protein families",
        "Missing protein family",
        "Fragment entries",
        "PTM annotations",
        "Toxic dose annotations",
        "Species count",
        "Order count",
    ]
    column_headers = years

    # Create figure and table
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.axis("tight")
    ax.axis("off")

    # Color scheme
    row_colors = plt.cm.BuPu(np.full(len(row_headers), 0.1))
    col_colors = plt.cm.BuPu(np.full(len(column_headers), 0.1))

    col_width = 0.22
    table = plt.table(
        cellText=cell_text,
        rowLabels=row_headers,
        rowColours=row_colors,
        rowLoc="center",
        colColours=col_colors,
        colLabels=column_headers,
        cellLoc="center",
        loc="center",
        colWidths=[col_width] * len(years),
    )

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(12)

    for cell in table.properties()["children"]:
        cell.set_height(0.12)

    # Make headers larger
    for key in table._cells:
        if key[0] == 0 or key[1] == -1:  # Column or row labels
            table._cells[key].set_fontsize(13)
            table._cells[key].set_text_props(weight="bold")

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    """Main function to analyze protein families."""
    parser = argparse.ArgumentParser(
        description="Analyze and visualize protein family distributions."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data/processed/toxprot"),
        help="Directory containing processed CSV files (default: data/processed/toxprot)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures"),
        help="Directory to save output figures (default: figures)",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Number of top protein families to display (default: 15)",
    )

    args = parser.parse_args()

    # Load all datasets
    years = ["2005", "2015", "2025"]
    print("Loading ToxProt datasets...")
    datasets = {}
    for year in years:
        filepath = args.data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)
            print(f"  → {year}: {len(datasets[year]):,} entries")
        else:
            print(f"  → {year}: not found, skipping")

    if len(datasets) < 2:
        print("Error: Need at least 2 datasets")
        return

    # Ensure output directories exist
    args.output_dir.mkdir(parents=True, exist_ok=True)
    protein_families_dir = args.output_dir / "protein_families"
    protein_families_dir.mkdir(parents=True, exist_ok=True)

    # Define output paths
    families_bar_path = protein_families_dir / "top_families_stacked_bar.png"
    summary_table_path = args.output_dir / "dataset_summary_statistics.png"

    # Generate visualizations
    print("\nGenerating visualizations...")

    # Use 2017 and 2025 for families chart
    if "2017" in datasets and "2025" in datasets:
        print(f"  → Top {args.top_n} protein families stacked bar chart...")
        plot_stacked_bar_protein_families(
            datasets["2017"], datasets["2025"], families_bar_path, top_n=args.top_n
        )
        print(f"    ✓ {families_bar_path}")

    print("  → Summary statistics table...")
    generate_summary_table(datasets, summary_table_path)
    print(f"    ✓ {summary_table_path}")

    print("\n" + "=" * 60)
    print("✓ All visualizations complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
