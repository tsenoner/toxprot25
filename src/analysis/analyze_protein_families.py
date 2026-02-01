#!/usr/bin/env python3
"""Analyze and visualize protein family distributions in ToxProt datasets.

Generates a 100% stacked bar chart comparing protein family distributions
between 2017 and 2025 ToxProt datasets.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
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


def main():
    """Main function for standalone execution."""
    parser = argparse.ArgumentParser(
        description="Generate protein family stacked bar chart."
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

    # Load 2017 and 2025 datasets
    print("Loading ToxProt datasets...")
    datasets = {}
    for year in ["2017", "2025"]:
        filepath = args.data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)
            print(f"  {year}: {len(datasets[year]):,} entries")
        else:
            print(f"Error: {filepath} not found")
            return

    if "2017" not in datasets or "2025" not in datasets:
        print("Error: Need both 2017 and 2025 datasets")
        return

    # Ensure output directory exists
    protein_families_dir = args.output_dir / "protein_families"
    protein_families_dir.mkdir(parents=True, exist_ok=True)

    output_path = protein_families_dir / "top_families_stacked_bar.png"
    plot_stacked_bar_protein_families(
        datasets["2017"], datasets["2025"], output_path, top_n=args.top_n
    )
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
