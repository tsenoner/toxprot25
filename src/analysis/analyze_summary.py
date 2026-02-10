#!/usr/bin/env python3
"""Generate summary statistics table for ToxProt datasets."""

import argparse
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import pandas as pd

from ..config import COMPARISON_YEARS, DATA_DIR, FIGURES_DIR
from .colors import YEAR_COLORS

# Transparency for column header colors
HEADER_ALPHA = 0.3


def _calculate_statistics(df: pd.DataFrame) -> dict:
    """Calculate summary statistics for a dataset.

    Args:
        df: DataFrame containing ToxProt data.

    Returns:
        Dictionary with calculated statistics.
    """
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


def generate_summary_table(datasets: dict[int, pd.DataFrame], output_path: Path):
    """Generate summary statistics table comparing all datasets.

    Args:
        datasets: Dictionary mapping year (int) to DataFrame.
        output_path: Path to save the figure.
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
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axis("off")

    # Color scheme - year colors with transparency for column headers
    col_colors = [mcolors.to_rgba(YEAR_COLORS[y], alpha=HEADER_ALPHA) for y in years]

    table = plt.table(
        cellText=cell_text,
        rowLabels=row_headers,
        rowLoc="center",
        colColours=col_colors,
        colLabels=column_headers,
        cellLoc="center",
        loc="center",
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

    plt.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.1)
    plt.close()


def main():
    """Main function for standalone execution."""
    parser = argparse.ArgumentParser(
        description="Generate summary statistics table for ToxProt datasets."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DATA_DIR,
        help="Directory containing processed CSV files (default: data/processed/toxprot)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=FIGURES_DIR,
        help="Directory to save output figure (default: figures)",
    )

    args = parser.parse_args()

    # Load datasets (default: venom_tissue definition)
    from .helpers import filter_by_definition

    years = COMPARISON_YEARS
    print("Loading ToxProt datasets...")
    datasets = {}
    for year in years:
        filepath = args.data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            df = filter_by_definition(pd.read_csv(filepath), "venom_tissue")
            datasets[year] = df
            print(f"  {year}: {len(df):,} entries")

    if len(datasets) < 2:
        print("Error: Need at least 2 datasets")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "dataset_summary_statistics.png"

    generate_summary_table(datasets, output_path)
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
