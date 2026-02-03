#!/usr/bin/env python3
"""
Analyze and visualize sequence length distributions in ToxProt datasets.

This script generates a histogram comparing sequence lengths across
multiple ToxProt time points (2005, 2015, 2025).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .colors import YEAR_COLORS

# Font sizes for consistent styling
AXIS_LABEL_FONTSIZE = 16
TICK_LABEL_FONTSIZE = 12
LEGEND_FONTSIZE = 16


def plot_sequence_length_histogram(datasets: dict[int, pd.DataFrame], output_path: Path) -> None:
    """
    Create overlaid histogram comparing sequence lengths across years.

    Args:
        datasets: Dictionary mapping year (int) to DataFrame with 'Length' column.
        output_path: Path to save the figure.
    """
    # Define bin edges: 25 AA bins up to 300, then 301+
    bin_edges_up_to_300 = np.arange(1, 302, 25)
    last_bin_edge = bin_edges_up_to_300[-1]  # 301
    all_hist_bins = np.append(bin_edges_up_to_300, last_bin_edge + 25)

    # Create bin labels
    bin_labels = [
        f"{start}-{end - 1}"
        for start, end in zip(all_hist_bins[:-2], all_hist_bins[1:-1], strict=False)
    ]
    bin_labels.append(f"{all_hist_bins[-2]}+")  # "301+"

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))

    # Plot in order: latest year (back) to earliest year (front)
    for year in sorted(datasets.keys(), reverse=True):
        df = datasets[year]
        lengths = df["Length"].dropna()

        # Map all lengths >= 301 to 301 for binning
        lengths_binned = lengths.copy()
        lengths_binned[lengths_binned >= last_bin_edge] = last_bin_edge

        ax.hist(
            lengths_binned,
            bins=all_hist_bins,
            color=YEAR_COLORS[year],
            label=str(year),
            alpha=1,
            edgecolor="black",
        )

    # Styling
    ax.set_xlabel("Sequence Length (amino acids)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Count", fontsize=AXIS_LABEL_FONTSIZE)

    # Set x-axis ticks to bin centers
    tick_centers = (all_hist_bins[:-1] + all_hist_bins[1:]) / 2.0
    ax.set_xticks(tick_centers)
    ax.set_xticklabels(bin_labels, fontsize=TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="y", labelsize=TICK_LABEL_FONTSIZE)

    ax.set_axisbelow(True)
    ax.grid(axis="y", linestyle="--", linewidth=2, alpha=0.5)
    ax.legend(fontsize=LEGEND_FONTSIZE)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def main():
    """Main function to run sequence length analysis."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze and visualize sequence length distributions."
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

    args = parser.parse_args()

    # Load datasets for years 2005, 2015, 2025
    years = [2005, 2015, 2025]
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

    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate histogram
    output_path = args.output_dir / "sequence_length_distribution.png"
    print("\nGenerating sequence length histogram...")
    plot_sequence_length_histogram(datasets, output_path)
    print(f"  Saved: {output_path}")

    print("\nDone!")


if __name__ == "__main__":
    main()
