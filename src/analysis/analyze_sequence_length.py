#!/usr/bin/env python3
"""
Analyze and visualize sequence length distributions in ToxProt datasets.

This script generates histograms comparing sequence lengths across
multiple ToxProt time points (2005, 2015, 2025), for three sequence
types: full precursor, mature (signal peptide removed), and active
peptide (signal peptide + propeptide removed).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..config import COMPARISON_YEARS, DATA_DIR, FIGURES_DIR
from .colors import YEAR_COLORS

# Font sizes for consistent styling
AXIS_LABEL_FONTSIZE = 16
TICK_LABEL_FONTSIZE = 12
LEGEND_FONTSIZE = 16


def plot_sequence_length_histogram(
    datasets: dict[int, pd.DataFrame],
    output_path: Path,
    length_column: str = "Length",
    xlabel: str = "Sequence Length (amino acids)",
) -> None:
    """
    Create overlaid histogram comparing sequence lengths across years.

    Args:
        datasets: Dictionary mapping year (int) to DataFrame.
        output_path: Path to save the figure.
        length_column: Column name to use for lengths.
        xlabel: Label for the x-axis.
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
        lengths = df[length_column].dropna()

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
    ax.set_xlabel(xlabel, fontsize=AXIS_LABEL_FONTSIZE)
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


def _filter_complete(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to complete sequences (no fragments)."""
    return df[df["Fragment"] != "fragment"]


def _filter_complete_with_sp(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to complete sequences with signal peptide annotation."""
    mask = df["Fragment"] != "fragment"
    if "Signal peptide" in df.columns:
        mask = mask & (df["Signal peptide"] == "Yes")
    return df[mask]


def _filter_complete_with_propep(df: pd.DataFrame) -> pd.DataFrame:
    """Filter to complete sequences with propeptide annotation."""
    mask = df["Fragment"] != "fragment"
    if "Propeptide" in df.columns:
        mask = mask & (df["Propeptide"] == "Yes")
    return df[mask]


def generate_all_length_figures(
    datasets: dict[int, pd.DataFrame], output_dir: Path
) -> list[Path]:
    """Generate three separate sequence length distribution figures.

    Args:
        datasets: Dictionary mapping year to DataFrame with Length,
            Mature_length, and Active_length columns.
        output_dir: Directory to save figures.

    Returns:
        List of paths to generated figures.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    generated = []

    # Figure 1: Full-length sequences (complete, no fragments)
    full_datasets = {
        year: _filter_complete(df)
        for year, df in datasets.items()
    }
    full_datasets = {y: df for y, df in full_datasets.items() if len(df) > 0}
    if len(full_datasets) >= 2:
        path = output_dir / "sequence_length_full.png"
        plot_sequence_length_histogram(
            full_datasets, path,
            length_column="Length",
            xlabel="Full-length Sequence Length (amino acids)",
        )
        generated.append(path)

    # Figure 2: Full-length sequences with SP annotation (confirmed prepropeptides)
    full_sp_datasets = {
        year: _filter_complete_with_sp(df)
        for year, df in datasets.items()
    }
    full_sp_datasets = {y: df for y, df in full_sp_datasets.items() if len(df) > 0}
    if len(full_sp_datasets) >= 2:
        path = output_dir / "sequence_length_full_sp.png"
        plot_sequence_length_histogram(
            full_sp_datasets, path,
            length_column="Length",
            xlabel="Full-length Sequence Length (amino acids)",
        )
        generated.append(path)

    # Figure 3: Mature lengths — SP removed where annotated (complete, no fragments)
    mature_datasets = {
        year: _filter_complete(df)
        for year, df in datasets.items()
    }
    mature_datasets = {y: df for y, df in mature_datasets.items() if len(df) > 0}
    has_mature = all("Mature_length" in df.columns for df in mature_datasets.values())
    if len(mature_datasets) >= 2 and has_mature:
        path = output_dir / "sequence_length_mature.png"
        plot_sequence_length_histogram(
            mature_datasets, path,
            length_column="Mature_length",
            xlabel="Mature Sequence Length (amino acids)",
        )
        generated.append(path)

    # Figure 4: Active peptide lengths — SP + propeptide removed (only entries with propeptide)
    active_datasets = {
        year: _filter_complete_with_propep(df)
        for year, df in datasets.items()
    }
    active_datasets = {y: df for y, df in active_datasets.items() if len(df) > 0}
    has_active = all("Active_length" in df.columns for df in active_datasets.values())
    if len(active_datasets) >= 2 and has_active:
        path = output_dir / "sequence_length_active.png"
        plot_sequence_length_histogram(
            active_datasets, path,
            length_column="Active_length",
            xlabel="Active Peptide Length (amino acids)",
        )
        generated.append(path)

    # Figure 5: Mature peptide lengths — SP + propeptide removed where annotated,
    # entries without annotations assumed already mature (all complete sequences)
    mature_pep_datasets = {
        year: _filter_complete(df)
        for year, df in datasets.items()
    }
    mature_pep_datasets = {y: df for y, df in mature_pep_datasets.items() if len(df) > 0}
    has_active = all("Active_length" in df.columns for df in mature_pep_datasets.values())
    if len(mature_pep_datasets) >= 2 and has_active:
        path = output_dir / "sequence_length_mature_peptide.png"
        plot_sequence_length_histogram(
            mature_pep_datasets, path,
            length_column="Active_length",
            xlabel="Mature Peptide Length (amino acids)",
        )
        generated.append(path)

    # Figure 6: Same as Figure 5 but including fragments
    all_datasets = {y: df for y, df in datasets.items() if len(df) > 0}
    has_active = all("Active_length" in df.columns for df in all_datasets.values())
    if len(all_datasets) >= 2 and has_active:
        path = output_dir / "sequence_length_mature_peptide_with_fragments.png"
        plot_sequence_length_histogram(
            all_datasets, path,
            length_column="Active_length",
            xlabel="Mature Peptide Length (amino acids)",
        )
        generated.append(path)

    return generated


def main():
    """Main function to run sequence length analysis."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyze and visualize sequence length distributions."
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
        help="Directory to save output figures (default: figures)",
    )

    args = parser.parse_args()

    # Load datasets for comparison years (default: venom_tissue definition)
    from .helpers import filter_by_definition

    years = COMPARISON_YEARS
    print("Loading ToxProt datasets...")
    datasets = {}
    for year in years:
        filepath = args.data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            df = filter_by_definition(pd.read_csv(filepath), "venom_tissue")
            datasets[year] = df
            print(f"  → {year}: {len(df):,} entries")
        else:
            print(f"  → {year}: not found, skipping")

    if len(datasets) < 2:
        print("Error: Need at least 2 datasets")
        return

    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Generate all figures
    print("\nGenerating sequence length histograms...")
    paths = generate_all_length_figures(datasets, args.output_dir)
    for p in paths:
        print(f"  Saved: {p}")

    print("\nDone!")


if __name__ == "__main__":
    main()
