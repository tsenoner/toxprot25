#!/usr/bin/env python3
"""Analyze and visualize source tissue evolution in ToxProt datasets.

This module creates an alluvial diagram showing how expression/source tissues
have changed over time (2005 -> 2015 -> 2025).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath

from ..config import COMPARISON_YEARS, FIGURES_DIR
from .colors import CATEGORICAL_PALETTE, OTHER_COLOR
from .helpers import load_datasets

# --- Configuration ---
YEARS = COMPARISON_YEARS


def explode_source_tissues(df: pd.DataFrame) -> pd.Series:
    """Split 'Source tissues' column and count individual tissues.

    Args:
        df: DataFrame with 'Source tissues' column containing semicolon-separated values.

    Returns:
        Series with tissue counts.
    """
    # Get non-empty source tissues
    tissues = df["Source tissues"].dropna()
    tissues = tissues[tissues != ""]

    # Split by "; " and explode to individual rows
    exploded = tissues.str.split("; ").explode()

    # Strip whitespace and filter empty
    exploded = exploded.str.strip()
    exploded = exploded[exploded != ""]

    return exploded.value_counts()


def get_top_tissues(
    datasets: dict[int, pd.DataFrame], n: int = 10, reference_year: int = 2025
) -> list[str]:
    """Get top N tissues by count in reference year.

    Args:
        datasets: Dictionary mapping year to DataFrame.
        n: Number of top tissues to return.
        reference_year: Year to use for ranking.

    Returns:
        List of top tissue names.
    """
    if reference_year not in datasets:
        reference_year = max(datasets.keys())

    counts = explode_source_tissues(datasets[reference_year])
    return counts.head(n).index.tolist()


def plot_source_tissue_alluvial(
    datasets: dict[int, pd.DataFrame],
    output_path: Path,
    top_n: int = 10,
    years: list[int] = YEARS,
) -> None:
    """Create alluvial diagram showing source tissue evolution.

    Args:
        datasets: Dictionary mapping year to DataFrame.
        output_path: Path to save the figure.
        top_n: Number of top tissues to show individually.
        years: Years to display as columns.
    """
    # Filter to requested years
    years = [y for y in years if y in datasets]
    if len(years) < 2:
        raise ValueError("Need at least 2 years of data")

    reference_year = years[-1]

    # Get tissue counts for each year
    tissue_counts_by_year = {}
    for year in years:
        tissue_counts_by_year[year] = explode_source_tissues(datasets[year])

    # Get top tissues from reference year
    top_tissues = get_top_tissues(datasets, n=top_n, reference_year=reference_year)

    # Group non-top tissues as "Other" and count how many
    data_by_year = {}
    other_tissue_counts = {}  # Track number of tissues grouped as "Other" per year
    for year in years:
        counts = tissue_counts_by_year[year]
        grouped = {}
        other_count = 0
        other_tissues = 0

        for tissue, count in counts.items():
            if tissue in top_tissues:
                grouped[tissue] = count
            else:
                other_count += count
                other_tissues += 1

        if other_count > 0:
            grouped["Other"] = other_count
        other_tissue_counts[year] = other_tissues

        data_by_year[year] = grouped

    # Determine tissue order: top tissues sorted by reference year count, then Other
    tissue_order = top_tissues + ["Other"]

    # Create label for "Other" showing count from reference year
    n_other = other_tissue_counts.get(reference_year, 0)
    other_label = f"Other ({n_other} tissues)" if n_other > 0 else "Other"

    # Assign colors
    tissue_colors = {}
    for i, tissue in enumerate(top_tissues):
        tissue_colors[tissue] = CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)]
    tissue_colors["Other"] = OTHER_COLOR

    # Figure setup
    n_tissues = len(tissue_order)
    fig_height = max(8, min(14, n_tissues * 0.8))
    fig, ax = plt.subplots(figsize=(12, fig_height))

    # Layout parameters
    n_years = len(years)
    x_positions = [i * 1.8 for i in range(n_years)]
    box_width = 0.4
    gap = 2

    def draw_alluvial_flow(
        x1: float,
        y1_bottom: float,
        y1_top: float,
        x2: float,
        y2_bottom: float,
        y2_top: float,
        color: str,
        alpha: float = 0.4,
    ):
        """Draw a curved flow between two vertical segments."""
        ctrl_offset = (x2 - x1) * 0.4
        verts = [
            (x1 + box_width / 2, y1_bottom),
            (x1 + box_width / 2 + ctrl_offset, y1_bottom),
            (x2 - box_width / 2 - ctrl_offset, y2_bottom),
            (x2 - box_width / 2, y2_bottom),
            (x2 - box_width / 2, y2_top),
            (x2 - box_width / 2 - ctrl_offset, y2_top),
            (x1 + box_width / 2 + ctrl_offset, y1_top),
            (x1 + box_width / 2, y1_top),
            (x1 + box_width / 2, y1_bottom),
        ]
        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.LINETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CLOSEPOLY,
        ]
        path = MplPath(verts, codes)
        patch = PathPatch(path, facecolor=color, edgecolor="none", alpha=alpha)
        ax.add_patch(patch)

    # Use reference year total for scaling
    reference_total = sum(data_by_year[reference_year].values())
    scale_factor = 100 / reference_total

    def needs_top_label(count):
        return count * scale_factor < 2.5

    # Calculate heights for centering
    year_heights = {}
    for year in years:
        total_height = 0
        counts = data_by_year[year]
        for tissue in tissue_order:
            count = counts.get(tissue, 0)
            if count > 0:
                total_height += count * scale_factor
        n_with_data = sum(1 for t in tissue_order if counts.get(t, 0) > 0)
        n_small = sum(
            1 for t in tissue_order if needs_top_label(counts.get(t, 0)) and counts.get(t, 0) > 0
        )
        total_height += (n_with_data - 1) * gap + n_small * 1.5
        year_heights[year] = total_height

    max_height = max(year_heights.values())
    center_y = max_height / 2

    # Track positions for drawing flows
    positions = {}

    # Draw bars for each year
    for i, year in enumerate(years):
        counts = data_by_year[year]
        x = x_positions[i]

        y_offset = center_y - year_heights[year] / 2
        positions[year] = {}

        for tissue in tissue_order:
            count = counts.get(tissue, 0)
            if count == 0:
                continue

            height = count * scale_factor
            is_small = needs_top_label(count)
            color = tissue_colors[tissue]

            positions[year][tissue] = {
                "bottom": y_offset,
                "top": y_offset + height,
                "count": count,
            }

            # Draw rectangle
            rect = Rectangle(
                (x - box_width / 2, y_offset),
                box_width,
                height,
                facecolor=color,
                edgecolor="black",
                linewidth=1,
            )
            ax.add_patch(rect)

            # Add count label
            if is_small:
                # Label on top
                ax.text(
                    x,
                    y_offset + height + 0.5,
                    f"{count:,}",
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    fontweight="bold",
                )
            else:
                # Label inside
                ax.text(
                    x,
                    y_offset + height / 2,
                    f"{count:,}",
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="bold",
                    color="white",
                )

            y_offset += height + gap
            if is_small:
                y_offset += 1.5

    # Draw flows between consecutive years
    for i in range(len(years) - 1):
        year1, year2 = years[i], years[i + 1]
        x1, x2 = x_positions[i], x_positions[i + 1]

        for tissue in tissue_order:
            if tissue in positions[year1] and tissue in positions[year2]:
                pos1 = positions[year1][tissue]
                pos2 = positions[year2][tissue]
                color = tissue_colors[tissue]

                draw_alluvial_flow(
                    x1,
                    pos1["bottom"],
                    pos1["top"],
                    x2,
                    pos2["bottom"],
                    pos2["top"],
                    color,
                    alpha=0.4,
                )

    # Add tissue labels on the right
    last_year = years[-1]
    last_x = x_positions[-1]
    for tissue in tissue_order:
        if tissue in positions[last_year]:
            pos = positions[last_year][tissue]
            color = tissue_colors[tissue]
            # Use descriptive label for "Other"
            label = other_label if tissue == "Other" else tissue
            ax.text(
                last_x + box_width / 2 + 0.15,
                (pos["bottom"] + pos["top"]) / 2,
                label,
                ha="left",
                va="center",
                fontsize=11,
                fontweight="bold",
                color=color,
            )

    # Add year labels at bottom
    for i, year in enumerate(years):
        total = sum(data_by_year[year].values())
        ax.text(
            x_positions[i],
            -5,
            f"{year}\n(n={total:,})",
            ha="center",
            va="top",
            fontsize=12,
            fontweight="bold",
        )

    # Title
    ax.set_title(
        "Source Tissue Evolution in ToxProt",
        fontsize=16,
        fontweight="bold",
        pad=20,
    )

    # Clean up axes
    ax.set_xlim(-0.5, x_positions[-1] + 2)
    y_min = min(positions[y][t]["bottom"] for y in years for t in positions[y])
    y_max = max(positions[y][t]["top"] for y in years for t in positions[y])
    ax.set_ylim(y_min - 10, y_max + 5)
    ax.axis("off")

    # Save
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    """Main analysis pipeline."""
    from .helpers import filter_by_definition

    output_dir = FIGURES_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_datasets(YEARS)
    print(f"Loaded {len(datasets)} datasets")

    # Filter for venom_tissue definition (includes "both")
    datasets = {year: filter_by_definition(df, "venom_tissue") for year, df in datasets.items()}

    plot_source_tissue_alluvial(
        datasets,
        output_dir / "source_tissue_evolution.png",
        top_n=10,
    )
    print(f"Saved to {output_dir}")


if __name__ == "__main__":
    main()
