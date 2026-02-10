#!/usr/bin/env python3
"""
Protein Evidence Alluvial Plot
==============================

Creates an alluvial-style diagram showing the flow of protein evidence categories
across ToxProt datasets from 2008 to 2015 to 2025.

Uses matplotlib with Bezier curve flows, consistent with other project figures.

Author: Tobias Senoner
"""

from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath

from .helpers import filter_by_definition, load_datasets

# --- Configuration ---
FIGURES_DIR = Path("figures/protein_evidence")
YEARS = [2008, 2015, 2025]

# Protein existence category colors
PE_COLORS = {
    "Evidence at protein level": "#2E86AB",  # Blue
    "Evidence at transcript level": "#A23B72",  # Purple
    "Inferred from homology": "#F18F01",  # Orange
    "Predicted": "#C73E1D",  # Red
    "Uncertain": "#592941",  # Dark purple
    "Removed": "#808080",  # Gray for removed entries
}

# Canonical order for PE categories (most certain to least certain)
PE_CATEGORIES = [
    "Evidence at protein level",
    "Evidence at transcript level",
    "Inferred from homology",
    "Predicted",
    "Uncertain",
]

# Display labels (two lines for long names)
PE_LABELS = {
    "Evidence at protein level": "Evidence at\nprotein level",
    "Evidence at transcript level": "Evidence at\ntranscript level",
    "Inferred from homology": "Inferred from\nhomology",
    "Predicted": "Predicted",
    "Uncertain": "Uncertain",
}


def normalize_pe_category(category: str) -> str:
    """
    Normalize protein existence categories to standard names.

    Args:
        category: Raw category from dataset

    Returns:
        Normalized category name
    """
    if pd.isna(category):
        return "Unknown"

    category = str(category).strip()

    # Map numbered categories to readable names
    if category.startswith("1:") or "protein level" in category.lower():
        return "Evidence at protein level"
    elif category.startswith("2:") or "transcript level" in category.lower():
        return "Evidence at transcript level"
    elif category.startswith("3:") or "homology" in category.lower():
        return "Inferred from homology"
    elif category.startswith("4:") or "predicted" in category.lower():
        return "Predicted"
    elif category.startswith("5:") or "uncertain" in category.lower():
        return "Uncertain"
    else:
        return category


def plot_protein_evidence_alluvial(
    datasets: dict[int, pd.DataFrame],
    output_path: Path,
    years: list[int] = YEARS,
    definition: str | None = "venom_tissue",
):
    """
    Create alluvial-style diagram showing protein evidence category transitions.

    Style matches plot_newcomers_alluvial() from analyze_taxa.py:
    - Centered nodes with black edges
    - White counts inside large bars, colored counts above small bars
    - Category names on right side of last column
    - Bezier curve flows between years
    - Year labels with totals at bottom
    - "Removed" entries shown midway between columns

    Args:
        datasets: Dictionary mapping year to DataFrame
        output_path: Path to save the figure
        years: Years to show as columns
        definition: ToxProt definition filter
    """
    # Apply definition filter
    filtered_datasets = {}
    for year in years:
        if year not in datasets:
            print(f"Warning: Dataset for {year} not found")
            continue
        filtered_datasets[year] = filter_by_definition(datasets[year], definition)

    # Prepare data: normalize PE categories and track protein IDs
    pe_by_year = {}
    for year in years:
        if year not in filtered_datasets:
            continue
        df = filtered_datasets[year].copy()

        # Determine protein ID column
        id_col = "Entry" if "Entry" in df.columns else df.columns[0]

        # Normalize PE categories
        df["PE_normalized"] = df["Protein existence"].apply(normalize_pe_category)
        pe_by_year[year] = df.set_index(id_col)["PE_normalized"]

    # Get category counts for each year
    data_by_year = {}
    for year in years:
        if year not in pe_by_year:
            continue
        counts = pe_by_year[year].value_counts().to_dict()
        # Filter to only canonical categories (ignore Unknown/other)
        data_by_year[year] = {
            cat: counts.get(cat, 0) for cat in PE_CATEGORIES if counts.get(cat, 0) > 0
        }

    # Calculate flows between consecutive years
    flows = {}
    for i in range(len(years) - 1):
        year1, year2 = years[i], years[i + 1]
        if year1 not in pe_by_year or year2 not in pe_by_year:
            continue

        pe1, pe2 = pe_by_year[year1], pe_by_year[year2]
        flow_counter = Counter()

        # Common proteins - track category transitions
        common_ids = set(pe1.index) & set(pe2.index)
        for pid in common_ids:
            src, tgt = pe1[pid], pe2[pid]
            if src in PE_CATEGORIES and tgt in PE_CATEGORIES:
                flow_counter[(src, tgt)] += 1

        # Removed proteins (in year1 but not in year2)
        removed_ids = set(pe1.index) - set(pe2.index)
        removed_counts = Counter()
        for pid in removed_ids:
            cat = pe1[pid]
            if cat in PE_CATEGORIES:
                removed_counts[cat] += 1

        flows[(year1, year2)] = {
            "transitions": flow_counter,
            "removed": removed_counts,
            "removed_total": sum(removed_counts.values()),
        }

    # Figure size (compact height)
    fig, ax = plt.subplots(figsize=(12, 7))

    # Layout parameters - add space for "Removed" nodes between years
    n_years = len(years)
    year_x_positions = [i * 2.4 for i in range(n_years)]
    removed_x_positions = [year_x_positions[i] + 1.2 for i in range(n_years - 1)]
    box_width = 0.7  # Wider to fit count + percentage
    removed_box_width = 0.3  # Slightly narrower for "Removed" nodes
    gap = 3  # Gap between stacked bars

    def draw_alluvial_flow(
        x1: float,
        y1_bottom: float,
        y1_top: float,
        x2: float,
        y2_bottom: float,
        y2_top: float,
        color: str,
        alpha: float = 0.4,
        width1: float = box_width,
        width2: float = box_width,
    ):
        """Draw a curved flow between two vertical segments."""
        ctrl_offset = (x2 - x1) * 0.4
        verts = [
            (x1 + width1 / 2, y1_bottom),
            (x1 + width1 / 2 + ctrl_offset, y1_bottom),
            (x2 - width2 / 2 - ctrl_offset, y2_bottom),
            (x2 - width2 / 2, y2_bottom),
            (x2 - width2 / 2, y2_top),
            (x2 - width2 / 2 - ctrl_offset, y2_top),
            (x1 + width1 / 2 + ctrl_offset, y1_top),
            (x1 + width1 / 2, y1_top),
            (x1 + width1 / 2, y1_bottom),
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

    # Use 2025 total as reference for scaling
    reference_year = years[-1]
    reference_total = sum(data_by_year[reference_year].values())
    scale_factor = 100 / reference_total

    # Determine which categories need labels on top (small bars)
    def needs_top_label(count):
        return count * scale_factor < 3.0

    # Calculate total heights for each year (for centering)
    year_heights = {}
    for year in years:
        if year not in data_by_year:
            continue
        total_height = 0
        counts = data_by_year[year]
        n_cats = len([c for c in counts.values() if c > 0])
        for count in counts.values():
            if count > 0:
                total_height += count * scale_factor
        # Add gaps
        n_small = sum(1 for c in counts.values() if c > 0 and needs_top_label(c))
        total_height += (n_cats - 1) * gap + n_small * 2.0
        year_heights[year] = total_height

    max_height = max(year_heights.values())
    center_y = max_height / 2

    # Track positions for each year and category
    positions = {}
    removed_positions = {}  # For "Removed" nodes

    # Draw bars for each year (centered)
    for i, year in enumerate(years):
        if year not in data_by_year:
            continue

        counts = data_by_year[year]
        x = year_x_positions[i]

        # Sort categories in canonical order
        sorted_cats = [(cat, counts.get(cat, 0)) for cat in PE_CATEGORIES if counts.get(cat, 0) > 0]

        # Calculate starting y offset for centering
        y_offset = center_y - year_heights[year] / 2
        positions[year] = {}

        for cat, count in sorted_cats:
            height = count * scale_factor
            is_small = needs_top_label(count)

            positions[year][cat] = {
                "bottom": y_offset,
                "top": y_offset + height,
                "count": count,
                "height": height,
            }

            color = PE_COLORS.get(cat, "#CCCCCC")

            # Draw rectangle with black edge
            rect = Rectangle(
                (x - box_width / 2, y_offset),
                box_width,
                height,
                facecolor=color,
                edgecolor="black",
                linewidth=1,
            )
            ax.add_patch(rect)

            # Calculate percentage for this year
            year_total = sum(counts.values())
            pct = (count / year_total * 100) if year_total > 0 else 0

            # Show count and percentage - inside for large bars, on top for small bars
            if not is_small:
                fontsize = 11 if height > 5 else 9
                label = f"{count:,}\n({pct:.1f}%)"
                ax.text(
                    x,
                    y_offset + height / 2,
                    label,
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    fontweight="bold",
                    color="white",
                )
            else:
                # Small bar - put count on top
                label = f"{count:,} ({pct:.1f}%)"
                ax.text(
                    x,
                    y_offset + height + 0.5,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    fontweight="bold",
                    color=color,
                )

            # Extra gap after small bars
            extra_gap = 2.0 if is_small else 0
            y_offset += height + gap + extra_gap

    # Draw "Removed" nodes between years
    for i in range(len(years) - 1):
        year1, year2 = years[i], years[i + 1]
        if (year1, year2) not in flows:
            continue

        flow_data = flows[(year1, year2)]
        removed_total = flow_data["removed_total"]

        if removed_total == 0:
            continue

        removed_x = removed_x_positions[i]
        removed_height = removed_total * scale_factor

        # Position "Removed" node - center it vertically
        removed_bottom = center_y - removed_height / 2
        removed_top = removed_bottom + removed_height

        removed_positions[(year1, year2)] = {
            "bottom": removed_bottom,
            "top": removed_top,
            "count": removed_total,
            "height": removed_height,
        }

        # Draw "Removed" rectangle
        rect = Rectangle(
            (removed_x - removed_box_width / 2, removed_bottom),
            removed_box_width,
            removed_height,
            facecolor=PE_COLORS["Removed"],
            edgecolor="black",
            linewidth=1,
        )
        ax.add_patch(rect)

        # Add count inside or above
        if needs_top_label(removed_total):
            ax.text(
                removed_x,
                removed_top + 0.5,
                f"{removed_total:,}",
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold",
                color=PE_COLORS["Removed"],
            )
        else:
            ax.text(
                removed_x,
                (removed_bottom + removed_top) / 2,
                f"{removed_total:,}",
                ha="center",
                va="center",
                fontsize=11,
                fontweight="bold",
                color="white",
            )

        # Add "Removed" label below the bar
        ax.text(
            removed_x,
            removed_bottom - 1.5,
            "Removed",
            ha="center",
            va="top",
            fontsize=9,
            fontweight="bold",
            fontstyle="italic",
            color=PE_COLORS["Removed"],
        )

    # Draw flows between consecutive years
    for i in range(len(years) - 1):
        year1, year2 = years[i], years[i + 1]
        if year1 not in positions or year2 not in positions:
            continue
        if (year1, year2) not in flows:
            continue

        x1, x2 = year_x_positions[i], year_x_positions[i + 1]
        removed_x = removed_x_positions[i]
        flow_data = flows[(year1, year2)]

        # Track cumulative positions for stacking flows within each bar
        src_offsets = {cat: positions[year1][cat]["bottom"] for cat in positions[year1]}
        tgt_offsets = {cat: positions[year2][cat]["bottom"] for cat in positions[year2]}

        # Draw flows for category transitions (existing proteins)
        for src_cat in PE_CATEGORIES:
            if src_cat not in positions[year1]:
                continue

            for tgt_cat in PE_CATEGORIES:
                if tgt_cat not in positions[year2]:
                    continue

                count = flow_data["transitions"].get((src_cat, tgt_cat), 0)
                if count == 0:
                    continue

                flow_height = count * scale_factor

                # Source position
                y1_bottom = src_offsets[src_cat]
                y1_top = y1_bottom + flow_height
                src_offsets[src_cat] = y1_top

                # Target position
                y2_bottom = tgt_offsets[tgt_cat]
                y2_top = y2_bottom + flow_height
                tgt_offsets[tgt_cat] = y2_top

                # Use source category color for the flow
                color = PE_COLORS.get(src_cat, "#CCCCCC")

                draw_alluvial_flow(
                    x1,
                    y1_bottom,
                    y1_top,
                    x2,
                    y2_bottom,
                    y2_top,
                    color,
                    alpha=0.4,
                )

        # Draw flows from PE categories to "Removed" node
        if (year1, year2) in removed_positions:
            removed_pos = removed_positions[(year1, year2)]
            removed_offset = removed_pos["bottom"]

            for src_cat in PE_CATEGORIES:
                if src_cat not in positions[year1]:
                    continue

                count = flow_data["removed"].get(src_cat, 0)
                if count == 0:
                    continue

                flow_height = count * scale_factor

                # Source position
                y1_bottom = src_offsets[src_cat]
                y1_top = y1_bottom + flow_height
                src_offsets[src_cat] = y1_top

                # Removed node position
                y2_bottom = removed_offset
                y2_top = y2_bottom + flow_height
                removed_offset = y2_top

                # Use source category color for the flow
                color = PE_COLORS.get(src_cat, "#CCCCCC")

                draw_alluvial_flow(
                    x1,
                    y1_bottom,
                    y1_top,
                    removed_x,
                    y2_bottom,
                    y2_top,
                    color,
                    alpha=0.4,
                    width1=box_width,
                    width2=removed_box_width,
                )

    # Add category names on the right side of the last column (two-line labels)
    last_year = years[-1]
    last_x = year_x_positions[-1]

    for cat in positions[last_year]:
        pos = positions[last_year][cat]
        color = PE_COLORS.get(cat, "#CCCCCC")
        label = PE_LABELS.get(cat, cat)

        ax.text(
            last_x + box_width / 2 + 0.1,
            (pos["bottom"] + pos["top"]) / 2,
            label,
            ha="left",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=color,
        )

    # Calculate totals for year labels
    totals = {year: sum(data_by_year[year].values()) for year in years if year in data_by_year}

    # Add year labels at bottom
    for i, year in enumerate(years):
        if year in totals:
            total = totals[year]
            ax.text(
                year_x_positions[i],
                -5,
                f"{year}\n(n={total:,})",
                ha="center",
                va="top",
                fontsize=14,
                fontweight="bold",
            )

    ax.set_xlim(-0.8, year_x_positions[-1] + 2.2)
    ax.set_ylim(-12, max_height + 5)
    ax.axis("off")

    # Save figure
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.1)
    plt.close()


def main():
    """Main function for standalone execution."""
    datasets = load_datasets(YEARS)

    if len(datasets) < 2:
        print("Error: Need at least 2 datasets to generate plot")
        return

    print(f"Loaded {len(datasets)} datasets")
    for year, df in datasets.items():
        print(f"  {year}: {len(df)} proteins")

    output_path = FIGURES_DIR / "protein_evidence_sankey.png"
    plot_protein_evidence_alluvial(datasets, output_path, years=YEARS)

    print("Done!")


if __name__ == "__main__":
    main()
