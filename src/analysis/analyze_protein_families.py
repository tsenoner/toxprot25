#!/usr/bin/env python3
"""Analyze and visualize protein family distributions in ToxProt datasets.

Generates visualizations comparing protein family distributions across
ToxProt datasets from different years (2005, 2015, 2025).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath

from .colors import CATEGORICAL_PALETTE, NAN_COLOR, OTHER_COLOR

# Comprehensive mapping of old family names to 2025 standardized names
# This handles name changes from 2005, 2015, and 2017 datasets
FAMILY_NAME_MAP = {
    # Snake toxins
    "Snake toxin family": "Snake three-finger toxin family",
    "Snake toxin myotoxin family": "Snake three-finger toxin family",
    # Scorpion toxins
    "Alpha/beta-scorpion toxin family": "Long (4 C-C) scorpion toxin superfamily",
    "Short scorpion toxin family": "Short scorpion toxin superfamily",
    "Scorpion toxin 6 family": "Long (3 C-C) scorpion toxin superfamily",
    # Neurotoxin standardization (spider families renamed in 2020s)
    "Huwentoxin-1 family": "Neurotoxin 10 (Hwtx-1) family",
    "Huwentoxin-2 family": "Neurotoxin 12 (Hwtx-2) family",
    "Magi-1 family": "Neurotoxin 14 (magi-1) family",
    "Spider toxin CSTX superfamily": "Neurotoxin 19 (CSTX) family",
    "Spider toxin Tx2 family": "Neurotoxin 03 (Tx2) family",
    "Spider toxin Tx3 family": "Neurotoxin 10 (Hwtx-1) family",  # Merged into Hwtx-1
    "Spider toxin Tx3-6 family": "Neurotoxin 09 (Tx3-6) family",
    "Omega-agatoxin family": "Neurotoxin 04 (omega-agtx) family",
    "Omega-agatoxin superfamily": "Neurotoxin 04 (omega-agtx) family",
    "Beta/delta-agatoxin family": "Neurotoxin 07 (Beta/delta-agtx) family",
    "Mu-agatoxin family": "Neurotoxin 04 (omega-agtx) family",  # Part of omega-agtx
    "U2-agatoxin family": "Neurotoxin 01 (U2-agtx) family",
    "Delta-atracotoxin family": "Neurotoxin 06 (delta-actx) family",
    "Omega-atracotoxin type 1 family": "Neurotoxin 04 (omega-agtx) family",
    "Omega-atracotoxin type 2 family": "Neurotoxin 04 (omega-agtx) family",
    "Plectoxin family": "Neurotoxin 02 (plectoxin) family",
    "Plectoxin superfamily": "Neurotoxin 02 (plectoxin) family",
    "Spider toxin SFI family": "Neurotoxin 16 (SFI) family",
    "Spider neurotoxin 21C2 family": "Neurotoxin 17 (21C2) family",
    "Shiva superfamily": "Neurotoxin 08 (Shiva) family",
    "Spider agouti family": "Neurotoxin 05 (agouti) family",
    "Insecticidal toxin ABC family": "Neurotoxin 13 (insecticidal toxin ABC) family",
    "Omega-lycotoxin family": "Neurotoxin omega-lctx family",
    # Metalloproteinases
    "Peptidase M12B family": "Venom metalloproteinase (M12B) family",
    # Conotoxin superfamilies (standardized naming)
    "Conotoxin O-superfamily": "Conotoxin O1 superfamily",
    "Conotoxin A-superfamily": "Conotoxin A superfamily",
    "Conotoxin M-superfamily": "Conotoxin M superfamily",
    "Conotoxin T-superfamily": "Conotoxin T superfamily",
    "Conotoxin I-superfamily": "Conotoxin I1 superfamily",
    "Conotoxin P-superfamily": "Conotoxin P superfamily",
    # Kunitz-type inhibitors
    "Contains 1 BPTI/Kunitz inhibitor domain": "Venom Kunitz-type family",
    # Spider phospholipases
    "Spider sphingomyelinase family": "Arthropod phospholipase D family",
    # Sea anemone toxins
    "Sea anemone potassium channel inhibitory toxin family": "Venom Kunitz-type family",
    "Sea anemone actinoporin family": "Actinoporin family",
    # Cationic peptides and antimicrobial
    "Latrotoxin family": "Cationic peptide 01 (latrotoxin) family",
    "Latarcin family": "Cationic peptide 06 (latarcin) family",
    "Cupiennin family": "Cationic peptide 08 (cupiennin) family",
    "Cytoinsectotoxin family": "Cationic peptide 10 (cytoinsectotoxin) family",
    "Antimicrobial peptide scorpion family": "Non-disulfide-bridged peptide (NDBP) superfamily",
    "Scorpion NDBP 5 family": "Non-disulfide-bridged peptide (NDBP) superfamily",
    "Scorpion BPP family": "Non-disulfide-bridged peptide (NDBP) superfamily",
    "Short cationic antimicrobial peptide family": "Non-disulfide-bridged peptide (NDBP) superfamily",
    # Other mappings
    "Spider potassium channel inhibitory toxin family": "Neurotoxin 10 (Hwtx-1) family",
    "Ergtoxin family": "Long (3 C-C) scorpion toxin superfamily",
    "Snake waprin family": "Waprin family",
    "Spider wap-1 family": "Waprin family",
}

# Font sizes for consistent styling
TITLE_FONTSIZE = 20
AXIS_LABEL_FONTSIZE = 16
TICK_LABEL_FONTSIZE = 12
LEGEND_FONTSIZE = 14

# Default years for comparison
DEFAULT_YEARS = [2005, 2015, 2025]

# Default ToxProt definition filter
DEFAULT_DEFINITION = "venom_tissue"


def filter_by_definition(
    df: pd.DataFrame,
    definition: str | None = None,
) -> pd.DataFrame:
    """Filter dataset by ToxProt definition.

    Args:
        df: DataFrame containing protein data
        definition: ToxProt definition to filter by ("venom_tissue", "kw_toxin", "both", or None for all)

    Returns:
        Filtered DataFrame
    """
    if definition is None:
        return df

    if "ToxProt definition" not in df.columns:
        return df

    if definition == "venom_tissue":
        # Include entries with "venom_tissue" or "both"
        return df[df["ToxProt definition"].isin(["venom_tissue", "both"])]
    elif definition == "kw_toxin":
        # Include entries with "kw_toxin" or "both"
        return df[df["ToxProt definition"].isin(["kw_toxin", "both"])]
    elif definition == "both":
        # Only entries that match both criteria
        return df[df["ToxProt definition"] == "both"]
    else:
        return df


def normalize_family_name(name: str) -> str:
    """Normalize a family name using the mapping dictionary."""
    if pd.isna(name):
        return name
    return FAMILY_NAME_MAP.get(name, name)


def get_reference_families(
    df: pd.DataFrame, column: str = "Protein families", top_n: int = 15
) -> list[str]:
    """Get top N families from a reference dataset (typically 2025).

    Args:
        df: DataFrame containing protein data (reference year)
        column: Column name containing family information
        top_n: Number of top families to return

    Returns:
        List of top N family names (normalized), sorted by count descending
    """
    family_col = df[column].apply(normalize_family_name)
    counts = family_col.value_counts()
    return counts.nlargest(top_n).index.tolist()


def get_family_counts_for_reference(
    df: pd.DataFrame,
    reference_families: list[str],
    column: str = "Protein families",
) -> pd.Series:
    """Count entries for reference families, grouping everything else as 'Other'.

    Args:
        df: DataFrame containing protein data
        reference_families: List of family names to track (from reference year)
        column: Column name containing family information

    Returns:
        Series with counts for each reference family, 'Other', and 'Unassigned'
    """
    # Normalize all family names
    family_col = df[column].apply(normalize_family_name)

    # Count entries for each reference family
    result = {}
    for family in reference_families:
        result[family] = (family_col == family).sum()

    # Count "Other" (families not in reference list, excluding NaN)
    non_nan_mask = family_col.notna()
    in_reference_mask = family_col.isin(reference_families)
    other_count = (~in_reference_mask & non_nan_mask).sum()

    # Count "Unassigned" (NaN entries)
    unassigned_count = family_col.isna().sum()

    if other_count > 0:
        result["Other"] = other_count
    if unassigned_count > 0:
        result["Unassigned"] = unassigned_count

    return pd.Series(result)


def _create_color_map(reference_families: list[str]) -> dict:
    """Create color mapping for protein families using centralized colors.

    Colors are assigned based on 2025 reference families. Old family names
    that map to a reference family get the same color.

    Args:
        reference_families: List of canonical family names from reference year

    Returns:
        Dictionary mapping family names to colors
    """
    color_map = {
        "Other": OTHER_COLOR,
        "Unassigned": NAN_COLOR,
        "NaN": NAN_COLOR,  # Keep for backward compatibility
    }

    # Assign colors to reference families
    for i, family in enumerate(reference_families):
        color_map[family] = CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)]

    # Also assign colors to old names that map to reference families
    for old_name, new_name in FAMILY_NAME_MAP.items():
        if new_name in color_map:
            color_map[old_name] = color_map[new_name]

    return color_map


def plot_stacked_bar_protein_families(
    datasets: dict[int, pd.DataFrame],
    output_path: Path,
    top_n: int = 15,
):
    """
    Create 100% stacked bar chart of top protein families.

    Uses the most recent year as the reference for which families to track.
    Families are sorted with the largest (in reference year) at the bottom.

    Args:
        datasets: Dictionary mapping year (int) to DataFrame
        output_path: Path to save the figure
        top_n: Number of top families to display
    """
    years = sorted(datasets.keys())
    reference_year = years[-1]  # Most recent year (e.g., 2025)

    # Get reference families from the most recent year
    reference_families = get_reference_families(
        datasets[reference_year], "Protein families", top_n=top_n
    )

    # Get counts for each year using reference families
    family_counts = {}
    for year in years:
        family_counts[year] = get_family_counts_for_reference(
            datasets[year], reference_families, "Protein families"
        )

    # Create ordered family list: ALL categories sorted by 2025 count (largest at bottom)
    # Get 2025 counts for sorting
    counts_2025 = family_counts[reference_year]

    # Build list of all families including Other and Unassigned
    all_categories = list(reference_families)
    if any("Other" in family_counts[y].index for y in years):
        all_categories.append("Other")
    if any("Unassigned" in family_counts[y].index for y in years):
        all_categories.append("Unassigned")

    # Sort by 2025 count descending (largest first = bottom of stack)
    ordered_families = sorted(all_categories, key=lambda f: counts_2025.get(f, 0), reverse=True)

    # Build data frame with counts and percentages
    df_counts = pd.DataFrame(index=ordered_families)
    for year in years:
        df_counts[str(year)] = family_counts[year].reindex(ordered_families).fillna(0)

    df_percentages = df_counts.div(df_counts.sum(axis=0), axis=1) * 100

    # Create color mapping based on reference families
    color_map = _create_color_map(reference_families)

    # Create stacked bar chart
    fig, ax = plt.subplots(figsize=(14, 9))

    def create_stacked_bars(year_label: str, year_col: str) -> list[dict]:
        """Create stacked bars for a single year column."""
        bars_info = []
        bottom = 0
        for family in ordered_families:
            pct = df_percentages.at[family, year_col]
            cnt = df_counts.at[family, year_col]
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

    # Create bars for all years
    all_bars = []
    for year in years:
        year_label = str(year)
        bars = create_stacked_bars(year_label, str(year))
        all_bars.extend(bars)

    # Add labels to all bars (use smaller font for narrow bars)
    for bar_info in all_bars:
        bar = bar_info["bar"]
        pct = bar_info["percentage"]
        if pct > 0:
            # Use smaller font for narrow bars
            fontsize = 8 if pct > 3 else 6 if pct > 1.5 else 5
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_y() + bar.get_height() / 2,
                f"{bar_info['family']}, {int(bar_info['count'])} ({pct:.1f}%)",
                ha="center",
                va="center",
                fontsize=fontsize,
            )

    # Styling
    ax.set_ylabel("Percentage", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylim(0, 100)
    ax.tick_params(axis="x", labelsize=TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="y", labelsize=TICK_LABEL_FONTSIZE)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


def _get_original_family_counts(
    df: pd.DataFrame,
    reference_families: list[str],
    column: str = "Protein families",
) -> dict[str, int]:
    """Get counts using ORIGINAL (un-normalized) family names.

    Groups families by their 2025 destination and returns original names with counts.

    Args:
        df: DataFrame containing protein data
        reference_families: List of canonical family names from reference year
        column: Column name containing family information

    Returns:
        Dictionary mapping original family names to counts
    """
    family_col = df[column]
    result = {}

    # Create reverse mapping: 2025 name -> list of old names
    reverse_map = {fam: [fam] for fam in reference_families}
    for old_name, new_name in FAMILY_NAME_MAP.items():
        if new_name in reverse_map:
            reverse_map[new_name].append(old_name)

    # Count each original family name
    for family_2025 in reference_families:
        for original_name in reverse_map[family_2025]:
            count = (family_col == original_name).sum()
            if count > 0:
                result[original_name] = count

    # Count "Other" (families not mapping to reference list, excluding NaN)
    all_original_names = set()
    for names in reverse_map.values():
        all_original_names.update(names)

    non_nan_mask = family_col.notna()
    in_reference_mask = family_col.isin(all_original_names)
    other_count = (~in_reference_mask & non_nan_mask).sum()

    # Count "Unassigned" (NaN entries)
    unassigned_count = family_col.isna().sum()

    if other_count > 0:
        result["Other"] = other_count
    if unassigned_count > 0:
        result["Unassigned"] = unassigned_count

    return result


def _get_family_destination(family_name: str, reference_families: list[str]) -> str:
    """Get the 2025 destination family for an original family name.

    Args:
        family_name: Original family name
        reference_families: List of canonical family names from reference year

    Returns:
        The 2025 family name this maps to, or the original if no mapping exists
    """
    if family_name in ("Other", "Unassigned"):
        return family_name
    normalized = normalize_family_name(family_name)
    if normalized in reference_families:
        return normalized
    return "Other"


def plot_alluvial_protein_families(
    datasets: dict[int, pd.DataFrame],
    output_path: Path,
    top_n: int = 15,
):
    """
    Create alluvial plot showing protein family name evolution across years.

    Shows original family names from each year, colored by their 2025 destination.
    Flows connect families based on name mapping.

    Args:
        datasets: Dictionary mapping year (int) to DataFrame
        output_path: Path to save the figure
        top_n: Number of top families to display
    """
    years = sorted(datasets.keys())
    reference_year = years[-1]  # Most recent year (e.g., 2025)

    # Get reference families from the most recent year
    reference_families = get_reference_families(
        datasets[reference_year], "Protein families", top_n=top_n
    )

    # Get original family counts for each year
    original_counts = {}
    for year in years:
        original_counts[year] = _get_original_family_counts(
            datasets[year], reference_families, "Protein families"
        )

    # Create color mapping based on reference families
    color_map = _create_color_map(reference_families)

    # Get 2025 counts for destination-based sorting
    counts_2025 = original_counts[reference_year]

    # Calculate total count per destination in 2025
    dest_totals_2025 = {}
    for fam, count in counts_2025.items():
        dest = _get_family_destination(fam, reference_families)
        dest_totals_2025[dest] = dest_totals_2025.get(dest, 0) + count

    # Sort families by their 2025 destination's total count (largest first = bottom)
    def sort_key(family_name: str, counts: dict) -> tuple:
        """Sort families by 2025 destination count (largest first), then by count within."""
        dest = _get_family_destination(family_name, reference_families)
        dest_count = dest_totals_2025.get(dest, 0)
        return (-dest_count, -counts.get(family_name, 0))

    # Layout parameters - wide and short aspect ratio
    fig, ax = plt.subplots(figsize=(24, 8))
    n_years = len(years)
    x_positions = [i * 1.8 for i in range(n_years)]
    box_width = 0.6

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

    # Calculate percentages and positions for each year
    positions = {}
    sorted_families_by_year = {}

    for i, year in enumerate(years):
        counts = original_counts[year]
        total = sum(counts.values())

        # Sort families for this year
        families_sorted = sorted(counts.keys(), key=lambda f: sort_key(f, counts))
        sorted_families_by_year[year] = families_sorted

        # Calculate percentages
        percentages = {f: (c / total * 100) for f, c in counts.items()}

        x = x_positions[i]
        y_offset = 0
        positions[year] = {}

        for family in families_sorted:
            pct = percentages.get(family, 0)
            count = counts.get(family, 0)

            positions[year][family] = {
                "bottom": y_offset,
                "top": y_offset + pct,
                "percentage": pct,
                "count": count,
            }

            if pct > 0:
                # Get color based on 2025 destination
                dest = _get_family_destination(family, reference_families)
                bar_color = color_map.get(dest, OTHER_COLOR)

                # Draw rectangle
                rect = Rectangle(
                    (x - box_width / 2, y_offset),
                    box_width,
                    pct,
                    facecolor=bar_color,
                    edgecolor="white",
                    linewidth=0.5,
                )
                ax.add_patch(rect)

                # Add labels: outside for first/last year, inside for middle years
                display_name = family.replace(" family", "").replace(" superfamily", "")

                is_first_year = i == 0
                is_last_year = i == n_years - 1

                # Get color for this family's destination
                dest = _get_family_destination(family, reference_families)
                label_color = color_map.get(dest, OTHER_COLOR)

                # Inside label: count and percentage for all years
                inside_fontsize = 10 if pct > 5 else 8 if pct > 2.5 else 6
                ax.text(
                    x,
                    y_offset + pct / 2,
                    f"{int(count)} ({pct:.1f}%)",
                    ha="center",
                    va="center",
                    fontsize=inside_fontsize,
                    fontweight="bold",
                    color="black",
                )

                # Outside labels for first and last year
                if is_first_year:
                    # Store label info for later (to handle overlaps)
                    if "left_labels" not in positions[year]:
                        positions[year]["left_labels"] = []
                    positions[year]["left_labels"].append(
                        {
                            "name": display_name,
                            "y": y_offset + pct / 2,
                            "color": label_color,
                            "bar_y": y_offset + pct / 2,
                        }
                    )
                elif is_last_year:
                    # Labels to the RIGHT of the bar - full name, colored text
                    ax.text(
                        x + box_width / 2 + 0.1,
                        y_offset + pct / 2,
                        display_name,
                        ha="left",
                        va="center",
                        fontsize=12,
                        fontweight="bold",
                        color=label_color,
                    )

            y_offset += pct

    # Draw left labels for first year with connecting lines to prevent overlap
    first_year = years[0]
    if "left_labels" in positions[first_year]:
        labels = positions[first_year]["left_labels"]
        x = x_positions[0]

        # Calculate minimum spacing between labels (in percentage units)
        min_spacing = 3.5  # Minimum vertical gap between label centers

        # Adjust label positions to prevent overlap - push DOWN instead of up
        adjusted_labels = []
        for label in labels:
            adjusted_labels.append(
                {
                    **label,
                    "adjusted_y": label["y"],
                }
            )

        # Sort by original y position (highest first for downward adjustment)
        adjusted_labels.sort(key=lambda lbl: lbl["y"], reverse=True)

        # Push labels down if they overlap (iterate from top to bottom)
        for i in range(1, len(adjusted_labels)):
            prev_y = adjusted_labels[i - 1]["adjusted_y"]
            curr_y = adjusted_labels[i]["adjusted_y"]
            if prev_y - curr_y < min_spacing:
                adjusted_labels[i]["adjusted_y"] = prev_y - min_spacing

        # Draw labels with connecting lines
        for label in adjusted_labels:
            # Draw connecting line from bar to label
            ax.plot(
                [x - box_width / 2 - 0.05, x - box_width / 2 - 0.4],
                [label["bar_y"], label["adjusted_y"]],
                color=label["color"],
                linewidth=1.5,
                alpha=0.6,
            )

            # Draw label text
            ax.text(
                x - box_width / 2 - 0.5,
                label["adjusted_y"],
                label["name"],
                ha="right",
                va="center",
                fontsize=12,
                fontweight="bold",
                color=label["color"],
            )

    # Draw flows between years
    for i in range(n_years - 1):
        year1 = years[i]
        year2 = years[i + 1]
        x1 = x_positions[i]
        x2 = x_positions[i + 1]

        # Group source families by their 2025 destination
        flows_by_dest = {}
        for family1 in sorted_families_by_year[year1]:
            dest = _get_family_destination(family1, reference_families)
            pos1 = positions[year1][family1]
            if pos1["percentage"] > 0:
                if dest not in flows_by_dest:
                    flows_by_dest[dest] = []
                flows_by_dest[dest].append(
                    {
                        "family": family1,
                        "bottom": pos1["bottom"],
                        "top": pos1["top"],
                        "count": pos1["count"],
                    }
                )

        # For each destination, calculate proportional landing positions
        for dest, sources in flows_by_dest.items():
            # Find destination families in year2
            dest_families = [
                f
                for f in sorted_families_by_year[year2]
                if _get_family_destination(f, reference_families) == dest
            ]

            if not dest_families:
                continue

            # Get destination bar bounds and total count
            y2_bottom = min(positions[year2][f]["bottom"] for f in dest_families)
            y2_top = max(positions[year2][f]["top"] for f in dest_families)
            dest_height = y2_top - y2_bottom
            dest_total_count = sum(positions[year2][f]["count"] for f in dest_families)

            # Draw each source flow - landing height proportional to destination total
            current_y2 = y2_bottom
            for source in sources:
                # Landing height = (source_count / dest_total_count) * dest_height
                # This shows what fraction of the destination this source represents
                landing_height = (
                    (source["count"] / dest_total_count) * dest_height
                    if dest_total_count > 0
                    else 0
                )

                bar_color = color_map.get(dest, OTHER_COLOR)
                draw_alluvial_flow(
                    x1,
                    source["bottom"],
                    source["top"],
                    x2,
                    current_y2,
                    current_y2 + landing_height,
                    bar_color,
                    alpha=0.3,
                )
                current_y2 += landing_height

    # Add year labels with totals
    for i, year in enumerate(years):
        total = sum(original_counts[year].values())
        ax.text(
            x_positions[i],
            -3,
            f"{year}\n(n={total:,})",
            ha="center",
            va="top",
            fontsize=14,
            fontweight="bold",
        )

    # Formatting - tight xlim to minimize whitespace
    ax.set_xlim(-2.6, x_positions[-1] + 2)
    ax.set_ylim(-8, 105)
    ax.axis("off")

    plt.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close()


def main():
    """Main function for standalone execution."""
    parser = argparse.ArgumentParser(
        description="Generate protein family visualizations comparing multiple years."
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
    parser.add_argument(
        "--years",
        type=int,
        nargs="+",
        default=DEFAULT_YEARS,
        help=f"Years to compare (default: {DEFAULT_YEARS})",
    )
    parser.add_argument(
        "--definition",
        type=str,
        default=DEFAULT_DEFINITION,
        choices=["venom_tissue", "kw_toxin", "both", "all"],
        help=f"ToxProt definition filter (default: {DEFAULT_DEFINITION})",
    )

    args = parser.parse_args()

    # Load datasets for all requested years
    definition_filter = None if args.definition == "all" else args.definition
    print(f"Loading ToxProt datasets (definition: {args.definition})...")
    datasets = {}
    for year in args.years:
        filepath = args.data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            df = pd.read_csv(filepath)
            df_filtered = filter_by_definition(df, definition_filter)
            datasets[year] = df_filtered
            print(f"  {year}: {len(df_filtered):,} entries (of {len(df):,} total)")
        else:
            print(f"Warning: {filepath} not found, skipping year {year}")

    if len(datasets) < 2:
        print("Error: Need at least 2 datasets for comparison")
        return

    # Ensure output directory exists
    protein_families_dir = args.output_dir / "protein_families"
    protein_families_dir.mkdir(parents=True, exist_ok=True)

    # Generate stacked bar chart
    stacked_bar_path = protein_families_dir / "top_families_stacked_bar.png"
    plot_stacked_bar_protein_families(datasets, stacked_bar_path, top_n=args.top_n)
    print(f"Saved {stacked_bar_path}")

    # Generate alluvial plot
    alluvial_path = protein_families_dir / "top_families_alluvial.png"
    plot_alluvial_protein_families(datasets, alluvial_path, top_n=args.top_n)
    print(f"Saved {alluvial_path}")


if __name__ == "__main__":
    main()
