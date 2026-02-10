#!/usr/bin/env python3
"""Analyze and visualize taxonomic distribution in ToxProt datasets."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from PIL import Image

from ..config import ALL_YEARS, COMPARISON_YEARS, FIGURES_DIR
from .colors import (
    OTHER_COLOR,
    TAXA_ORDER_COLORS,
    TAXA_ORDER_COLORS_LIST,
)
from .helpers import filter_by_definition, load_datasets

# --- Configuration ---
SILHOUETTE_DIR = Path("data/raw/PhyloPic/png")
YEARS = ALL_YEARS

# Silhouette mapping: taxa order -> image base name
SILHOUETTE_MAP = {
    "Squamata": "Cobra",
    "Araneae": "Spider",
    "Neogastropoda": "Conus",
    "Scorpiones": "Scorpion",
    "Hymenoptera": "Hymenoptera",
}

# Silhouettes for "Others" category
OTHERS_SILHOUETTES = ["Aedes", "Annelida", "Scolopendra", "Moth"]

# Common names for newcomer taxa orders
ORDER_COMMON_NAMES = {
    "Lepidoptera": "Butterflies and Moths",
    "Scutigeromorpha": "House Centipedes",
    "Rhynchobdellida": "Proboscis Leeches",
    "Zoantharia": "Colonial Anemones",
    "Chiroptera": "Bats",
    "Hirudinida": "Leeches",
    "Xiphosura": "Horseshoe Crabs",
    "Suberitida": "Sponges",
    "Semaeostomeae": "Jellyfish",
    "Rhabditida": "Roundworms",
    "Nectiopoda": "Remipedes",
    "Euphausiacea": "Krill",
}

# Default ToxProt definition filter
DEFAULT_DEFINITION = "venom_tissue"

# Taxa levels for visualization
TAXA_LEVELS = ["Phylum", "Class", "Order"]


def load_silhouette(name: str, color: str = "grey") -> Image.Image | None:
    """Load a silhouette PNG image by name."""
    # Prefer color-suffixed version (e.g., Spider_grey.png) over base file
    path = SILHOUETTE_DIR / f"{name}_{color}.png"
    if path.exists():
        return Image.open(path)
    return None


def get_silhouette(order: str) -> Image.Image | None:
    """Get silhouette image for a taxa order."""
    if order not in SILHOUETTE_MAP:
        return None
    return load_silhouette(SILHOUETTE_MAP[order])


def get_taxa_newcomers(
    df_old: pd.DataFrame, df_new: pd.DataFrame, taxa_level: str = "Order"
) -> pd.Series:
    """Identify taxa present in new dataset but not in old dataset."""
    old_taxa = set(df_old[taxa_level].unique())
    new_taxa = set(df_new[taxa_level].unique())
    newcomers = new_taxa - old_taxa
    return df_new[df_new[taxa_level].isin(newcomers)][taxa_level].value_counts()


def plot_top_taxa_trend(
    datasets: dict[int, pd.DataFrame],
    output_path: Path,
    reference_year: int = 2025,
):
    """Create line plot showing taxa trends over time with silhouettes."""
    ref_df = datasets[reference_year]
    top_orders = ref_df["Order"].value_counts().nlargest(5).index.tolist()

    taxa_data = {}
    for year, df in datasets.items():
        counts = df["Order"].value_counts()
        taxa_data[year] = {order: counts.get(order, 0) for order in top_orders}
        taxa_data[year]["Others"] = counts[~counts.index.isin(top_orders)].sum()

    plot_df = pd.DataFrame(taxa_data).T
    plot_df.index = plot_df.index.astype(int)

    fig, ax = plt.subplots(figsize=(10, 7))
    fig.subplots_adjust(right=0.72)

    # =================================================================
    # SILHOUETTE CONFIGURATION
    # x_offset/y_offset: position in points, zoom: scale, rotation: degrees
    # =================================================================
    silhouette_config = {
        "Squamata": {"x_offset": 90, "y_offset": -10, "zoom": 0.02, "rotation": 0},
        "Araneae": {"x_offset": 45, "y_offset": 35, "zoom": 0.13, "rotation": -90},
        "Neogastropoda": {"x_offset": 45, "y_offset": -22, "zoom": 0.018, "rotation": 90},
        "Scorpiones": {"x_offset": 110, "y_offset": 0, "zoom": 0.09, "rotation": 0},
        "Hymenoptera": {"x_offset": 45, "y_offset": 35, "zoom": 0.18, "rotation": -15},
        # Others silhouettes
        "Scolopendra": {"x_offset": 15, "y_offset": -35, "zoom": 0.018, "rotation": 0},
        "Annelida": {"x_offset": 38, "y_offset": -35, "zoom": 0.009, "rotation": 0},
        "Aedes": {"x_offset": 80, "y_offset": -22, "zoom": 0.035, "rotation": 0},
        "Moth": {"x_offset": 80, "y_offset": -50, "zoom": 0.013, "rotation": 2},
    }

    for i, order in enumerate(plot_df.columns):
        color = TAXA_ORDER_COLORS.get(
            order, TAXA_ORDER_COLORS_LIST[i % len(TAXA_ORDER_COLORS_LIST)]
        )
        ax.plot(plot_df.index, plot_df[order], marker="o", markersize=8, linewidth=2.5, color=color)

        last_value = plot_df[order].iloc[-1]
        config = silhouette_config.get(order, {"x_offset": 45, "y_offset": -35, "zoom": 0.045})

        # Add text label
        ax.annotate(
            order,
            xy=(1.0, last_value),
            xycoords=("axes fraction", "data"),
            xytext=(8, 0),
            textcoords="offset points",
            fontsize=13,
            va="center",
            color=color,
            fontweight="bold",
            clip_on=False,
        )

        # Add silhouette
        silhouette = get_silhouette(order)
        if silhouette is not None:
            rotation = config.get("rotation", 0)
            if rotation != 0:
                silhouette = silhouette.rotate(rotation, expand=True, resample=Image.BICUBIC)
            imagebox = OffsetImage(silhouette, zoom=config["zoom"])
            ab = AnnotationBbox(
                imagebox,
                (1.0, last_value),
                xycoords=("axes fraction", "data"),
                xybox=(config["x_offset"], config["y_offset"]),
                boxcoords="offset points",
                frameon=False,
                box_alignment=(0.5, 0.5),
                clip_on=False,
            )
            ax.add_artist(ab)
        elif order == "Others":
            # Display multiple silhouettes for "Others"
            for name in OTHERS_SILHOUETTES:
                cfg = silhouette_config.get(name, {"x_offset": 45, "y_offset": -30, "zoom": 0.025})
                other_img = load_silhouette(name)
                if other_img is not None:
                    rotation = cfg.get("rotation", 0)
                    if rotation != 0:
                        other_img = other_img.rotate(rotation, expand=True, resample=Image.BICUBIC)
                    imagebox = OffsetImage(other_img, zoom=cfg["zoom"])
                    ab = AnnotationBbox(
                        imagebox,
                        (1.0, last_value),
                        xycoords=("axes fraction", "data"),
                        xybox=(cfg["x_offset"], cfg["y_offset"]),
                        boxcoords="offset points",
                        frameon=False,
                        box_alignment=(0.5, 0.5),
                        clip_on=False,
                    )
                    ax.add_artist(ab)

    ax.set_title("Top Taxa Orders in ToxProt Over Time", fontsize=16)
    ax.set_xlabel("Year", fontsize=14)
    ax.set_ylabel("Number of Entries", fontsize=14)
    ax.set_xticks(plot_df.index[::2])
    ax.set_ylim(bottom=0)
    ax.tick_params(axis="both", labelsize=12)
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    plt.savefig(output_path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close()


def plot_taxa_newcomers(
    df_old: pd.DataFrame,
    df_new: pd.DataFrame,
    taxa_level: str,
    output_path: Path,
    color: str,
    max_items: int = 15,
):
    """Create horizontal bar chart of newcomer taxa."""
    newcomers = get_taxa_newcomers(df_old, df_new, taxa_level)

    if len(newcomers) == 0:
        return

    plot_data = newcomers.head(max_items)
    title_suffix = f" (Top {max_items} Shown)" if len(newcomers) > max_items else ""

    fig, ax = plt.subplots(figsize=(12, max(6, len(plot_data) * 0.5)))
    plot_data.plot(kind="barh", color=color, ax=ax, zorder=3)

    if taxa_level == "Order":
        labels = [f"{name} ({ORDER_COMMON_NAMES.get(name, '')})" for name in plot_data.index]
        ax.set_yticks(range(len(labels)))
        ax.set_yticklabels(labels, fontsize=11)

    ax.set_title(f"New Taxa {taxa_level}s in ToxProt 2025 (Not in 2017){title_suffix}", fontsize=16)
    ax.set_xlabel("Number of Entries", fontsize=14)
    ax.set_ylabel(taxa_level, fontsize=14)

    for i, v in enumerate(plot_data):
        ax.text(v + 0.5, i, str(v), va="center", fontsize=10)

    ax.grid(axis="x", linestyle="--", alpha=0.7, zorder=0)
    fig.tight_layout()

    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def plot_newcomers_alluvial(
    datasets: dict[int, pd.DataFrame],
    output_path: Path,
    taxa_level: str = "Phylum",
    years: list[int] | None = None,
    definition: str | None = DEFAULT_DEFINITION,
):
    if years is None:
        years = list(COMPARISON_YEARS)
    """Create alluvial-style diagram showing taxa at decade steps.

    Style matches definition_comparison.png:
    - Centered nodes
    - Only numbers inside nodes (on top for small nodes)
    - Names on the right side in appropriate colors
    - New taxa appear with appropriate height

    Args:
        datasets: Dictionary mapping year string to DataFrame
        output_path: Path to save the figure
        taxa_level: Taxonomic level to display ("Phylum", "Class", "Order")
        years: Years to show as columns (decade steps)
        definition: ToxProt definition filter
    """
    # Custom color palette (different from default)
    ALLUVIAL_COLORS = [
        "#4C78A8",  # Steel blue
        "#F58518",  # Orange
        "#54A24B",  # Green
        "#E45756",  # Red
        "#72B7B2",  # Teal
        "#EECA3B",  # Yellow
        "#B279A2",  # Purple
        "#FF9DA6",  # Pink
        "#9D755D",  # Brown
        "#79706E",  # Gray-brown
        "#D67195",  # Rose
        "#BAB0AC",  # Light gray
        "#1B9E77",  # Dark cyan
        "#D95F02",  # Dark orange
        "#7570B3",  # Slate blue
        "#E7298A",  # Magenta
        "#66A61E",  # Olive green
        "#E6AB02",  # Gold
        "#A6761D",  # Sienna
        "#666666",  # Dark gray
    ]

    # Get taxa counts for each year at specified level
    data_by_year = {}
    for year in years:
        if year not in datasets:
            continue
        df = filter_by_definition(datasets[year], definition)
        data_by_year[year] = df[taxa_level].value_counts().to_dict()

    # Get all unique taxa across all years, sorted by 2025 count
    reference_year = years[-1]
    all_taxa = set()
    for year in years:
        all_taxa.update(data_by_year[year].keys())

    # Sort by 2025 count (descending), taxa not in 2025 go at end
    taxa_order = sorted(all_taxa, key=lambda t: -data_by_year[reference_year].get(t, 0))

    # Create color map
    taxa_colors = {}
    for i, taxon in enumerate(taxa_order):
        taxa_colors[taxon] = ALLUVIAL_COLORS[i % len(ALLUVIAL_COLORS)]

    # Figure size depends on number of taxa
    n_taxa = len(all_taxa)
    fig_height = max(10, min(16, n_taxa * 0.5))
    fig, ax = plt.subplots(figsize=(14, fig_height))

    # Layout parameters
    n_years = len(years)
    x_positions = [i * 1.8 for i in range(n_years)]
    box_width = 0.4
    gap = 2  # Gap between stacked bars

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

    # Track which taxa are new at each year
    new_taxa_by_year = {}
    prev_taxa = set()
    for year in years:
        current_taxa = set(data_by_year[year].keys())
        new_taxa_by_year[year] = current_taxa - prev_taxa
        prev_taxa = current_taxa

    # Use 2025 total as reference for scaling
    reference_total = sum(data_by_year[reference_year].values())
    scale_factor = 100 / reference_total

    # Determine which taxa need labels on top (small bars)
    def needs_top_label(count):
        return count * scale_factor < 2.5

    # Calculate total heights for each year (for centering)
    # Add extra gap for taxa with top labels
    year_heights = {}
    for year in years:
        total_height = 0
        taxa_counts = data_by_year[year]
        n_taxa_year = len(taxa_counts)
        for _taxon, count in taxa_counts.items():
            total_height += count * scale_factor
        # Add gaps - extra gap for small bars with top labels
        n_small = sum(1 for c in taxa_counts.values() if needs_top_label(c))
        total_height += (n_taxa_year - 1) * gap + n_small * 1.5  # Extra space for top labels
        year_heights[year] = total_height

    max_height = max(year_heights.values())
    center_y = max_height / 2

    # Track positions for each year
    positions = {}

    # Draw bars for each year (centered)
    for i, year in enumerate(years):
        if year not in data_by_year:
            continue

        taxa_counts = data_by_year[year]
        x = x_positions[i]

        # Sort ALL taxa by 2025 order (consistent across all years)
        sorted_taxa = sorted(
            taxa_counts.items(),
            key=lambda x: taxa_order.index(x[0]) if x[0] in taxa_order else len(taxa_order),
        )

        # Calculate starting y offset for centering
        y_offset = center_y - year_heights[year] / 2
        positions[year] = {}
        positions[year]["right_labels"] = []  # For collecting labels

        for taxon, count in sorted_taxa:
            height = count * scale_factor
            is_small = needs_top_label(count)

            positions[year][taxon] = {
                "bottom": y_offset,
                "top": y_offset + height,
                "count": count,
                "height": height,
                "is_new": taxon in new_taxa_by_year[year],
            }

            color = taxa_colors.get(taxon, OTHER_COLOR)

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

            # Show count number - inside for large bars, on top for small bars
            if not is_small:
                fontsize = 11 if height > 5 else 9
                ax.text(
                    x,
                    y_offset + height / 2,
                    str(count),
                    ha="center",
                    va="center",
                    fontsize=fontsize,
                    fontweight="bold",
                    color="white",
                )
            else:
                # Small bar - put count on top with larger font
                ax.text(
                    x,
                    y_offset + height + 0.5,
                    str(count),
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    fontweight="bold",
                    color=color,
                )

            # Extra gap after small bars
            extra_gap = 1.5 if is_small else 0
            y_offset += height + gap + extra_gap

    # Draw flows between years (only for existing taxa, not new ones)
    for i in range(len(years) - 1):
        year1 = years[i]
        year2 = years[i + 1]

        if year1 not in positions or year2 not in positions:
            continue

        x1 = x_positions[i]
        x2 = x_positions[i + 1]

        # For each taxon that appears in both years (not new in year2)
        # Exclude the "right_labels" key which is used for label positioning
        taxa_year1 = {k for k in positions[year1].keys() if k != "right_labels"}
        taxa_year2 = {k for k in positions[year2].keys() if k != "right_labels"}
        common_taxa = taxa_year1 & taxa_year2

        for taxon in common_taxa:
            if positions[year2][taxon]["is_new"]:
                continue  # Don't draw flow to new taxa

            pos1 = positions[year1][taxon]
            pos2 = positions[year2][taxon]
            color = taxa_colors.get(taxon, OTHER_COLOR)

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

        # Add labels for taxa that disappear (in year1 but not in year2)
        disappeared_taxa = taxa_year1 - taxa_year2
        for taxon in disappeared_taxa:
            pos1 = positions[year1][taxon]
            color = taxa_colors.get(taxon, OTHER_COLOR)
            # Add label to the right of the bar (since there's no follow-up)
            ax.text(
                x1 + box_width / 2 + 0.1,
                (pos1["bottom"] + pos1["top"]) / 2,
                taxon,
                ha="left",
                va="center",
                fontsize=10,
                fontweight="bold",
                fontstyle="italic",
                color=color,
            )

    # Add taxa names on the right side of the last column only
    last_year = years[-1]
    last_x = x_positions[-1]

    for taxon in positions[last_year]:
        if taxon == "right_labels":
            continue
        pos = positions[last_year][taxon]
        color = taxa_colors.get(taxon, OTHER_COLOR)
        is_newcomer = taxon not in data_by_year[years[0]]

        ax.text(
            last_x + box_width / 2 + 0.1,
            (pos["bottom"] + pos["top"]) / 2,
            taxon,
            ha="left",
            va="center",
            fontsize=10,
            fontweight="bold",
            fontstyle="italic" if is_newcomer else "normal",
            color=color,
        )

    # Calculate totals for year labels
    totals = {year: sum(data_by_year[year].values()) for year in years if year in data_by_year}

    # Add year labels at bottom
    for i, year in enumerate(years):
        if year in totals:
            total = totals[year]
            ax.text(
                x_positions[i],
                -5,
                f"{year}\n(n={total:,})",
                ha="center",
                va="top",
                fontsize=14,
                fontweight="bold",
            )

    ax.set_xlim(-0.8, x_positions[-1] + 2.5)
    ax.set_ylim(-12, max_height + 5)
    ax.set_title(f"{taxa_level}-Level Taxa Evolution", fontsize=16, pad=20)
    ax.axis("off")

    plt.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.1)
    plt.close()


def main():
    """Generate all taxa analysis plots (for standalone execution)."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading ToxProt datasets...")
    datasets = load_datasets(YEARS)

    if len(datasets) < 2:
        print("Error: Need at least 2 datasets to generate plots")
        return

    print(f"Loaded {len(datasets)} datasets")

    print("Generating top taxa trend plot...")
    plot_top_taxa_trend(datasets, FIGURES_DIR / "top_taxa_trend")

    for level in ["Phylum", "Class", "Order", "Family"]:
        print(f"Generating newcomers alluvial plot ({level})...")
        plot_newcomers_alluvial(
            datasets,
            FIGURES_DIR / f"taxa_newcomers_alluvial_{level.lower()}.png",
            taxa_level=level,
        )

    print("Taxa analysis complete.")


if __name__ == "__main__":
    main()
