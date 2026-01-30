#!/usr/bin/env python3
"""
ToxProt Definition Comparison Analysis

Creates a two-panel figure comparing the ToxProt selection definitions:
- Panel A: Alluvial diagram showing taxonomic coverage
- Panel B: Venn diagram showing entry overlap
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib_venn import venn2

from .colors import DEFINITION_COLORS

# Display labels
CRITERIA_LABELS = {
    "venom_tissue": "Venom Tissue",
    "kw_toxin": "Toxin Keyword",
}

DEFINITION_LABELS = {
    "venom_tissue": "Venom Tissue",
    "kw_toxin": "Toxin Keyword",
    "both": "Both",
}


def get_criterion_subset(df: pd.DataFrame, criterion: str) -> pd.DataFrame:
    """Get entries matching a criterion (includes 'both' entries)."""
    if criterion == "venom_tissue":
        return df[df["ToxProt definition"].isin(["venom_tissue", "both"])]
    elif criterion == "kw_toxin":
        return df[df["ToxProt definition"].isin(["kw_toxin", "both"])]
    else:
        raise ValueError(f"Unknown criterion: {criterion}")


def get_definition_counts(df: pd.DataFrame) -> dict[str, int]:
    """Get counts for each ToxProt definition category."""
    counts = df["ToxProt definition"].value_counts().to_dict()
    return {
        "venom_tissue": counts.get("venom_tissue", 0),
        "kw_toxin": counts.get("kw_toxin", 0),
        "both": counts.get("both", 0),
    }


def get_taxa_exclusivity_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate unique taxa at each level that are shared vs exclusive."""
    results = []

    for level in ["Phylum", "Class", "Order", "Family"]:
        venom_df = get_criterion_subset(df, "venom_tissue")
        keyword_df = get_criterion_subset(df, "kw_toxin")

        venom_taxa = set(venom_df[level].dropna().unique())
        keyword_taxa = set(keyword_df[level].dropna().unique())

        shared = len(venom_taxa & keyword_taxa)
        venom_exclusive = len(venom_taxa - keyword_taxa)
        keyword_exclusive = len(keyword_taxa - venom_taxa)

        results.append(
            {
                "Level": level,
                "Shared": shared,
                "Venom": venom_exclusive,
                "Keyword": keyword_exclusive,
                "Total": shared + venom_exclusive + keyword_exclusive,
            }
        )

    return pd.DataFrame(results)


def get_taxa_by_exclusivity(df: pd.DataFrame, level: str) -> dict[str, set]:
    """Get taxa names grouped by exclusivity category."""
    venom_df = get_criterion_subset(df, "venom_tissue")
    keyword_df = get_criterion_subset(df, "kw_toxin")

    venom_taxa = set(venom_df[level].dropna().unique())
    keyword_taxa = set(keyword_df[level].dropna().unique())

    return {
        "Shared": venom_taxa & keyword_taxa,
        "Venom": venom_taxa - keyword_taxa,
        "Keyword": keyword_taxa - venom_taxa,
    }


def plot_alluvial_panel(df: pd.DataFrame, ax: plt.Axes) -> None:
    """Draw alluvial diagram showing taxonomic coverage."""
    counts_df = get_taxa_exclusivity_counts(df)

    levels = ["Phylum", "Order", "Family"]
    categories = ["Shared", "Venom", "Keyword"]
    colors = {
        "Shared": DEFINITION_COLORS["both"],
        "Venom": DEFINITION_COLORS["venom_tissue"],
        "Keyword": DEFINITION_COLORS["kw_toxin"],
    }

    # Get values as dict
    data = {}
    for _, row in counts_df.iterrows():
        data[row["Level"]] = {
            "Shared": row["Shared"],
            "Venom": row["Venom"],
            "Keyword": row["Keyword"],
        }

    phylum_taxa = get_taxa_by_exclusivity(df, "Phylum")

    # Layout parameters
    x_positions = [0, 1.2, 2.4]
    box_width = 0.3
    gap = 3

    def draw_alluvial_flow(x1, y1_bottom, y1_top, x2, y2_bottom, y2_top, color, alpha=0.4):
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

    # Normalize heights
    scale_factor = 100 / max(sum(data[level].values()) for level in levels)

    # Calculate total height for centering
    level_heights = {}
    for level in levels:
        total_height = 0
        for cat in categories:
            value = data[level][cat]
            total_height += value * scale_factor
            if value > 0:
                total_height += gap
        level_heights[level] = total_height - gap

    max_height = max(level_heights.values())
    center_y = max_height / 2

    positions = {level: {} for level in levels}

    for i, level in enumerate(levels):
        x = x_positions[i]
        y_offset = center_y - level_heights[level] / 2

        for cat in categories:
            value = data[level][cat]
            height = value * scale_factor

            positions[level][cat] = {
                "bottom": y_offset,
                "top": y_offset + height,
                "value": value,
            }

            if value > 0:
                rect = Rectangle(
                    (x - box_width / 2, y_offset),
                    box_width,
                    height,
                    facecolor=colors[cat],
                    edgecolor="black",
                    linewidth=1.5,
                )
                ax.add_patch(rect)

                ax.text(
                    x,
                    y_offset + height / 2,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=11,
                    fontweight="bold",
                    color="white",
                )

            y_offset += height + gap

    # Add phylum names
    phylum_x = x_positions[0]

    if phylum_taxa["Keyword"]:
        taxa_names = sorted(phylum_taxa["Keyword"])
        names_str = "\n".join(taxa_names)
        pos = positions["Phylum"]["Keyword"]
        ax.text(
            phylum_x,
            pos["top"] + 2,
            names_str,
            ha="center",
            va="bottom",
            fontsize=12,
            style="italic",
            color=colors["Keyword"],
        )

    if phylum_taxa["Shared"]:
        taxa_names = sorted(phylum_taxa["Shared"])
        names_str = "\n".join(taxa_names)
        pos = positions["Phylum"]["Shared"]
        ax.text(
            phylum_x,
            pos["bottom"] - 2,
            names_str,
            ha="center",
            va="top",
            fontsize=12,
            style="italic",
            color=colors["Shared"],
        )

    # Draw flows between levels
    for i in range(len(levels) - 1):
        level1 = levels[i]
        level2 = levels[i + 1]
        x1 = x_positions[i]
        x2 = x_positions[i + 1]

        for cat in categories:
            val1 = positions[level1][cat]["value"]
            val2 = positions[level2][cat]["value"]
            if val1 > 0 or val2 > 0:
                draw_alluvial_flow(
                    x1,
                    positions[level1][cat]["bottom"],
                    positions[level1][cat]["top"],
                    x2,
                    positions[level2][cat]["bottom"],
                    positions[level2][cat]["top"],
                    colors[cat],
                    alpha=0.4,
                )

    # Add level labels with totals
    for i, level in enumerate(levels):
        total = sum(data[level].values())
        ax.text(
            x_positions[i],
            -8,
            f"{level}\n(n={total})",
            ha="center",
            va="top",
            fontsize=12,
            fontweight="bold",
        )

    # Add legend
    legend_labels = {
        "Shared": "Both",
        "Venom": "Venom Tissue",
        "Keyword": "Toxin Keyword",
    }
    legend_order = ["Keyword", "Venom", "Shared"]
    legend_elements = [
        plt.Rectangle(
            (0, 0), 1, 1, facecolor=colors[cat], edgecolor="black", label=legend_labels[cat]
        )
        for cat in legend_order
    ]
    ax.legend(
        handles=legend_elements,
        loc="upper left",
        fontsize=12,
        framealpha=0.9,
    )

    # Formatting
    ax.set_xlim(-0.3, 2.7)
    all_bottoms = [positions[level][cat]["bottom"] for level in levels for cat in categories]
    all_tops = [positions[level][cat]["top"] for level in levels for cat in categories]
    y_min = min(all_bottoms)
    y_max = max(all_tops)
    ax.set_ylim(y_min - 1, y_max + 1)
    ax.axis("off")
    ax.set_title("Taxonomic Coverage", fontsize=14, fontweight="bold", pad=10)


def plot_venn_panel(df: pd.DataFrame, ax: plt.Axes) -> None:
    """Draw Venn diagram showing entry overlap."""
    counts = get_definition_counts(df)

    venn = venn2(
        subsets=(counts["venom_tissue"], counts["kw_toxin"], counts["both"]),
        set_labels=("", ""),
        ax=ax,
    )

    # Style circles (no border)
    if venn.get_patch_by_id("10"):
        venn.get_patch_by_id("10").set_color(DEFINITION_COLORS["venom_tissue"])
        venn.get_patch_by_id("10").set_alpha(0.7)
        venn.get_patch_by_id("10").set_edgecolor("none")
    if venn.get_patch_by_id("01"):
        venn.get_patch_by_id("01").set_color(DEFINITION_COLORS["kw_toxin"])
        venn.get_patch_by_id("01").set_alpha(0.7)
        venn.get_patch_by_id("01").set_edgecolor("none")
    if venn.get_patch_by_id("11"):
        venn.get_patch_by_id("11").set_color(DEFINITION_COLORS["both"])
        venn.get_patch_by_id("11").set_alpha(0.7)
        venn.get_patch_by_id("11").set_edgecolor("none")

    # Style numbers
    for text in venn.subset_labels:
        if text:
            text.set_fontsize(14)
            text.set_fontweight("bold")

    # Add colored labels below
    ax.text(
        0.25,
        -0.05,
        "Venom Tissue",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        color=DEFINITION_COLORS["venom_tissue"],
        ha="center",
    )
    ax.text(
        0.75,
        -0.05,
        "Toxin Keyword",
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        color=DEFINITION_COLORS["kw_toxin"],
        ha="center",
    )

    ax.set_title("Entry Overlap", fontsize=14, fontweight="bold", pad=10)


def create_definition_comparison_figure(df: pd.DataFrame, output_dir: Path) -> None:
    """Create combined two-panel figure with alluvial and Venn diagrams."""
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={"width_ratios": [1.5, 1]})

    plot_alluvial_panel(df, ax_a)
    ax_a.text(-0.05, 1.05, "A", transform=ax_a.transAxes, fontsize=18, fontweight="bold")

    plot_venn_panel(df, ax_b)
    ax_b.text(-0.1, 1.05, "B", transform=ax_b.transAxes, fontsize=18, fontweight="bold")

    fig.suptitle(
        "ToxProt Selection Definition Comparison (2025)",
        fontsize=16,
        fontweight="bold",
        y=1.02,
    )

    plt.tight_layout()
    plt.savefig(output_dir / "definition_comparison.png", dpi=300, bbox_inches="tight")
    plt.savefig(output_dir / "definition_comparison.pdf", bbox_inches="tight")
    plt.close()


def main():
    """Main analysis pipeline."""
    data_path = Path("data/processed/toxprot/toxprot_2025.csv")
    output_dir = Path("figures/definitions")

    df = pd.read_csv(data_path)
    print(f"Loaded {len(df):,} entries")

    create_definition_comparison_figure(df, output_dir)
    print(f"Saved to {output_dir}")


if __name__ == "__main__":
    main()
