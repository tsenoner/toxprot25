#!/usr/bin/env python3
"""Analyze and visualize taxonomic distribution in ToxProt datasets."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from PIL import Image

from .colors import TAXA_ORDER_COLORS, TAXA_ORDER_COLORS_LIST

# --- Configuration ---
DATA_DIR = Path("data/processed/toxprot")
FIGURES_DIR = Path("figures/taxa")
SILHOUETTE_DIR = Path("data/raw/PhyloPic/png")
YEARS = [str(y) for y in range(2005, 2026)]

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


def load_datasets(years: list[str], data_dir: Path = DATA_DIR) -> dict[str, pd.DataFrame]:
    """Load ToxProt datasets for specified years."""
    datasets = {}
    for year in years:
        filepath = data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)
    return datasets


def get_taxa_newcomers(
    df_old: pd.DataFrame, df_new: pd.DataFrame, taxa_level: str = "Order"
) -> pd.Series:
    """Identify taxa present in new dataset but not in old dataset."""
    old_taxa = set(df_old[taxa_level].unique())
    new_taxa = set(df_new[taxa_level].unique())
    newcomers = new_taxa - old_taxa
    return df_new[df_new[taxa_level].isin(newcomers)][taxa_level].value_counts()


def plot_top_taxa_trend(
    datasets: dict[str, pd.DataFrame],
    output_path: Path,
    reference_year: str = "2025",
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


def main():
    """Generate all taxa analysis plots."""
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    datasets = load_datasets(YEARS)

    if len(datasets) < 2:
        print("Error: Need at least 2 datasets to generate plots")
        return

    plot_top_taxa_trend(datasets, FIGURES_DIR / "top_taxa_trend")

    if "2017" in datasets and "2025" in datasets:
        plot_taxa_newcomers(
            datasets["2017"],
            datasets["2025"],
            "Order",
            FIGURES_DIR / "taxa_newcomers_order.png",
            "tab:green",
        )
        plot_taxa_newcomers(
            datasets["2017"],
            datasets["2025"],
            "Family",
            FIGURES_DIR / "taxa_newcomers_family.png",
            "tab:blue",
        )

    print("Taxa analysis complete.")


if __name__ == "__main__":
    main()
