#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import os

# --- Configuration ---
DATA_DIR = "data/processed/"
FIGURES_DIR = "figures/taxa/"
TOXPROT_2017_FILE = os.path.join(DATA_DIR, "toxprot_2017.csv")
TOXPROT_2025_FILE = os.path.join(DATA_DIR, "toxprot_2025.csv")

# Ensure the figures directory exists
os.makedirs(FIGURES_DIR, exist_ok=True)

# --- Function Definitions ---


def get_top_taxa_counts(df, taxa_column="Order", top_n=5):
    """Function to get top N taxa and group the rest as "Others"."""
    taxa_counts = df[taxa_column].value_counts()
    top_taxa = taxa_counts.nlargest(top_n)
    result = top_taxa.to_dict()
    others_count = taxa_counts[~taxa_counts.index.isin(top_taxa.index)].sum()
    result["Others"] = others_count
    return result


def get_taxa_newcomers(df_2017, df_2025, taxa_level="Order"):
    """Identify taxa present in 2025 dataset but not in 2017 dataset."""
    taxa_2017_set = set(df_2017[taxa_level].unique())
    taxa_2025_set = set(df_2025[taxa_level].unique())
    newcomers = taxa_2025_set - taxa_2017_set
    return df_2025[df_2025[taxa_level].isin(newcomers)][taxa_level].value_counts()


def main():
    """Main function to load data and generate plots."""
    # Load the data
    try:
        toxprot_2017 = pd.read_csv(TOXPROT_2017_FILE)
        toxprot_2025 = pd.read_csv(TOXPROT_2025_FILE)
    except FileNotFoundError as e:
        print(f"Error loading data: {e}")
        print(
            f"Please ensure the files '{TOXPROT_2017_FILE}' and '{TOXPROT_2025_FILE}' exist."
        )
        return

    # --- Plot 1: Top 5 Taxa Distribution ---
    print("Generating Top 5 Taxa Distribution plot...")
    top_taxa_2017_counts = get_top_taxa_counts(toxprot_2017, "Order", 5)
    top_taxa_2025_counts = get_top_taxa_counts(toxprot_2025, "Order", 5)

    plot_data_top_taxa = pd.DataFrame(
        {"2017": top_taxa_2017_counts, "2025": top_taxa_2025_counts}
    )
    plot_data_top_taxa = plot_data_top_taxa.fillna(0)
    plot_data_top_taxa = plot_data_top_taxa.T

    fig1, ax_top_taxa = plt.subplots(figsize=(10, 8))
    plot_data_top_taxa.plot(
        kind="bar", stacked=True, colormap="tab10", ax=ax_top_taxa, zorder=3
    )

    ax_top_taxa.set_title(
        "Distribution of Top 5 Taxa Orders in ToxProt 2017 vs 2025", fontsize=20
    )
    ax_top_taxa.set_xlabel("Year", fontsize=18)
    ax_top_taxa.set_ylabel("Count", fontsize=18)
    ax_top_taxa.tick_params(axis="x", rotation=0, labelsize=16)
    ax_top_taxa.tick_params(axis="y", labelsize=16)  # Added y-axis tick label size
    # Get the legend handles and labels, then reverse them
    handles, labels = ax_top_taxa.get_legend_handles_labels()
    ax_top_taxa.legend(
        handles[::-1],
        labels[::-1],  # Reverse the order
        title="Order",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        fontsize=16,
        title_fontsize=18,
    )

    for container in ax_top_taxa.containers:
        ax_top_taxa.bar_label(
            container, label_type="center", fmt="%d", fontsize=16
        )  # Added fontsize to bar labels

    ax_top_taxa.grid(axis="y", linestyle="--", alpha=0.7, zorder=0)
    fig1.tight_layout()
    plt.savefig(
        os.path.join(FIGURES_DIR, "top_taxa_distribution.png"),
        dpi=300,
        bbox_inches="tight",
    )
    print(f"Saved {os.path.join(FIGURES_DIR, 'top_taxa_distribution.png')}")

    # --- Plot 2: Newcomer Taxa Orders ---
    print("Generating Newcomer Taxa Orders plot...")
    order_common_names = {
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

    order_newcomers = get_taxa_newcomers(toxprot_2017, toxprot_2025, "Order")

    fig2, ax_order_newcomers = plt.subplots(figsize=(12, 8))
    plot_data_order_newcomers = (
        order_newcomers.head(15) if len(order_newcomers) > 15 else order_newcomers
    )
    title_suffix_order = " (Top 15 Shown)" if len(order_newcomers) > 15 else ""

    plot_data_order_newcomers.plot(
        kind="barh", color="tab:green", ax=ax_order_newcomers, zorder=3
    )

    labels_order = [
        f"{order} ({order_common_names.get(order, 'Unknown')})"
        for order in plot_data_order_newcomers.index
    ]
    ax_order_newcomers.set_yticks(range(len(labels_order)))
    ax_order_newcomers.set_yticklabels(labels_order, fontsize=11)

    ax_order_newcomers.set_title(
        f"New Taxa Orders in ToxProt 2025 (Not Present in 2017){title_suffix_order}",
        fontsize=16,
    )
    ax_order_newcomers.set_xlabel("Number of Entries", fontsize=14)
    ax_order_newcomers.set_ylabel("Order (Common Name)", fontsize=14)

    for i, v_val in enumerate(
        plot_data_order_newcomers
    ):  # Renamed v to v_val to avoid confusion
        ax_order_newcomers.text(v_val + 0.5, i, str(v_val), va="center", fontsize=10)

    ax_order_newcomers.grid(axis="x", linestyle="--", alpha=0.7, zorder=0)
    fig2.tight_layout()
    plt.savefig(
        os.path.join(FIGURES_DIR, "taxa_newcomers_order.png"),
        dpi=300,
        bbox_inches="tight",
    )
    print(f"Saved {os.path.join(FIGURES_DIR, 'taxa_newcomers_order.png')}")

    # --- Plot 3: Newcomer Taxa Families ---
    print("Generating Newcomer Taxa Families plot...")
    family_newcomers = get_taxa_newcomers(toxprot_2017, toxprot_2025, "Family")

    fig3, ax_family_newcomers = plt.subplots(figsize=(12, 10))
    plot_data_family_newcomers = (
        family_newcomers.head(15) if len(family_newcomers) > 15 else family_newcomers
    )
    title_suffix_family = " (Top 15 Shown)" if len(family_newcomers) > 15 else ""

    plot_data_family_newcomers.plot(
        kind="barh", color="tab:blue", ax=ax_family_newcomers, zorder=3
    )
    ax_family_newcomers.set_title(
        f"New Taxa Families in ToxProt 2025 (Not Present in 2017){title_suffix_family}",
        fontsize=16,
    )
    ax_family_newcomers.set_xlabel("Number of Entries", fontsize=14)
    ax_family_newcomers.set_ylabel("Family", fontsize=14)
    ax_family_newcomers.tick_params(axis="y", labelsize=12)

    for i, v_val in enumerate(plot_data_family_newcomers):  # Renamed v to v_val
        ax_family_newcomers.text(v_val + 0.5, i, str(v_val), va="center", fontsize=10)

    ax_family_newcomers.grid(axis="x", linestyle="--", alpha=0.7, zorder=0)
    fig3.tight_layout()
    plt.savefig(
        os.path.join(FIGURES_DIR, "taxa_newcomers_family.png"),
        dpi=300,
        bbox_inches="tight",
    )
    print(f"Saved {os.path.join(FIGURES_DIR, 'taxa_newcomers_family.png')}")

    print("\nAll plots generated and saved.")


if __name__ == "__main__":
    main()
