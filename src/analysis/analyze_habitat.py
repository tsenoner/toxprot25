"""Analyze protein family distributions by habitat (terrestrial vs marine).

Generates visualizations comparing dual-habitat protein families across
ToxProt datasets from different years (2005, 2015, 2025).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, PathPatch
from matplotlib.path import Path as MplPath

from ..config import COMPARISON_YEARS
from .analyze_protein_families import FAMILY_NAME_MAP

# --- Configuration ---
TOP_N_FAMILIES = 15
YEARS = COMPARISON_YEARS


def normalize_family_name(name: str) -> str:
    """Normalize a family name using the mapping dictionary."""
    if pd.isna(name) or name == "":
        return name
    return FAMILY_NAME_MAP.get(name, name)


def get_protein_family_by_habitat(df: pd.DataFrame, habitat_type: str) -> pd.Series:
    """Extract protein family counts for a specific habitat type."""
    habitat_data = df[df["Habitat"] == habitat_type].copy()
    habitat_data["Protein families"] = habitat_data["Protein families"].fillna("").astype(str)
    habitat_data["Protein families"] = habitat_data["Protein families"].apply(normalize_family_name)
    return (
        habitat_data[habitat_data["Protein families"] != ""]["Protein families"]
        .value_counts()
        .sort_index()
    )


def get_dual_habitat_families(datasets: dict[int, pd.DataFrame]) -> set[str]:
    """Identify protein families present in both terrestrial and marine habitats."""
    all_terrestrial = set()
    all_marine = set()

    for df in datasets.values():
        pf_terrestrial = get_protein_family_by_habitat(df, "terrestrial")
        pf_marine = get_protein_family_by_habitat(df, "marine")
        all_terrestrial.update(pf_terrestrial.index)
        all_marine.update(pf_marine.index)

    return all_terrestrial & all_marine


def filter_df_by_families(df: pd.DataFrame, families_to_keep: set[str]) -> pd.DataFrame:
    """Filter DataFrame to keep only rows belonging to specified protein families."""
    if not families_to_keep:
        return pd.DataFrame(columns=df.columns)
    temp_families_col = (
        df["Protein families"]
        .fillna("__NAN_PLACEHOLDER__")
        .astype(str)
        .apply(normalize_family_name)
    )
    mask = temp_families_col.isin(families_to_keep)
    return df[mask].copy()


def plot_habitat_protein_families(
    datasets: dict[int, pd.DataFrame],
    ax: plt.Axes,
    top_n: int = TOP_N_FAMILIES,
) -> None:
    """Plot dual-habitat protein family evolution across years.

    Args:
        datasets: Dictionary mapping years to DataFrames
        ax: Axes to draw on.
        top_n: Number of top families to display
    """
    years = sorted(datasets.keys())
    if len(years) < 3:
        raise ValueError("Need at least 3 years of data")

    year_early, year_mid, year_late = years[0], years[1], years[-1]

    # Identify dual-habitat families
    dual_habitat_families = get_dual_habitat_families(datasets)
    if not dual_habitat_families:
        print("No dual-habitat protein families found.")
        return

    # Filter datasets
    filtered_datasets = {
        year: filter_df_by_families(df, dual_habitat_families) for year, df in datasets.items()
    }

    # Get protein family counts by habitat for each year
    pf_data = {}
    for year, df in filtered_datasets.items():
        pf_data[f"terrestrial_{year}"] = get_protein_family_by_habitat(df, "terrestrial")
        pf_data[f"marine_{year}"] = get_protein_family_by_habitat(df, "marine")

    # Calculate total counts for ranking
    total_counts = pd.Series(dtype=float)
    for counts in pf_data.values():
        total_counts = total_counts.add(counts, fill_value=0)
    top_families = total_counts.sort_values(ascending=False).head(top_n).index.tolist()

    # Build summary DataFrame
    summary = pd.DataFrame(index=top_families)
    for key, counts in pf_data.items():
        summary[key] = counts.reindex(top_families, fill_value=0)

    # Sort by total representatives in latest year (most on top)
    summary["total_latest"] = summary[f"terrestrial_{year_late}"] + summary[f"marine_{year_late}"]
    summary = summary.sort_values(by="total_latest", ascending=False)

    y_pos = np.arange(len(summary))

    # Colors for the three periods (lighter to darker)
    colors = {
        "terrestrial": ["#90EE90", "#31c42f", "forestgreen"],
        "marine": ["#ADD8E6", "#56a0dd", "steelblue"],
    }

    def abs_formatter(x, _pos):
        return f"{abs(x):.0f}"

    # Get values for each year
    t_early = summary[f"terrestrial_{year_early}"]
    t_mid = summary[f"terrestrial_{year_mid}"]
    t_late = summary[f"terrestrial_{year_late}"]
    m_early = summary[f"marine_{year_early}"]
    m_mid = summary[f"marine_{year_mid}"]
    m_late = summary[f"marine_{year_late}"]

    # Terrestrial bars (right side)
    t_base = np.minimum(t_early, t_late)
    t_mid_cumulative = np.minimum(t_mid, t_late)
    t_added_early_to_mid = np.maximum(0, t_mid_cumulative - t_base)
    t_added_mid_to_late = np.maximum(0, t_late - t_mid_cumulative)

    ax.barh(y_pos, t_base, color=colors["terrestrial"][0], label=str(year_early))
    ax.barh(
        y_pos,
        t_added_early_to_mid,
        left=t_base,
        color=colors["terrestrial"][1],
        label=str(year_mid),
    )
    ax.barh(
        y_pos,
        t_added_mid_to_late,
        left=t_base + t_added_early_to_mid,
        color=colors["terrestrial"][2],
        label=str(year_late),
    )

    # Marine bars (left side, negative values)
    m_base = np.minimum(m_early, m_late)
    m_mid_cumulative = np.minimum(m_mid, m_late)
    m_added_early_to_mid = np.maximum(0, m_mid_cumulative - m_base)
    m_added_mid_to_late = np.maximum(0, m_late - m_mid_cumulative)

    ax.barh(y_pos, -m_base, color=colors["marine"][0])
    ax.barh(y_pos, -m_added_early_to_mid, left=-m_base, color=colors["marine"][1])
    ax.barh(
        y_pos,
        -m_added_mid_to_late,
        left=-(m_base + m_added_early_to_mid),
        color=colors["marine"][2],
    )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(summary.index)
    ax.invert_yaxis()
    ax.set_xlabel("Number of Entries (Marine ← | → Terrestrial)")

    # Custom legend with Marine/Terrestrial columns and years
    legend_elements = [
        # Marine column (left)
        Patch(facecolor=colors["marine"][0], edgecolor="none", label=str(year_early)),
        Patch(facecolor=colors["marine"][1], edgecolor="none", label=str(year_mid)),
        Patch(facecolor=colors["marine"][2], edgecolor="none", label=str(year_late)),
        # Terrestrial column (right)
        Patch(facecolor=colors["terrestrial"][0], edgecolor="none", label=str(year_early)),
        Patch(facecolor=colors["terrestrial"][1], edgecolor="none", label=str(year_mid)),
        Patch(facecolor=colors["terrestrial"][2], edgecolor="none", label=str(year_late)),
    ]
    ax.legend(
        handles=legend_elements,
        loc="lower right",
        fontsize=10,
        ncol=2,
        columnspacing=1.5,
        handletextpad=0.5,
        handlelength=1.5,
        title="Marine      Terrestrial",
        title_fontproperties={"weight": "bold", "size": 10},
    )
    ax.axvline(0, color="grey", lw=0.8)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(abs_formatter))

    # Expand x-limits for annotations (balance space for annotations on both sides)
    xlim = ax.get_xlim()
    ax.set_xlim(xlim[0] * 2.3, xlim[1] * 1.1)
    xlim = ax.get_xlim()
    offset = abs(xlim[1]) * 0.02

    # Add annotations showing evolution: early→mid→late
    for i, family in enumerate(summary.index):
        t_e, t_m, t_l = int(t_early[family]), int(t_mid[family]), int(t_late[family])
        m_e, m_m, m_l = int(m_early[family]), int(m_mid[family]), int(m_late[family])

        # Terrestrial annotation (right side, at bar tip)
        ax.text(
            t_l + offset,
            y_pos[i],
            f"{t_e}→{t_m}→{t_l}",
            va="center",
            ha="left",
            fontsize=8,
        )
        # Marine annotation (left side, at bar tip)
        ax.text(
            -m_l - offset,
            y_pos[i],
            f"{m_e}→{m_m}→{m_l}",
            va="center",
            ha="right",
            fontsize=8,
        )


def plot_taxa_by_habitat(
    df: pd.DataFrame,
    ax: plt.Axes,
) -> None:
    """Create single-panel figure showing taxa distribution by habitat with flow connections.

    Shows stacked bars for Entries, Species, and Protein Families by habitat,
    with the Protein Families bar "expanding" into Shared/Exclusive breakdown
    via connecting flow lines.

    Args:
        df: DataFrame with columns: Entry, Species, Habitat, Protein families.
            Should already be filtered to venom_tissue criterion.
        ax: Axes to draw on.
    """
    # Filter to venom_tissue criterion
    venom_df = df[df["ToxProt definition"].isin(["venom_tissue", "both"])].copy()

    if venom_df.empty:
        print("No data matching venom_tissue criterion.")
        return

    # --- Data preparation ---
    # Habitat-level summary
    habitat_summary = (
        venom_df.groupby("Habitat")
        .agg(
            entry_count=("Entry", "count"),
            species_count=("Species", "nunique"),
            family_count=("Protein families", "nunique"),
        )
        .reset_index()
    )

    # === Identify shared protein families ===
    habitats = ["terrestrial", "marine"]
    terr_families = set(
        venom_df[venom_df["Habitat"] == "terrestrial"]["Protein families"].dropna().unique()
    )
    marine_families = set(
        venom_df[venom_df["Habitat"] == "marine"]["Protein families"].dropna().unique()
    )
    shared_families = terr_families & marine_families
    n_shared = len(shared_families)

    # Get counts per habitat
    data_by_habitat = {}
    for hab in habitats:
        hab_data = habitat_summary[habitat_summary["Habitat"] == hab]
        if len(hab_data) > 0:
            data_by_habitat[hab] = {
                "species": hab_data["species_count"].values[0],
                "entries": hab_data["entry_count"].values[0],
                "families": hab_data["family_count"].values[0],
            }
        else:
            data_by_habitat[hab] = {"species": 0, "entries": 0, "families": 0}

    # Calculate totals
    total_species = sum(d["species"] for d in data_by_habitat.values())
    total_entries = sum(d["entries"] for d in data_by_habitat.values())
    total_families = len(terr_families | marine_families)  # Union of all families

    # Exclusive family counts
    n_terr_exclusive = len(terr_families - shared_families)
    n_marine_exclusive = len(marine_families - shared_families)

    # Use 2005 colors from dual_habitat_protein_families figure (darker, more readable)
    habitat_colors = {
        "terrestrial": "forestgreen",  # Dark green (2005)
        "marine": "steelblue",  # Dark blue (2005)
    }

    # === Layout configuration ===
    # Left side: Entries, Species, Protein Families (3 bars)
    # Right side: Stacked vertically - Shared (top 50%), Exclusive (bottom 50%)
    metrics_left = ["Entries", "Species", "Protein\nFamilies"]

    x_left = np.array([0, 1, 2])
    x_right = np.array([3.7, 4.3])  # Just Entries, Species columns (closer together)
    bar_width = 0.5

    # === Calculate percentages for left-side bars ===
    # Order: Entries, Species, Protein Families
    terr_pcts_left = [
        data_by_habitat["terrestrial"]["entries"] / total_entries * 100 if total_entries > 0 else 0,
        data_by_habitat["terrestrial"]["species"] / total_species * 100 if total_species > 0 else 0,
        n_terr_exclusive / total_families * 100 if total_families > 0 else 0,
    ]
    marine_pcts_left = [
        data_by_habitat["marine"]["entries"] / total_entries * 100 if total_entries > 0 else 0,
        data_by_habitat["marine"]["species"] / total_species * 100 if total_species > 0 else 0,
        n_marine_exclusive / total_families * 100 if total_families > 0 else 0,
    ]
    shared_pcts_left = [
        0,  # Entries can't be "shared" in same way
        0,  # Species can't be "shared" in same way
        n_shared / total_families * 100 if total_families > 0 else 0,
    ]

    # Counts for labels
    terr_counts_left = [
        data_by_habitat["terrestrial"]["entries"],
        data_by_habitat["terrestrial"]["species"],
        n_terr_exclusive,
    ]
    marine_counts_left = [
        data_by_habitat["marine"]["entries"],
        data_by_habitat["marine"]["species"],
        n_marine_exclusive,
    ]

    # === Draw left-side stacked bars ===
    # Terrestrial (bottom)
    ax.bar(
        x_left,
        terr_pcts_left,
        bar_width,
        label="Terrestrial",
        color=habitat_colors["terrestrial"],
        edgecolor="black",
        linewidth=0.5,
    )

    # Shared families (middle, only for Protein Families)
    ax.bar(
        x_left,
        shared_pcts_left,
        bar_width,
        bottom=terr_pcts_left,
        label=f"Shared (n={n_shared})",
        color="#7570b3",  # Purple
        edgecolor="black",
        linewidth=0.5,
    )

    # Marine (top)
    ax.bar(
        x_left,
        marine_pcts_left,
        bar_width,
        bottom=[t + s for t, s in zip(terr_pcts_left, shared_pcts_left, strict=True)],
        label="Marine",
        color=habitat_colors["marine"],
        edgecolor="black",
        linewidth=0.5,
    )

    # === Calculate data for right-side bars (Shared/Exclusive by Entries/Species) ===
    venom_df.loc[:, "is_shared"] = venom_df["Protein families"].isin(shared_families)

    panel_right_data = {}
    for hab in habitats:
        hab_df = venom_df[venom_df["Habitat"] == hab]
        panel_right_data[hab] = {
            "shared_entries": hab_df[hab_df["is_shared"]]["Entry"].count(),
            "exclusive_entries": hab_df[~hab_df["is_shared"]]["Entry"].count(),
            "shared_species": hab_df[hab_df["is_shared"]]["Species"].nunique(),
            "exclusive_species": hab_df[~hab_df["is_shared"]]["Species"].nunique(),
        }

    # Calculate totals for each metric (Entries, Species)
    shared_totals = [
        panel_right_data["terrestrial"]["shared_entries"]
        + panel_right_data["marine"]["shared_entries"],
        panel_right_data["terrestrial"]["shared_species"]
        + panel_right_data["marine"]["shared_species"],
    ]
    exclusive_totals = [
        panel_right_data["terrestrial"]["exclusive_entries"]
        + panel_right_data["marine"]["exclusive_entries"],
        panel_right_data["terrestrial"]["exclusive_species"]
        + panel_right_data["marine"]["exclusive_species"],
    ]

    # Shared: Terrestrial counts and percentages (scaled to 50% height)
    shared_terr_counts = [
        panel_right_data["terrestrial"]["shared_entries"],
        panel_right_data["terrestrial"]["shared_species"],
    ]
    shared_terr_pcts = [
        (t / total * 50) if total > 0 else 0
        for t, total in zip(shared_terr_counts, shared_totals, strict=True)
    ]
    shared_marine_counts = [
        panel_right_data["marine"]["shared_entries"],
        panel_right_data["marine"]["shared_species"],
    ]
    shared_marine_pcts = [
        (m / total * 50) if total > 0 else 0
        for m, total in zip(shared_marine_counts, shared_totals, strict=True)
    ]

    # Exclusive: Terrestrial counts and percentages (scaled to 50% height)
    exclusive_terr_counts = [
        panel_right_data["terrestrial"]["exclusive_entries"],
        panel_right_data["terrestrial"]["exclusive_species"],
    ]
    exclusive_terr_pcts = [
        (t / total * 50) if total > 0 else 0
        for t, total in zip(exclusive_terr_counts, exclusive_totals, strict=True)
    ]
    exclusive_marine_counts = [
        panel_right_data["marine"]["exclusive_entries"],
        panel_right_data["marine"]["exclusive_species"],
    ]
    exclusive_marine_pcts = [
        (m / total * 50) if total > 0 else 0
        for m, total in zip(exclusive_marine_counts, exclusive_totals, strict=True)
    ]

    # === Draw right-side bars (vertically stacked with gap) ===
    # Layout: Exclusive (0-45), gap (45-55), Shared (55-100)
    gap_bottom = 45
    gap_top = 55
    exclusive_height = gap_bottom  # 0-45
    shared_height = 100 - gap_top  # 55-100 = 45

    # Scale percentages to fit in their respective sections
    exclusive_terr_scaled = [(t / 50) * exclusive_height for t in exclusive_terr_pcts]
    exclusive_marine_scaled = [(m / 50) * exclusive_height for m in exclusive_marine_pcts]
    shared_terr_scaled = [(t / 50) * shared_height for t in shared_terr_pcts]
    shared_marine_scaled = [(m / 50) * shared_height for m in shared_marine_pcts]

    # Exclusive group (bottom: y=0-45)
    ax.bar(
        x_right,
        exclusive_terr_scaled,
        bar_width,
        bottom=0,
        color=habitat_colors["terrestrial"],
        edgecolor="black",
        linewidth=0.5,
    )
    ax.bar(
        x_right,
        exclusive_marine_scaled,
        bar_width,
        bottom=exclusive_terr_scaled,
        color=habitat_colors["marine"],
        edgecolor="black",
        linewidth=0.5,
    )

    # Shared group (top: y=55-100)
    ax.bar(
        x_right,
        shared_terr_scaled,
        bar_width,
        bottom=gap_top,
        color=habitat_colors["terrestrial"],
        edgecolor="black",
        linewidth=0.5,
    )
    ax.bar(
        x_right,
        shared_marine_scaled,
        bar_width,
        bottom=[gap_top + t for t in shared_terr_scaled],
        color=habitat_colors["marine"],
        edgecolor="black",
        linewidth=0.5,
    )

    # === Add value labels inside left-side bars ===
    for i in range(len(x_left)):
        # Terrestrial label
        t_pct = terr_pcts_left[i]
        ax.text(
            x_left[i],
            t_pct / 2,
            f"{terr_counts_left[i]:,}\n({t_pct:.1f}%)",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            fontweight="bold",
        )

        # Shared label (only for Protein Families bar)
        s_pct = shared_pcts_left[i]
        if s_pct > 0:
            ax.text(
                x_left[i],
                terr_pcts_left[i] + s_pct / 2,
                f"{n_shared}\n({s_pct:.1f}%)",
                ha="center",
                va="center",
                fontsize=7,
                color="white",
                fontweight="bold",
            )

        # Marine label
        m_pct = marine_pcts_left[i]
        bottom_m = terr_pcts_left[i] + shared_pcts_left[i]
        ax.text(
            x_left[i],
            bottom_m + m_pct / 2,
            f"{marine_counts_left[i]:,}\n({m_pct:.1f}%)",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            fontweight="bold",
        )

    # === Add value labels inside right-side bars ===
    # Exclusive group labels (bottom: 0-45)
    for i in range(len(x_right)):
        # Terrestrial
        t_scaled = exclusive_terr_scaled[i]
        t_pct_actual = (
            exclusive_terr_counts[i] / exclusive_totals[i] * 100 if exclusive_totals[i] > 0 else 0
        )
        ax.text(
            x_right[i],
            t_scaled / 2,
            f"{exclusive_terr_counts[i]:,}\n({t_pct_actual:.1f}%)",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            fontweight="bold",
        )
        # Marine
        m_scaled = exclusive_marine_scaled[i]
        m_pct_actual = (
            exclusive_marine_counts[i] / exclusive_totals[i] * 100 if exclusive_totals[i] > 0 else 0
        )
        ax.text(
            x_right[i],
            exclusive_terr_scaled[i] + m_scaled / 2,
            f"{exclusive_marine_counts[i]:,}\n({m_pct_actual:.1f}%)",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            fontweight="bold",
        )

    # Shared group labels (top: 55-100)
    min_height_for_label = 8  # Minimum scaled height to fit label inside
    for i in range(len(x_right)):
        # Terrestrial
        t_scaled = shared_terr_scaled[i]
        t_pct_actual = shared_terr_counts[i] / shared_totals[i] * 100 if shared_totals[i] > 0 else 0
        ax.text(
            x_right[i],
            gap_top + t_scaled / 2,
            f"{shared_terr_counts[i]:,}\n({t_pct_actual:.1f}%)",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            fontweight="bold",
        )
        # Marine - place outside if section too small
        m_scaled = shared_marine_scaled[i]
        m_pct_actual = (
            shared_marine_counts[i] / shared_totals[i] * 100 if shared_totals[i] > 0 else 0
        )
        if m_scaled >= min_height_for_label:
            ax.text(
                x_right[i],
                gap_top + shared_terr_scaled[i] + m_scaled / 2,
                f"{shared_marine_counts[i]:,}\n({m_pct_actual:.1f}%)",
                ha="center",
                va="center",
                fontsize=7,
                color="white",
                fontweight="bold",
            )
        else:
            # Place label outside, directly on top of the bar
            bar_top = gap_top + shared_terr_scaled[i] + m_scaled
            ax.text(
                x_right[i],
                bar_top + 1,
                f"{shared_marine_counts[i]:,} ({m_pct_actual:.1f}%)",
                ha="center",
                va="bottom",
                fontsize=7,
                fontweight="bold",
                color=habitat_colors["marine"],
            )

    # === Draw flow connections using Bézier curves ===
    def draw_flow(
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
        half_width = bar_width / 2
        ctrl_offset = (x2 - x1) * 0.4
        verts = [
            (x1 + half_width, y1_bottom),
            (x1 + half_width + ctrl_offset, y1_bottom),
            (x2 - ctrl_offset, y2_bottom),
            (x2, y2_bottom),
            (x2, y2_top),
            (x2 - ctrl_offset, y2_top),
            (x1 + half_width + ctrl_offset, y1_top),
            (x1 + half_width, y1_top),
            (x1 + half_width, y1_bottom),
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

    # Protein Families bar positions (source)
    pf_x = x_left[2]
    pf_terr_bottom = 0
    pf_terr_top = terr_pcts_left[2]
    pf_shared_bottom = pf_terr_top
    pf_shared_top = pf_shared_bottom + shared_pcts_left[2]
    pf_marine_bottom = pf_shared_top
    pf_marine_top = pf_marine_bottom + marine_pcts_left[2]

    # Flow target x position (stop before the bars to give sense of direction)
    flow_target_x = x_right[0] - bar_width / 2 - 0.15

    # Calculate proportional targets for exclusive section
    # Average terrestrial percentage across exclusive entries/species
    avg_terr_pct = sum(exclusive_terr_pcts) / len(exclusive_terr_pcts)  # ~72.7%
    terr_height_in_exclusive = (avg_terr_pct / 50) * exclusive_height  # Scale to 0-45 range

    # Flow 1: Terrestrial exclusive → bottom portion of Exclusive group (0 to ~72.7%)
    draw_flow(
        pf_x,
        pf_terr_bottom,
        pf_terr_top,
        flow_target_x,
        0,
        terr_height_in_exclusive,
        habitat_colors["terrestrial"],
        alpha=0.35,
    )

    # Flow 2: Shared → Shared group (55-100)
    draw_flow(
        pf_x,
        pf_shared_bottom,
        pf_shared_top,
        flow_target_x,
        gap_top,
        100,
        "#7570b3",
        alpha=0.35,
    )

    # Flow 3: Marine exclusive → top portion of Exclusive group (~72.7% to 100%)
    draw_flow(
        pf_x,
        pf_marine_bottom,
        pf_marine_top,
        flow_target_x,
        terr_height_in_exclusive,
        exclusive_height,
        habitat_colors["marine"],
        alpha=0.35,
    )

    # === Axis configuration ===
    all_x = list(x_left) + list(x_right)
    all_labels = metrics_left + ["Entries", "Species"]
    ax.set_xticks(all_x)
    ax.set_xticklabels(all_labels)
    ax.set_ylabel("Percentage (%)")
    ax.set_ylim(0, 105)
    ax.set_xlim(-0.3, 5.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Hide bottom spine (no x-axis line)
    ax.spines["bottom"].set_visible(False)

    # Add two separate y-axes for right-side bars (Exclusive bottom, Shared top)
    right_axis_x = x_right[0] - bar_width / 2 - 0.15
    tick_length = 0.02

    # Bottom y-axis for Exclusive section (0-45 maps to 0-100%)
    ax.plot([right_axis_x, right_axis_x], [0, gap_bottom], color="black", lw=0.8, clip_on=False)
    for tick_val, label in [(0, "0"), (gap_bottom / 2, "50"), (gap_bottom, "100")]:
        ax.plot(
            [right_axis_x - tick_length, right_axis_x],
            [tick_val, tick_val],
            color="black",
            lw=0.8,
            clip_on=False,
        )
        ax.text(
            right_axis_x - tick_length - 0.05,
            tick_val,
            label,
            ha="right",
            va="center",
            fontsize=9,
        )

    # Top y-axis for Shared section (55-100 maps to 0-100%)
    ax.plot([right_axis_x, right_axis_x], [gap_top, 100], color="black", lw=0.8, clip_on=False)
    for tick_val, label in [(gap_top, "0"), ((gap_top + 100) / 2, "50"), (100, "100")]:
        ax.plot(
            [right_axis_x - tick_length, right_axis_x],
            [tick_val, tick_val],
            color="black",
            lw=0.8,
            clip_on=False,
        )
        ax.text(
            right_axis_x - tick_length - 0.05,
            tick_val,
            label,
            ha="right",
            va="center",
            fontsize=9,
        )

    # Add dashed horizontal line between Shared and Exclusive groups on the right
    ax.hlines(
        y=(gap_bottom + gap_top) / 2,
        xmin=x_right[0] - bar_width / 2 - 0.1,
        xmax=x_right[1] + bar_width / 2 + 0.1,
        colors="gray",
        linestyles="dashed",
        linewidth=1.0,
    )

    # Add group labels on the right side (centered in each section)
    shared_center_y = gap_top + shared_height / 2  # Center of shared section
    exclusive_center_y = exclusive_height / 2  # Center of exclusive section
    label_x = x_right[1] + bar_width / 2 + 0.3
    ax.text(
        label_x,
        shared_center_y,
        f"Shared\n(n={n_shared})",
        ha="left",
        va="center",
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        label_x,
        exclusive_center_y,
        f"Exclusive\n(n={n_terr_exclusive + n_marine_exclusive})",
        ha="left",
        va="center",
        fontsize=9,
        fontweight="bold",
    )

    # Move legend to top right, above Shared label
    ax.legend(loc="lower left", bbox_to_anchor=(0.85, 0.82), fontsize=9, framealpha=1.0)


def plot_habitat_combined(
    df: pd.DataFrame,
    datasets: dict[int, pd.DataFrame],
    output_path: Path,
    top_n: int = TOP_N_FAMILIES,
) -> None:
    """Create combined two-panel habitat figure with taxa distribution and protein families.

    Args:
        df: DataFrame for Panel A (taxa by habitat). Should be 2025 data.
        datasets: Dictionary mapping years to DataFrames for Panel B (protein families).
        output_path: Path for combined output figure.
        top_n: Number of top protein families to display in Panel B.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Create figure with independent panel positioning
    fig = plt.figure(figsize=(14, 12))

    # Panel A: positioned with minimal left margin (can extend further left)
    ax_a = fig.add_axes([0.06, 0.55, 0.88, 0.40])  # [left, bottom, width, height]

    # Panel B: positioned with larger left margin to accommodate long y-axis labels
    ax_b = fig.add_axes([0.22, 0.08, 0.72, 0.42])

    # Panel A: Taxa by habitat
    plot_taxa_by_habitat(df, ax_a)

    # Panel B: Dual-habitat protein families
    plot_habitat_protein_families(datasets, ax_b, top_n=top_n)

    # Add panel labels at absolute left edge of figure
    fig.text(0.01, 0.96, "A", fontsize=18, fontweight="bold")
    fig.text(0.01, 0.50, "B", fontsize=18, fontweight="bold")

    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Generated: {output_path}")
    plt.close(fig)
