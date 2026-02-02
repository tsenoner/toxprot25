"""Analyze protein family distributions by habitat (terrestrial vs marine).

Generates visualizations comparing dual-habitat protein families across
ToxProt datasets from different years (2005, 2015, 2025).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from .analyze_protein_families import FAMILY_NAME_MAP

# --- Configuration ---
TOP_N_FAMILIES = 15
YEARS = [2005, 2015, 2025]


def normalize_family_name(name: str) -> str:
    """Normalize a family name using the mapping dictionary."""
    if pd.isna(name) or name == "":
        return name
    return FAMILY_NAME_MAP.get(name, name)


def load_datasets(
    years: list[int],
    data_dir: Path = Path("data/processed/toxprot"),
) -> dict[int, pd.DataFrame]:
    """Load ToxProt datasets for specified years."""
    datasets = {}
    for year in years:
        filepath = data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)
    return datasets


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
    output_path: Path,
    top_n: int = TOP_N_FAMILIES,
) -> None:
    """Plot dual-habitat protein family evolution across years.

    Args:
        datasets: Dictionary mapping years to DataFrames
        output_path: Path for output figure (without extension)
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

    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(6, len(summary) * 0.45)))
    y_pos = np.arange(len(summary))

    # Colors for the three periods (darker to lighter)
    colors = {
        "terrestrial": ["forestgreen", "#31c42f", "#90EE90"],
        "marine": ["steelblue", "#56a0dd", "#ADD8E6"],
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

    # Custom legend showing years with colors for both habitats (Marine first, then Terrestrial)
    legend_elements = [
        Patch(facecolor=colors["marine"][0], edgecolor="none", label=f"Marine {year_early}"),
        Patch(facecolor=colors["marine"][1], edgecolor="none", label=f"Marine {year_mid}"),
        Patch(facecolor=colors["marine"][2], edgecolor="none", label=f"Marine {year_late}"),
        Patch(facecolor=colors["terrestrial"][0], edgecolor="none", label=f"Terrestrial {year_early}"),
        Patch(facecolor=colors["terrestrial"][1], edgecolor="none", label=f"Terrestrial {year_mid}"),
        Patch(facecolor=colors["terrestrial"][2], edgecolor="none", label=f"Terrestrial {year_late}"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=10, ncol=2)
    ax.axvline(0, color="grey", lw=0.8)
    ax.xaxis.set_major_formatter(mticker.FuncFormatter(abs_formatter))

    # Expand x-limits for annotations (balance space for annotations on both sides)
    xlim = ax.get_xlim()
    ax.set_xlim(xlim[0] * 2.5, xlim[1] * 1.1)
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

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Generated: {output_path}")
    plt.close(fig)
