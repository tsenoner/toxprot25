import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from pathlib import Path
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker

# --- Configuration ---
TOP_N_FAMILIES = 15
SORT_DIVERGING_BARS_BY_SUMMED_INCREASE = True

# --- Path Definitions ---
try:
    BASE_PATH = Path(__file__).resolve().parent.parent
except NameError:  # Fallback for interactive environments
    BASE_PATH = (Path.cwd()) if Path.cwd().name == "toxprot25" else Path.cwd() / ".."

DATA_PATH = BASE_PATH / "data" / "processed"
FIGURE_OUTPUT_PATH = BASE_PATH / "figures" / "habitat"
FIGURE_OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

print(f"Project BASE_PATH resolved to: {BASE_PATH}")
print(f"Data will be loaded from: {DATA_PATH}")
print(f"Figures will be saved to: {FIGURE_OUTPUT_PATH}")

# --- Data Loading ---
try:
    toxprot_2017_df_orig = pd.read_csv(DATA_PATH / "toxprot_2017.csv")
    toxprot_2025_df_orig = pd.read_csv(DATA_PATH / "toxprot_2025.csv")
except FileNotFoundError as e:
    print(
        f"Error loading original data: {e}. Please ensure CSV files are in {DATA_PATH}"
    )
    exit()


# --- Helper Functions ---
def get_protein_family_by_habitat(df, habitat_type):
    """Extract protein family counts for a specific habitat type from a given DataFrame."""
    habitat_data = df[df["Habitat"] == habitat_type].copy()
    habitat_data["Protein families"] = (
        habitat_data["Protein families"].fillna("").astype(str)
    )
    return (
        habitat_data[habitat_data["Protein families"] != ""]["Protein families"]
        .value_counts()
        .sort_index()
    )


def calculate_percentage_increase(old_val, new_val):
    """Safely calculate percentage increase, handling division by zero."""
    if old_val == 0:
        return np.inf if new_val > 0 else 0.0
    return ((new_val - old_val) / old_val) * 100


# --- Step 1: Identify families present in both terrestrial and marine habitats globally ---
print(
    "Identifying protein families present in both terrestrial and marine habitats globally..."
)
pf_terrestrial_2017_orig = get_protein_family_by_habitat(
    toxprot_2017_df_orig, "terrestrial"
)
pf_marine_2017_orig = get_protein_family_by_habitat(toxprot_2017_df_orig, "marine")
pf_terrestrial_2025_orig = get_protein_family_by_habitat(
    toxprot_2025_df_orig, "terrestrial"
)
pf_marine_2025_orig = get_protein_family_by_habitat(toxprot_2025_df_orig, "marine")

all_terrestrial_families_orig = set(pf_terrestrial_2017_orig.index) | set(
    pf_terrestrial_2025_orig.index
)
all_marine_families_orig = set(pf_marine_2017_orig.index) | set(
    pf_marine_2025_orig.index
)
families_in_both_habitats_globally = (
    all_terrestrial_families_orig & all_marine_families_orig
)

if not families_in_both_habitats_globally:
    print(
        "CRITICAL WARNING: No protein families found with representatives in both terrestrial AND marine habitats across both datasets."
    )
elif len(families_in_both_habitats_globally) < TOP_N_FAMILIES:
    print(
        f"INFO: Found {len(families_in_both_habitats_globally)} families with representatives in both habitats globally."
    )
else:
    print(
        f"INFO: Found {len(families_in_both_habitats_globally)} families with representatives in both habitats globally. Will proceed with these."
    )


# --- Step 2: Filter original DataFrames ---
def filter_df_by_families(df, families_to_keep):
    """Filters DataFrame to keep only rows belonging to specified protein families."""
    if not families_to_keep:
        return pd.DataFrame(columns=df.columns)
    temp_families_col = (
        df["Protein families"]
        .fillna(
            "__NAN_PLACEHOLDER__"
        )  # Retained placeholder logic as it might handle specific edge cases.
        .astype(str)
    )
    mask = temp_families_col.isin(families_to_keep)
    return df[mask].copy()


print("Filtering DataFrames to include only dual-habitat protein families...")
toxprot_2017_df = filter_df_by_families(
    toxprot_2017_df_orig, families_in_both_habitats_globally
)
toxprot_2025_df = filter_df_by_families(
    toxprot_2025_df_orig, families_in_both_habitats_globally
)

print(f"Shape of toxprot_2017_df after filtering: {toxprot_2017_df.shape}")
print(f"Shape of toxprot_2025_df after filtering: {toxprot_2025_df.shape}")

if (
    toxprot_2017_df.empty
    and toxprot_2025_df.empty
    and families_in_both_habitats_globally
):
    print("Warning: Filtered DataFrames are empty.")

# --- Step 3: Data Preparation (using filtered DataFrames) ---
print("Preparing data using filtered DataFrames...")
total_terrestrial_2017 = toxprot_2017_df[
    toxprot_2017_df["Habitat"] == "terrestrial"
].shape[0]
total_marine_2017 = toxprot_2017_df[toxprot_2017_df["Habitat"] == "marine"].shape[0]
total_2017_sum_filtered = total_terrestrial_2017 + total_marine_2017

total_terrestrial_2025 = toxprot_2025_df[
    toxprot_2025_df["Habitat"] == "terrestrial"
].shape[0]
total_marine_2025 = toxprot_2025_df[toxprot_2025_df["Habitat"] == "marine"].shape[0]
total_2025_sum_filtered = total_terrestrial_2025 + total_marine_2025

pf_terrestrial_2017 = get_protein_family_by_habitat(toxprot_2017_df, "terrestrial")
pf_marine_2017 = get_protein_family_by_habitat(toxprot_2017_df, "marine")
pf_terrestrial_2025 = get_protein_family_by_habitat(toxprot_2025_df, "terrestrial")
pf_marine_2025 = get_protein_family_by_habitat(toxprot_2025_df, "marine")

s_2017 = pf_terrestrial_2017.add(pf_marine_2017, fill_value=0)
s_2025 = pf_terrestrial_2025.add(pf_marine_2025, fill_value=0)
overall_total_counts_pf_filtered = s_2017.add(s_2025, fill_value=0).sort_values(
    ascending=False
)
top_families = overall_total_counts_pf_filtered.head(TOP_N_FAMILIES).index.tolist()

if not top_families:
    print("Warning: No top families could be selected after filtering.")
elif len(top_families) < TOP_N_FAMILIES and not overall_total_counts_pf_filtered.empty:
    print(
        f"INFO: Selected {len(top_families)} top families out of {len(overall_total_counts_pf_filtered)} dual-habitat families with counts."
    )


# --- Plotting Functions (Matplotlib) ---
def plot_100_stacked_bar_habitat_vs_year_mpl():
    title = "Terrestrial vs. Marine Distribution (Dual-Habitat Families Only)"
    labels = ["2017", "2025"]
    terrestrial_abs_counts = [total_terrestrial_2017, total_terrestrial_2025]
    marine_abs_counts = [total_marine_2017, total_marine_2025]

    totals_for_percentages = np.array(
        [total_2017_sum_filtered, total_2025_sum_filtered]
    )

    safe_totals_for_percentages = np.where(
        totals_for_percentages == 0, 1, totals_for_percentages
    )

    terrestrial_percentages = (
        np.divide(
            terrestrial_abs_counts,
            safe_totals_for_percentages,
            out=np.zeros_like(terrestrial_abs_counts, dtype=float),
            where=totals_for_percentages != 0,
        )
        * 100
    )
    marine_percentages = (
        np.divide(
            marine_abs_counts,
            safe_totals_for_percentages,
            out=np.zeros_like(marine_abs_counts, dtype=float),
            where=totals_for_percentages != 0,
        )
        * 100
    )

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.bar(labels, terrestrial_percentages, label="Terrestrial", color="forestgreen")
    ax.bar(
        labels,
        marine_percentages,
        bottom=terrestrial_percentages,
        label="Marine",
        color="steelblue",
    )

    ax.set_ylabel("Percentage of Sequences (%)")
    ax.set_xlabel("Year")
    ax.set_title(title, pad=20)
    ax.legend()
    ax.set_ylim(0, 100)

    for i, (t_perc, m_perc, t_abs, m_abs) in enumerate(
        zip(
            terrestrial_percentages,
            marine_percentages,
            terrestrial_abs_counts,
            marine_abs_counts,
        )
    ):
        if t_perc > 0 or t_abs > 0:
            ax.text(
                i,
                t_perc / 2,
                f"{t_perc:.1f}%\n({t_abs})",
                ha="center",
                va="center",
                color="white",
                fontweight="bold",
            )
        if m_perc > 0 or m_abs > 0:
            ax.text(
                i,
                t_perc + (m_perc / 2),
                f"{m_perc:.1f}%\n({m_abs})",
                ha="center",
                va="center",
                color="white",
                fontweight="bold",
            )

    plt.tight_layout()
    filepath = FIGURE_OUTPUT_PATH / "1_stacked_bar_habitat_vs_year_mpl.png"
    plt.savefig(filepath)
    plt.close(fig)


def plot_percentage_increase_by_habitat_mpl():
    title = "Percentage Increase (2017-2025) by Habitat (Dual-Habitat Families Only)"
    habitats = ["Terrestrial", "Marine"]

    abs_terrestrial_2017 = total_terrestrial_2017
    abs_terrestrial_2025 = total_terrestrial_2025
    abs_marine_2017 = total_marine_2017
    abs_marine_2025 = total_marine_2025

    increase_terrestrial = calculate_percentage_increase(
        abs_terrestrial_2017, abs_terrestrial_2025
    )
    increase_marine = calculate_percentage_increase(abs_marine_2017, abs_marine_2025)

    increases = [increase_terrestrial, increase_marine]
    abs_counts_list = [
        (abs_terrestrial_2017, abs_terrestrial_2025),
        (abs_marine_2017, abs_marine_2025),
    ]
    colors = ["forestgreen", "steelblue"]

    fig, ax = plt.subplots(figsize=(10, 7))
    bars = ax.bar(habitats, increases, color=colors)

    ax.set_ylabel("Percentage Increase (%)")
    ax.set_title(title, pad=20)

    for bar, inc, counts in zip(bars, increases, abs_counts_list):
        perc_text = (
            f"{inc:.1f}%" if np.isfinite(inc) else "Inf" if inc == np.inf else "0%"
        )
        abs_text = f"({counts[0]} → {counts[1]})"
        full_text = f"{perc_text}\n{abs_text}"

        bar_height = bar.get_height()
        if bar_height > 0:
            text_y = bar_height * 0.5
            va_align = "center"
        elif bar_height < 0:
            text_y = bar_height * 0.5
            va_align = "center"
        else:
            text_y = 0.1
            va_align = "bottom"

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            text_y,
            full_text,
            ha="center",
            va=va_align,
            fontsize=12,
            linespacing=1.3,
            color="white",
            fontweight="bold",
        )

    plt.xticks(rotation=0, ha="center")
    plt.tight_layout()
    filepath = FIGURE_OUTPUT_PATH / "2_percentage_increase_by_habitat_mpl.png"
    plt.savefig(filepath)
    plt.close(fig)


def plot_heatmap_percentage_increase_mpl():
    title = f"Heatmap: % Increase by Habitat & Top {TOP_N_FAMILIES} Protein Families (Dual-Habitat Only)"

    if not top_families:
        print("Skipping heatmap: No top families identified.")
        return

    heatmap_data = pd.DataFrame(index=top_families, columns=["Terrestrial", "Marine"])

    for family in top_families:
        t_old = pf_terrestrial_2017.get(family, 0)
        t_new = pf_terrestrial_2025.get(family, 0)
        m_old = pf_marine_2017.get(family, 0)
        m_new = pf_marine_2025.get(family, 0)

        heatmap_data.loc[family, "Terrestrial"] = calculate_percentage_increase(
            t_old, t_new
        )
        heatmap_data.loc[family, "Marine"] = calculate_percentage_increase(m_old, m_new)

    heatmap_data_numeric = heatmap_data.apply(pd.to_numeric, errors="coerce").fillna(0)
    heatmap_plot_data = heatmap_data_numeric.copy()

    data_for_imshow = heatmap_plot_data.replace([np.inf, -np.inf], [100.0, -100.0])

    num_rows, num_cols = heatmap_plot_data.shape

    bounds_for_norm = [0, 20, 40, 60, 80, np.inf]
    colors = ["lightblue", "cornflowerblue", "royalblue", "mediumblue", "darkblue"]
    cmap_categorical = mcolors.ListedColormap(colors)
    norm_for_imshow = mcolors.BoundaryNorm(bounds_for_norm, cmap_categorical.N)

    cell_size = 0.9
    max_family_name_len = (
        max(len(name) for name in heatmap_plot_data.index)
        if not heatmap_plot_data.empty
        else 10
    )
    estimated_y_label_width = max_family_name_len * 0.15
    fig_width = num_cols * cell_size + estimated_y_label_width + 5.0
    fig_height = num_rows * cell_size + 2.5

    fig, ax = plt.subplots(figsize=(max(8, fig_width), max(6, fig_height)))

    ax.imshow(
        data_for_imshow,
        cmap=cmap_categorical,
        norm=norm_for_imshow,
        aspect="equal",
    )

    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes("right", size="20%", pad=0.8)

    cbar_display_bounds = [0, 20, 40, 60, 80, 100]
    norm_for_cbar = mcolors.BoundaryNorm(cbar_display_bounds, cmap_categorical.N)
    sm = cm.ScalarMappable(cmap=cmap_categorical, norm=norm_for_cbar)

    cb = fig.colorbar(sm, cax=cbar_ax, ticks=[10, 30, 50, 70, 90])
    cb.set_label("% Increase Categories", fontsize=20)
    cb.set_ticklabels(["0-20%", "20-40%", "40-60%", "60-80%", "80%-Inf"])
    cb.ax.tick_params(labelsize=18)

    ax.set_xticks(np.arange(num_cols))
    ax.set_yticks(np.arange(num_rows))
    ax.set_xticklabels(heatmap_plot_data.columns)
    ax.set_yticklabels(heatmap_plot_data.index)
    ax.set_title(title, pad=20, fontsize=20, loc="center")

    plt.setp(ax.get_xticklabels(), rotation=0, ha="center", rotation_mode="anchor")
    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=18)

    for i in range(num_rows):
        for j in range(num_cols):
            family = heatmap_data.index[i]
            habitat = heatmap_data.columns[j]
            perc_val = heatmap_data.iloc[i, j]

            if habitat == "Terrestrial":
                old_count = pf_terrestrial_2017.get(family, 0)
                new_count = pf_terrestrial_2025.get(family, 0)
            else:
                old_count = pf_marine_2017.get(family, 0)
                new_count = pf_marine_2025.get(family, 0)

            text_perc = (
                f"{perc_val:.0f}%"
                if np.isfinite(perc_val)
                else "Inf"
                if perc_val > 0
                else "-Inf"
                if perc_val < 0
                else "0%"
            )
            text_abs = f"({old_count}→{new_count})"
            full_text = f"{text_perc}\n{text_abs}"

            if perc_val == np.inf or perc_val >= 80:
                text_color = "white"
            else:
                text_color = "black"

            ax.text(
                j,
                i,
                full_text,
                ha="center",
                va="center",
                color=text_color,
                fontsize=12,
                linespacing=1.3,
            )

    # Apply a general tight layout first.
    fig.tight_layout()
    # Now, adjust subplot parameters to shift content.
    fig.subplots_adjust(
        left=0.25,  # Start plots very close to left edge
        right=0.9,  # Allow plots to extend, leaving whitespace on far right
        top=0.95,  # Space for title
        bottom=0.05,
    )  # Space for x-axis labels if any are very long

    filepath = FIGURE_OUTPUT_PATH / "5_heatmap_percentage_increase_mpl.png"
    plt.savefig(filepath)
    print(f"Generated: {filepath}")
    plt.close(fig)


def plot_diverging_bar_charts_protein_family_mpl():
    if not top_families:
        print("Skipping diverging bar charts: No top families identified.")
        return

    protein_family_summary = pd.DataFrame(index=top_families)
    protein_family_summary["Terrestrial_2017"] = pf_terrestrial_2017.reindex(
        top_families, fill_value=0
    )
    protein_family_summary["Marine_2017"] = pf_marine_2017.reindex(
        top_families, fill_value=0
    )
    protein_family_summary["Terrestrial_2025"] = pf_terrestrial_2025.reindex(
        top_families, fill_value=0
    )
    protein_family_summary["Marine_2025"] = pf_marine_2025.reindex(
        top_families, fill_value=0
    )

    protein_family_summary["Terrestrial_Change_Abs"] = (
        protein_family_summary["Terrestrial_2025"]
        - protein_family_summary["Terrestrial_2017"]
    )
    protein_family_summary["Marine_Change_Abs"] = (
        protein_family_summary["Marine_2025"] - protein_family_summary["Marine_2017"]
    )

    protein_family_summary["Terrestrial_Change_Perc"] = [
        calculate_percentage_increase(old, new)
        for old, new in zip(
            protein_family_summary["Terrestrial_2017"],
            protein_family_summary["Terrestrial_2025"],
        )
    ]
    protein_family_summary["Marine_Change_Perc"] = [
        calculate_percentage_increase(old, new)
        for old, new in zip(
            protein_family_summary["Marine_2017"], protein_family_summary["Marine_2025"]
        )
    ]

    if SORT_DIVERGING_BARS_BY_SUMMED_INCREASE:
        protein_family_summary["Summed_Abs_Increase"] = (
            protein_family_summary["Terrestrial_Change_Abs"]
            + protein_family_summary["Marine_Change_Abs"]
        )
        protein_family_summary = protein_family_summary.sort_values(
            by="Summed_Abs_Increase", ascending=False
        )

    title_abs = f"Top {len(protein_family_summary)} Protein Families: Absolute Change (2017-2025, Dual-Habitat)"
    y_labels_abs = protein_family_summary.index
    terrestrial_values_abs = protein_family_summary["Terrestrial_Change_Abs"]
    marine_values_abs = protein_family_summary["Marine_Change_Abs"]

    fig_abs, ax_abs = plt.subplots(figsize=(12, max(6, len(y_labels_abs) * 0.4)))
    y_pos_abs = np.arange(len(y_labels_abs))

    ax_abs.barh(
        y_pos_abs,
        terrestrial_values_abs,
        color="forestgreen",
        label="Terrestrial",
    )
    ax_abs.barh(
        y_pos_abs,
        -marine_values_abs,
        color="steelblue",
        label="Marine",
    )

    ax_abs.set_yticks(y_pos_abs)
    ax_abs.set_yticklabels(y_labels_abs)
    ax_abs.invert_yaxis()
    ax_abs.set_xlabel("Absolute Change (Marine Magnitude | Terrestrial Magnitude)")
    ax_abs.set_title(title_abs, pad=20)
    ax_abs.legend(loc="lower right")
    ax_abs.axvline(0, color="grey", lw=0.8)

    # Custom formatter to show absolute values on x-axis
    def abs_formatter(x, pos):
        return f"{abs(x):.0f}"

    ax_abs.xaxis.set_major_formatter(mticker.FuncFormatter(abs_formatter))

    for i, (t_val, m_val) in enumerate(zip(terrestrial_values_abs, marine_values_abs)):
        # Always draw text for t_val, even if 0
        ax_abs.text(
            t_val
            + (
                0.01 * ax_abs.get_xlim()[1]
                if t_val >= 0  # Position to the right for 0 and positive
                else -0.03 * ax_abs.get_xlim()[1]  # Further left for negative
            ),
            y_pos_abs[i],
            f"{t_val}",
            va="center",
            ha="left" if t_val >= 0 else "right",  # Adjust ha for 0
            fontsize=8,
        )
        # Always draw text for m_val, even if 0
        ax_abs.text(
            -m_val  # Plotting position for marine bar
            - (
                0.01 * ax_abs.get_xlim()[1]
                if m_val
                >= 0  # Position to the left for 0 and positive marine values (plotted on negative side)
                else -0.03
                * ax_abs.get_xlim()[1]  # Further right for negative marine values
            ),
            y_pos_abs[i],
            f"{m_val}",
            va="center",
            ha="right" if m_val >= 0 else "left",  # Adjust ha for 0 on marine side
            fontsize=8,
        )

    current_xlim_abs = ax_abs.get_xlim()
    expansion_factor_abs = 0.25
    ax_abs.set_xlim(
        current_xlim_abs[0] - abs(current_xlim_abs[0]) * expansion_factor_abs,
        current_xlim_abs[1] + abs(current_xlim_abs[1]) * expansion_factor_abs,
    )

    plt.tight_layout()
    filepath_abs = (
        FIGURE_OUTPUT_PATH / "6a_diverging_bar_abs_change_protein_family_mpl.png"
    )
    plt.savefig(filepath_abs)
    print(f"Generated: {filepath_abs}")
    plt.close(fig_abs)

    title_perc = f"Top {len(protein_family_summary)} Protein Families: Percentage Change (2017-2025, Dual-Habitat)"
    y_labels_perc = protein_family_summary.index

    terrestrial_values_perc = protein_family_summary["Terrestrial_Change_Perc"].replace(
        [np.inf, -np.inf],
        100,
    )
    marine_values_perc = protein_family_summary["Marine_Change_Perc"].replace(
        [np.inf, -np.inf],
        100,
    )

    terrestrial_original_perc = protein_family_summary["Terrestrial_Change_Perc"]
    marine_original_perc = protein_family_summary["Marine_Change_Perc"]

    fig_perc, ax_perc = plt.subplots(figsize=(12, max(6, len(y_labels_perc) * 0.4)))
    y_pos_perc = np.arange(len(y_labels_perc))

    ax_perc.barh(
        y_pos_perc,
        terrestrial_values_perc,
        color="forestgreen",
        label="Terrestrial",
    )
    ax_perc.barh(
        y_pos_perc,
        -marine_values_perc,
        color="steelblue",
        label="Marine",
    )

    ax_perc.set_yticks(y_pos_perc)
    ax_perc.set_yticklabels(y_labels_perc)
    ax_perc.invert_yaxis()
    ax_perc.set_xlabel(
        "% Change (Marine Magnitude | Terrestrial Magnitude) [Plot capped at +/-100%]"
    )
    ax_perc.set_title(title_perc, pad=20)
    ax_perc.legend(loc="lower right")
    ax_perc.axvline(0, color="grey", lw=0.8)

    # Apply custom formatter to percentage plot as well
    # For percentages, we might want one decimal place if numbers are small
    def abs_perc_formatter(x, pos):
        return f"{abs(x):.0f}"  # Keeping .0f for consistency with absolute, adjust if .1f desired for perc

    ax_perc.xaxis.set_major_formatter(mticker.FuncFormatter(abs_perc_formatter))

    for i, (t_orig_perc, m_orig_perc, t_plot_perc, m_plot_perc) in enumerate(
        zip(
            terrestrial_original_perc,
            marine_original_perc,
            terrestrial_values_perc,
            marine_values_perc,
        )
    ):
        if pd.notna(
            t_orig_perc
        ):  # Draw if t_orig_perc is a valid number (including 0.0)
            t_text = (
                "Inf"
                if np.isinf(t_orig_perc)
                else f"{t_orig_perc:.1f}%"
                if pd.notna(
                    t_orig_perc
                )  # This condition is somewhat redundant due to outer if, but safe
                else "N/A"  # Should not be reached if outer if pd.notna(t_orig_perc) is true
            )
            ax_perc.text(
                t_plot_perc  # Use plot_perc for position (handles capping for Inf)
                + (
                    0.01 * ax_perc.get_xlim()[1]
                    if t_plot_perc >= 0  # Position to the right for 0 and positive
                    else -0.03 * ax_perc.get_xlim()[1]  # Further left for negative
                ),
                y_pos_perc[i],
                t_text,
                va="center",
                ha="left" if t_plot_perc >= 0 else "right",  # Adjust ha for 0
                fontsize=8,
            )

        if pd.notna(
            m_orig_perc
        ):  # Draw if m_orig_perc is a valid number (including 0.0)
            m_text = (
                "Inf"
                if np.isinf(m_orig_perc)
                else f"{m_orig_perc:.1f}%"
                if pd.notna(m_orig_perc)
                else "N/A"
            )
            ax_perc.text(
                -m_plot_perc  # Use plot_perc for position
                - (
                    0.01 * ax_perc.get_xlim()[1]
                    if m_plot_perc
                    >= 0  # Position to the left for 0 and positive marine values
                    else -0.03
                    * ax_perc.get_xlim()[1]  # Further right for negative marine values
                ),
                y_pos_perc[i],
                m_text,
                va="center",
                ha="right"
                if m_plot_perc >= 0
                else "left",  # Adjust ha for 0 on marine side
                fontsize=8,
            )

    current_xlim_perc = ax_perc.get_xlim()
    expansion_factor_perc = 0.25
    ax_perc.set_xlim(
        current_xlim_perc[0] - abs(current_xlim_perc[0]) * expansion_factor_perc,
        current_xlim_perc[1] + abs(current_xlim_perc[1]) * expansion_factor_perc,
    )

    plt.tight_layout()
    filepath_perc = (
        FIGURE_OUTPUT_PATH / "6b_diverging_bar_perc_change_protein_family_mpl.png"
    )
    plt.savefig(filepath_perc)
    print(f"Generated: {filepath_perc}")
    plt.close(fig_perc)

    # --- Plot 6c: Current Absolute Numbers (2025) with 2017 context ---
    title_c = f"Top {len(protein_family_summary)} Protein Families: Absolute Numbers 2025 (Dual-Habitat)"
    y_labels_c = protein_family_summary.index  # Same sorted order

    terrestrial_values_2017 = protein_family_summary["Terrestrial_2017"]
    terrestrial_values_2025 = protein_family_summary["Terrestrial_2025"]
    marine_values_2017 = protein_family_summary["Marine_2017"]
    marine_values_2025 = protein_family_summary["Marine_2025"]

    fig_c, ax_c = plt.subplots(figsize=(12, max(6, len(y_labels_c) * 0.4)))
    y_pos_c = np.arange(len(y_labels_c))

    # Define bright colors for added portions
    bright_green = "#31c42f"  # A bright green
    bright_blue = "#56a0dd"  # A bright blue

    # Terrestrial bars (right side)
    # Portion of 2025 bar that was effectively present in 2017 (or full 2025 bar if count decreased)
    t_existing_or_total_if_decrease = np.minimum(
        terrestrial_values_2017, terrestrial_values_2025
    )
    # Portion of 2025 bar that is new/added since 2017 (will be 0 if count decreased)
    t_added_since_2017 = np.maximum(
        0, terrestrial_values_2025 - terrestrial_values_2017
    )

    ax_c.barh(
        y_pos_c,
        t_existing_or_total_if_decrease,
        color="forestgreen",
        label="Terrestrial (2025)",  # Main label for terrestrial
    )
    ax_c.barh(
        y_pos_c,
        t_added_since_2017,
        left=t_existing_or_total_if_decrease,  # Stack on top of the existing part
        color=bright_green,
        # No label for this part to keep legend clean
    )

    # Marine bars (left side, plotted with negative values)
    # Portion of 2025 bar that was effectively present in 2017 (or full 2025 bar if count decreased) - positive value
    m_existing_or_total_if_decrease_positive = np.minimum(
        marine_values_2017, marine_values_2025
    )
    # Portion of 2025 bar that is new/added since 2017 - positive value
    m_added_since_2017_positive = np.maximum(0, marine_values_2025 - marine_values_2017)

    ax_c.barh(
        y_pos_c,
        -m_existing_or_total_if_decrease_positive,  # Plot as negative
        color="steelblue",
        label="Marine (2025)",  # Main label for marine
    )
    ax_c.barh(
        y_pos_c,
        -m_added_since_2017_positive,  # Plot as negative
        left=-m_existing_or_total_if_decrease_positive,  # Stack on the 'left' end of the existing part
        color=bright_blue,
        # No label for this part
    )

    ax_c.set_yticks(y_pos_c)
    ax_c.set_yticklabels(y_labels_c)
    ax_c.invert_yaxis()  # Highest summed increase (from original sort) at the top
    ax_c.set_xlabel(
        "Absolute Number of Sequences (Marine Magnitude | Terrestrial Magnitude) - 2025 Data"
    )
    ax_c.set_title(title_c, pad=20)
    ax_c.legend(loc="lower right")
    ax_c.axvline(0, color="grey", lw=0.8)

    # Apply custom formatter for absolute values on x-axis (reusing abs_formatter)
    ax_c.xaxis.set_major_formatter(
        mticker.FuncFormatter(abs_formatter)
    )  # abs_formatter defined in Plot 6a section

    # Add text annotations: (2017 count -> 2025 count)
    for i, family_name in enumerate(y_labels_c):
        t_2017 = terrestrial_values_2017[family_name]
        t_2025 = terrestrial_values_2025[family_name]
        m_2017 = marine_values_2017[family_name]
        m_2025 = marine_values_2025[family_name]

        # Terrestrial text (right)
        ax_c.text(
            t_2025
            + (
                0.01 * ax_c.get_xlim()[1] if t_2025 >= 0 else -0.03 * ax_c.get_xlim()[1]
            ),
            y_pos_c[i],
            f"{t_2017}→{t_2025}",
            va="center",
            ha="left" if t_2025 >= 0 else "right",
            fontsize=8,
            # Ensure no specific color is set here, to use default
        )
        # Marine text (left)
        ax_c.text(
            -m_2025
            - (
                0.01 * ax_c.get_xlim()[1] if m_2025 >= 0 else -0.03 * ax_c.get_xlim()[1]
            ),
            y_pos_c[i],
            f"{m_2017}→{m_2025}",
            va="center",
            ha="right" if m_2025 >= 0 else "left",
            fontsize=8,
            # Ensure no specific color is set here, to use default
        )

    # Adjust x-limits to ensure text is visible
    current_xlim_c = ax_c.get_xlim()
    expansion_factor_c = 0.25
    ax_c.set_xlim(
        current_xlim_c[0] - abs(current_xlim_c[0]) * expansion_factor_c,
        current_xlim_c[1] + abs(current_xlim_c[1]) * expansion_factor_c,
    )

    plt.tight_layout()
    filepath_c = (
        FIGURE_OUTPUT_PATH / "6c_diverging_bar_abs_numbers_protein_family_mpl.png"
    )
    plt.savefig(filepath_c)
    print(f"Generated: {filepath_c}")
    plt.close(fig_c)


def generate_dual_habitat_protein_counts_summary():
    """
    Generates a CSV file summarizing protein counts per habitat, year, and protein family,
    considering only families present in both terrestrial and marine habitats (dual-habitat).
    """
    print("Generating summary file of dual-habitat protein counts in wide format...")

    s_t_2017 = pf_terrestrial_2017.rename("2017_Terrestrial")
    s_m_2017 = pf_marine_2017.rename("2017_Marine")
    s_t_2025 = pf_terrestrial_2025.rename("2025_Terrestrial")
    s_m_2025 = pf_marine_2025.rename("2025_Marine")

    summary_df = pd.concat([s_t_2017, s_m_2017, s_t_2025, s_m_2025], axis=1)
    summary_df = summary_df.fillna(0).astype(int)
    summary_df = summary_df.reset_index().rename(
        columns={"Protein families": "ProteinFamily"}
    )

    if summary_df.empty:
        print(
            "No data to summarize for dual-habitat protein counts file. Skipping file generation."
        )
        return

    summary_df = summary_df.sort_values(by="ProteinFamily").reset_index(drop=True)

    output_filename = "dual_habitat_protein_counts_summary.csv"
    output_filepath = FIGURE_OUTPUT_PATH / output_filename
    summary_df.to_csv(output_filepath, index=False)
    print(f"Generated summary file: {output_filepath}")


# --- Main Execution ---
if __name__ == "__main__":
    if not families_in_both_habitats_globally:
        print("Halting script: No dual-habitat protein families identified.")
    elif (
        toxprot_2017_df.empty
        and toxprot_2025_df.empty
        and families_in_both_habitats_globally
    ):
        print(
            "Halting script: Filtered DataFrames for dual-habitat families are empty. Cannot generate meaningful plots."
        )
    else:
        generate_dual_habitat_protein_counts_summary()
        plot_100_stacked_bar_habitat_vs_year_mpl()
        plot_percentage_increase_by_habitat_mpl()
        plot_heatmap_percentage_increase_mpl()
        plot_diverging_bar_charts_protein_family_mpl()
        print(
            f"All Matplotlib plots and summary file attempted to be saved to {FIGURE_OUTPUT_PATH}"
        )
