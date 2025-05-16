import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from pathlib import Path

# --- Configuration ---
TOP_N_FAMILIES = 15
SORT_DIVERGING_BARS_BY_SUMMED_INCREASE = True  # New flag for sorting diverging bars

# --- Path Definitions ---
try:
    BASE_PATH = Path(__file__).resolve().parent.parent
except NameError:  # Fallback for interactive environments
    BASE_PATH = (
        (
            Path.cwd()  # Assuming script is run from project root or notebook dir directly if not __file__
        )
        if Path.cwd().name == "toxprot25"
        else Path.cwd() / ".."
    )

DATA_PATH = BASE_PATH / "data" / "processed"
FIGURE_OUTPUT_PATH = BASE_PATH / "figures" / "taxa" / "habitat_analysis_exploration"
FIGURE_OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

print(f"Project BASE_PATH resolved to: {BASE_PATH}")
print(f"Data will be loaded from: {DATA_PATH}")
print(f"Figures will be saved to: {FIGURE_OUTPUT_PATH}")

# --- Data Loading (Original) ---
try:
    toxprot_2017_df_orig = pd.read_csv(DATA_PATH / "toxprot_2017.csv")
    toxprot_2025_df_orig = pd.read_csv(DATA_PATH / "toxprot_2025.csv")
except FileNotFoundError as e:
    print(
        f"Error loading original data: {e}. Please ensure CSV files are in {DATA_PATH}"
    )
    exit()


# --- Helper Functions ---
def get_protein_family_by_habitat(df, habitat_type, split_char):
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


def get_all_protein_families_global(df, protein_family_col_name, split_char):
    """Extract global protein family counts from a given DataFrame."""
    if protein_family_col_name not in df.columns:
        print(
            f"Warning: Column '{protein_family_col_name}' not found in DataFrame for global family counting. Returning empty Series."
        )
        return pd.Series(dtype="int")

    temp_df = df.copy()
    temp_df[protein_family_col_name] = (
        temp_df[protein_family_col_name].fillna("").astype(str)
    )
    valid_families = temp_df[protein_family_col_name][
        temp_df[protein_family_col_name] != ""
    ]
    family_counts = valid_families.value_counts()
    return family_counts.sort_index()


# --- Step 1: Identify families present in both terrestrial and marine habitats globally ---
print(
    "Identifying protein families present in both terrestrial and marine habitats globally..."
)
pf_terrestrial_2017_orig = get_protein_family_by_habitat(
    toxprot_2017_df_orig, "terrestrial", None
)
pf_marine_2017_orig = get_protein_family_by_habitat(
    toxprot_2017_df_orig, "marine", None
)
pf_terrestrial_2025_orig = get_protein_family_by_habitat(
    toxprot_2025_df_orig, "terrestrial", None
)
pf_marine_2025_orig = get_protein_family_by_habitat(
    toxprot_2025_df_orig, "marine", None
)

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
    # exit()
elif len(families_in_both_habitats_globally) < TOP_N_FAMILIES:
    print(
        f"INFO: Found {len(families_in_both_habitats_globally)} families with representatives in both habitats globally."
    )
else:
    print(
        f"INFO: Found {len(families_in_both_habitats_globally)} families with representatives in both habitats globally. Will proceed with these."
    )


# --- Step 2: Filter original DataFrames ---
def filter_df_by_families(df, families_to_keep, split_char):
    """Filters DataFrame to keep only rows belonging to specified protein families."""
    if not families_to_keep:
        return pd.DataFrame(columns=df.columns)
    temp_families_col = df["Protein families"].fillna("__NAN_PLACEHOLDER__").astype(str)
    mask = temp_families_col.isin(families_to_keep)
    return df[mask].copy()


print("Filtering DataFrames to include only dual-habitat protein families...")
toxprot_2017_df = filter_df_by_families(
    toxprot_2017_df_orig, families_in_both_habitats_globally, None
)
toxprot_2025_df = filter_df_by_families(
    toxprot_2025_df_orig, families_in_both_habitats_globally, None
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
total_2017_sum_filtered = (
    total_terrestrial_2017 + total_marine_2017
)  # Renamed to avoid confusion

total_terrestrial_2025 = toxprot_2025_df[
    toxprot_2025_df["Habitat"] == "terrestrial"
].shape[0]
total_marine_2025 = toxprot_2025_df[toxprot_2025_df["Habitat"] == "marine"].shape[0]
total_2025_sum_filtered = total_terrestrial_2025 + total_marine_2025  # Renamed

pf_terrestrial_2017 = get_protein_family_by_habitat(
    toxprot_2017_df, "terrestrial", None
)
pf_marine_2017 = get_protein_family_by_habitat(toxprot_2017_df, "marine", None)
pf_terrestrial_2025 = get_protein_family_by_habitat(
    toxprot_2025_df, "terrestrial", None
)
pf_marine_2025 = get_protein_family_by_habitat(toxprot_2025_df, "marine", None)

s_2017 = pf_terrestrial_2017.add(pf_marine_2017, fill_value=0)
s_2025 = pf_terrestrial_2025.add(pf_marine_2025, fill_value=0)
overall_total_counts_pf_filtered = s_2017.add(s_2025, fill_value=0).sort_values(
    ascending=False
)
top_families = overall_total_counts_pf_filtered.head(TOP_N_FAMILIES).index.tolist()

if not top_families:
    print(f"Warning: No top families could be selected after filtering.")
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

    # Use renamed total_sum_filtered for clarity
    totals_for_percentages = np.array(
        [total_2017_sum_filtered, total_2025_sum_filtered]
    )

    # Ensure totals_for_percentages are not zero to avoid runtime warnings with np.divide
    safe_totals_for_percentages = np.where(
        totals_for_percentages == 0, 1, totals_for_percentages
    )

    terrestrial_percentages = (
        np.divide(
            terrestrial_abs_counts,
            safe_totals_for_percentages,  # Use safe version
            out=np.zeros_like(terrestrial_abs_counts, dtype=float),
            where=totals_for_percentages
            != 0,  # Original condition still applies for output
        )
        * 100
    )
    marine_percentages = (
        np.divide(
            marine_abs_counts,
            safe_totals_for_percentages,  # Use safe version
            out=np.zeros_like(marine_abs_counts, dtype=float),
            where=totals_for_percentages
            != 0,  # Original condition still applies for output
        )
        * 100
    )

    fig, ax = plt.subplots(figsize=(8, 7))
    bars1 = ax.bar(
        labels, terrestrial_percentages, label="Terrestrial", color="forestgreen"
    )
    bars2 = ax.bar(
        labels,
        marine_percentages,
        bottom=terrestrial_percentages,
        label="Marine",
        color="steelblue",
    )

    ax.set_ylabel("Percentage of Sequences (%)")
    ax.set_xlabel("Year")
    ax.set_title(title, pad=20)  # Added padding to title
    ax.legend()
    ax.set_ylim(0, 100)

    # Add text annotations (Percentage and Absolute Count)
    for i, (t_perc, m_perc, t_abs, m_abs) in enumerate(
        zip(
            terrestrial_percentages,
            marine_percentages,
            terrestrial_abs_counts,
            marine_abs_counts,
        )
    ):
        # Terrestrial annotation
        if t_perc > 0 or t_abs > 0:  # Add text if there's something to show
            ax.text(
                i,
                t_perc / 2,
                f"{t_perc:.1f}%\n({t_abs})",
                ha="center",
                va="center",
                color="white",
                fontweight="bold",
            )
        # Marine annotation
        if m_perc > 0 or m_abs > 0:  # Add text if there's something to show
            ax.text(
                i,
                t_perc + (m_perc / 2),
                f"{m_perc:.1f}%\n({m_abs})",
                ha="center",
                va="center",
                color="white",
                fontweight="bold",
            )

    plt.tight_layout()  # Adjust layout
    filepath = FIGURE_OUTPUT_PATH / "1_stacked_bar_habitat_vs_year_mpl.png"
    plt.savefig(filepath)


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
        abs_text = f"({counts[0]} → {counts[1]})"  # Arrow: →
        full_text = f"{perc_text}\n{abs_text}"

        bar_height = bar.get_height()
        # Determine text_y position to be inside the bar
        if bar_height > 0:
            text_y = bar_height * 0.5  # Middle of positive bar
            va_align = "center"
        elif bar_height < 0:
            text_y = bar_height * 0.5  # Middle of negative bar (plots downwards)
            va_align = "center"
        else:  # Zero height bar
            text_y = 0.1  # Slightly above the baseline if bar height is zero
            va_align = "bottom"

        ax.text(
            bar.get_x() + bar.get_width() / 2,
            text_y,  # Use new y position
            full_text,
            ha="center",
            va=va_align,  # Use dynamic va
            fontsize=11,  # Increased fontsize
            linespacing=1.3,
            color="white",  # Assuming white is generally readable on these bar colors
        )

    plt.xticks(rotation=90, ha="center")  # Rotated x-axis ticks to vertical
    plt.tight_layout()
    filepath = FIGURE_OUTPUT_PATH / "2_percentage_increase_by_habitat_mpl.png"
    plt.savefig(filepath)


def plot_percentage_increase_by_protein_family_mpl():
    title = f"Top {TOP_N_FAMILIES} Protein Families: % Incr. & Abs. Change (Dual-Habitat)"  # Shortened title

    _pf_t_2017 = pf_terrestrial_2017.reindex(top_families, fill_value=0)
    _pf_m_2017 = pf_marine_2017.reindex(top_families, fill_value=0)
    _pf_t_2025 = pf_terrestrial_2025.reindex(top_families, fill_value=0)
    _pf_m_2025 = pf_marine_2025.reindex(top_families, fill_value=0)

    total_pf_2017 = _pf_t_2017.add(_pf_m_2017, fill_value=0)
    total_pf_2025 = _pf_t_2025.add(_pf_m_2025, fill_value=0)

    percentage_increase_pf = pd.Series(
        [
            calculate_percentage_increase(old, new)
            for old, new in zip(total_pf_2017, total_pf_2025)
        ],
        index=top_families,
    )

    abs_change_t = (_pf_t_2025 - _pf_t_2017).reindex(top_families, fill_value=0)
    abs_change_m = (_pf_m_2025 - _pf_m_2017).reindex(top_families, fill_value=0)

    # Sort by overall percentage increase (descending) for horizontal bar display
    # This makes the largest increases appear at the top typically for barh
    sorted_families_idx = percentage_increase_pf.sort_values(
        ascending=True
    ).index  # Ascending True for barh from top

    percentage_increase_pf_sorted = percentage_increase_pf.reindex(sorted_families_idx)
    abs_change_t_sorted = abs_change_t.reindex(sorted_families_idx)
    abs_change_m_sorted = abs_change_m.reindex(sorted_families_idx)

    fig, ax = plt.subplots(
        figsize=(12, max(8, len(sorted_families_idx) * 0.5))
    )  # Adjusted figsize for horizontal

    y_pos = np.arange(len(sorted_families_idx))
    bars = ax.barh(
        y_pos, percentage_increase_pf_sorted.values, color="purple", align="center"
    )

    ax.set_xlabel("Overall Percentage Increase (%)", fontsize=14)
    ax.set_ylabel("Protein Family", fontsize=14)
    ax.set_title(title, pad=20, fontsize=16)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(percentage_increase_pf_sorted.index, fontsize=12)
    ax.tick_params(axis="x", labelsize=12)
    # ax.invert_yaxis() # Keep if sort_values ascending=False, not needed for ascending=True for barh

    for i, bar in enumerate(bars):
        perc_inc = percentage_increase_pf_sorted.iloc[i]
        t_abs_c = abs_change_t_sorted.iloc[i]
        m_abs_c = abs_change_m_sorted.iloc[i]

        bar_width = bar.get_width()  # This is the percentage increase value

        text_perc = (
            f"{perc_inc:.1f}%"
            if np.isfinite(perc_inc)
            else "Inf"
            if perc_inc == np.inf
            else "-Inf"
            if perc_inc == -np.inf
            else "N/A"  # For NaN or other non-finite cases
        )
        annotation_text = f"{text_perc}\n(T: {t_abs_c:+.0f}, M: {m_abs_c:+.0f})"

        text_color = "white"
        ha_align = "right"

        if np.isfinite(bar_width):
            if bar_width >= 0:
                text_x_pos = bar_width - (
                    0.05 * bar_width if bar_width > 0 else 0.01
                )  # Offset from end
                ha_align = "right"
            else:  # Negative bar_width
                text_x_pos = bar_width + (
                    0.05 * abs(bar_width) if bar_width < 0 else 0.01
                )  # Offset from end (left)
                ha_align = "left"
        elif bar_width == np.inf:
            # Place text near the right edge of the plot if bar goes to infinity
            text_x_pos = ax.get_xlim()[1] * 0.95
            ha_align = "right"
        elif bar_width == -np.inf:
            # Place text near the left edge of the plot if bar goes to -infinity
            text_x_pos = ax.get_xlim()[0] * 0.95
            ha_align = "left"
        else:  # NaN or other non-finite, place at 0
            text_x_pos = 0
            ha_align = "center"

        # Ensure text is not placed exactly at the bar_width if it is 0, to avoid ha issues.
        if (
            bar_width == 0 and text_x_pos == 0 and ha_align == "right"
        ):  # common case from calculation
            ha_align = (
                "left"  # shift it slightly to the right of the y-axis for 0 values
            )
            text_x_pos = 0.01 * ax.get_xlim()[1] if ax.get_xlim()[1] > 0 else 0.01

        ax.text(
            text_x_pos,
            bar.get_y() + bar.get_height() / 2,
            annotation_text,
            ha=ha_align,
            va="center",
            color=text_color,
            fontsize=10,  # Increased, but check for fit
            linespacing=1.3,
        )
    # Add a line at x=0 if there can be negative values (decreases)
    ax.axvline(0, color="grey", lw=0.8)

    plt.tight_layout(rect=[0.05, 0, 1, 0.95])  # Adjust left for long family names
    filepath = FIGURE_OUTPUT_PATH / "3_percentage_increase_by_protein_family_mpl.png"
    plt.savefig(filepath)
    print(f"Generated: {filepath}")
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
    # Replace inf with a value slightly larger than vmax for coloring purposes if needed,
    # or rely on text annotation. For 0-100% scale, infs will get max color.
    heatmap_plot_data = heatmap_data_numeric.copy()

    # Define categories for the heatmap
    bounds_for_norm = [
        0,
        20,
        40,
        60,
        80,
        np.inf,
    ]  # np.inf will handle actual infinities
    colors = ["lightblue", "cornflowerblue", "royalblue", "mediumblue", "darkblue"]
    cmap_categorical = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(bounds_for_norm, cmap_categorical.N)

    # heatmap_plot_data already contains the correct data, including np.inf where applicable
    # Values > 100 will be colored as 100. Inf will also be colored as 100.
    # The text will show the original value.

    num_rows, num_cols = heatmap_plot_data.shape

    # Adjust figsize for square cells and potentially long labels. Cell size determines base.
    cell_size = 0.8  # Increased cell size for clarity
    # Width needs to account for protein family labels and colorbar
    # Max length of family names for y-tick label space estimation
    max_family_name_len = (
        max(len(name) for name in heatmap_plot_data.index)
        if not heatmap_plot_data.empty
        else 10
    )
    estimated_y_label_width = (
        max_family_name_len * 0.1
    )  # Rough estimate, adjust factor as needed

    fig_width = (
        num_cols * cell_size + estimated_y_label_width + 2.5
    )  # Space for y-labels, colorbar, padding
    fig_height = num_rows * cell_size + 2  # Space for title, x-labels, padding

    fig, ax = plt.subplots(figsize=(max(8, fig_width), max(6, fig_height)))

    cax = ax.imshow(
        heatmap_plot_data,  # Original data with potential np.inf
        cmap=cmap_categorical,  # Use categorical cmap
        norm=norm,  # Use boundary norm
        aspect="equal",
        # vmin and vmax are controlled by the norm now
    )

    # Adjust colorbar to match heatmap height
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes("right", size="5%", pad=0.1)

    # Use finite boundaries for drawing the colorbar segments to avoid issues with np.inf
    # The actual data mapping is handled by the norm passed to imshow.
    cbar_display_bounds = [
        0,
        20,
        40,
        60,
        80,
        100,
    ]  # Last segment (80-100) represents 80-Inf

    cb = fig.colorbar(
        cax,
        cax=cbar_ax,
        boundaries=cbar_display_bounds,  # Finite boundaries for drawing
        ticks=[10, 30, 50, 70, 90],
    )  # Midpoints for labels
    cb.set_label("% Increase Categories")
    cb.set_ticklabels(["0-20%", "20-40%", "40-60%", "60-80%", "80%-Inf"])

    ax.set_xticks(np.arange(num_cols))
    ax.set_yticks(np.arange(num_rows))
    ax.set_xticklabels(heatmap_plot_data.columns)
    ax.set_yticklabels(
        heatmap_plot_data.index
    )  # Should be sorted if top_families is sorted
    ax.set_title(title, pad=20)

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    # Add text annotations (original values and absolute counts)
    for i in range(num_rows):
        for j in range(num_cols):
            family = heatmap_data.index[i]
            habitat = heatmap_data.columns[j]

            perc_val = heatmap_data.iloc[i, j]

            # Get original counts for this cell
            if habitat == "Terrestrial":
                old_count = pf_terrestrial_2017.get(family, 0)
                new_count = pf_terrestrial_2025.get(family, 0)
            else:  # Marine
                old_count = pf_marine_2017.get(family, 0)
                new_count = pf_marine_2025.get(family, 0)

            text_perc = (
                f"{perc_val:.0f}%"
                if np.isfinite(perc_val)
                else "Inf"
                if perc_val == np.inf
                else "0%"
            )
            text_abs = f"({old_count}→{new_count})"
            full_text = f"{text_perc}\n{text_abs}"

            cell_color_val_for_text = heatmap_plot_data.iloc[i, j]
            # Determine text color based on the actual value and chosen colormap
            # We need to find which category the value falls into to pick text color
            # This is a bit more complex with BoundaryNorm. A simple threshold on normalized value might work.
            # Or, determine the color of the cell and then pick contrasting text.

            # For simplicity, let's try a heuristic based on the raw value or its category.
            # If the value is in the 'darkblue' category (last one), use white text.
            # Otherwise, use black (assuming other colors are light enough).
            # The raw percentage value is `perc_val`
            if (
                perc_val >= 80 or perc_val == np.inf
            ):  # Corresponds to the last category (darkblue)
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
                fontsize=8,
                linespacing=1.3,
            )

    plt.tight_layout(
        rect=[0, 0, 0.9, 1]
    )  # Initial tight_layout before specific adjustments
    # Further adjust left margin if y-tick labels are long, after tight_layout has done its job
    fig.subplots_adjust(left=max(0.15, estimated_y_label_width / fig_width + 0.05))

    filepath = FIGURE_OUTPUT_PATH / "5_heatmap_percentage_increase_mpl.png"
    plt.savefig(filepath)
    print(f"Generated: {filepath}")


def plot_diverging_bar_charts_protein_family_mpl():
    # Ensure top_families is not empty
    if not top_families:
        print("Skipping diverging bar charts: No top families identified.")
        return

    # Prepare data for diverging bar chart
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

    # Sort by largest summed absolute increase if flag is set
    if SORT_DIVERGING_BARS_BY_SUMMED_INCREASE:
        protein_family_summary["Summed_Abs_Increase"] = (
            protein_family_summary["Terrestrial_Change_Abs"]
            + protein_family_summary["Marine_Change_Abs"]
        )
        # For diverging, we typically sort by the magnitude of one side or a combined metric.
        # Here, sorting by summed absolute increase.
        # To ensure families with largest positive terrestrial change are on top (if that's desired for right side)
        # or largest marine change (if that's desired for left side)
        # Let's sort by the sum of absolute changes.
        protein_family_summary = protein_family_summary.sort_values(
            by="Summed_Abs_Increase", ascending=False
        )

    # --- Plot 1: Absolute Numbers Change ---
    title_abs = f"Top {len(protein_family_summary)} Protein Families: Absolute Change (2017-2025, Dual-Habitat)"
    y_labels_abs = protein_family_summary.index
    terrestrial_values_abs = protein_family_summary["Terrestrial_Change_Abs"]
    marine_values_abs = protein_family_summary[
        "Marine_Change_Abs"
    ]  # Positive values for left side plotting

    fig_abs, ax_abs = plt.subplots(
        figsize=(12, max(6, len(y_labels_abs) * 0.4))
    )  # Dynamic height
    y_pos_abs = np.arange(len(y_labels_abs))

    # Terrestrial (right side, positive)
    bars_t_abs = ax_abs.barh(
        y_pos_abs,
        terrestrial_values_abs,
        color="forestgreen",
        label="Terrestrial Increase/Decrease",
    )
    # Marine (left side, plot positive values but on negative axis by negating them)
    bars_m_abs = ax_abs.barh(
        y_pos_abs,
        -marine_values_abs,
        color="steelblue",
        label="Marine Increase/Decrease (shown left)",
    )

    ax_abs.set_yticks(y_pos_abs)
    ax_abs.set_yticklabels(y_labels_abs)
    ax_abs.invert_yaxis()  # Highest summed increase at the top
    ax_abs.set_xlabel("Absolute Change in Number of Sequences")
    ax_abs.set_title(title_abs, pad=20)
    ax_abs.legend(loc="lower right")
    ax_abs.axvline(0, color="grey", lw=0.8)  # Add a line at x=0

    # Add text annotations for absolute numbers
    for i, (t_val, m_val) in enumerate(zip(terrestrial_values_abs, marine_values_abs)):
        # Terrestrial text (right)
        if t_val != 0:
            ax_abs.text(
                t_val
                + (
                    0.01 * ax_abs.get_xlim()[1]
                    if t_val > 0
                    else -0.03 * ax_abs.get_xlim()[1]
                ),  # Dynamic offset
                y_pos_abs[i],
                f"{t_val}",
                va="center",
                ha="left" if t_val > 0 else "right",
                fontsize=8,
            )
        # Marine text (left, around the negated bar)
        if m_val != 0:
            ax_abs.text(
                -m_val
                - (
                    0.01 * ax_abs.get_xlim()[1]
                    if m_val > 0
                    else -0.03 * ax_abs.get_xlim()[1]
                ),  # Dynamic offset for negated bar
                y_pos_abs[i],
                f"{m_val}",
                va="center",
                ha="right"
                if m_val > 0
                else "left",  # ha is 'right' for positive m_val on left
                fontsize=8,
            )

    # Adjust x-limits to ensure text is visible
    current_xlim_abs = ax_abs.get_xlim()
    expansion_factor_abs = 0.25  # Increase to 25%
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

    # --- Plot 2: Percentage Change ---
    title_perc = f"Top {len(protein_family_summary)} Protein Families: Percentage Change (2017-2025, Dual-Habitat)"
    y_labels_perc = protein_family_summary.index  # Same sorted order

    # Handle inf/-inf for plotting by replacing with a large number or NaN
    # Let's use a large number and note it in text. np.inf won't plot on barh directly.
    # Or, we can cap them. For now, let's plot what we have and handle text.
    terrestrial_values_perc = protein_family_summary["Terrestrial_Change_Perc"].replace(
        [np.inf, -np.inf],
        100,  # Cap inf for plotting at 100%
    )
    marine_values_perc = protein_family_summary["Marine_Change_Perc"].replace(
        [np.inf, -np.inf],
        100,  # Cap inf for plotting at 100%
    )

    # Store original values for text display
    terrestrial_original_perc = protein_family_summary["Terrestrial_Change_Perc"]
    marine_original_perc = protein_family_summary["Marine_Change_Perc"]

    fig_perc, ax_perc = plt.subplots(figsize=(12, max(6, len(y_labels_perc) * 0.4)))
    y_pos_perc = np.arange(len(y_labels_perc))

    bars_t_perc = ax_perc.barh(
        y_pos_perc,
        terrestrial_values_perc,
        color="lightcoral",
        label="Terrestrial % Change",
    )
    bars_m_perc = ax_perc.barh(
        y_pos_perc,
        -marine_values_perc,
        color="lightskyblue",
        label="Marine % Change (shown left)",
    )

    ax_perc.set_yticks(y_pos_perc)
    ax_perc.set_yticklabels(y_labels_perc)
    ax_perc.invert_yaxis()
    ax_perc.set_xlabel(
        "Percentage Change (%) [Values > +/-100% capped for vis]"
    )  # Updated label
    ax_perc.set_title(title_perc, pad=20)
    ax_perc.legend(loc="lower right")

    # Add text annotations for percentage change
    for i, (t_orig_perc, m_orig_perc, t_plot_perc, m_plot_perc) in enumerate(
        zip(
            terrestrial_original_perc,
            marine_original_perc,
            terrestrial_values_perc,
            marine_values_perc,
        )
    ):
        # Terrestrial text (right)
        t_text = (
            "Inf"
            if np.isinf(t_orig_perc)
            else f"{t_orig_perc:.1f}%"
            if pd.notna(t_orig_perc)
            else "N/A"
        )
        if t_plot_perc != 0 or np.isinf(
            t_orig_perc
        ):  # Check plot value or if original was inf
            ax_perc.text(
                t_plot_perc
                + (
                    0.01 * ax_perc.get_xlim()[1]
                    if t_plot_perc > 0
                    else -0.03 * ax_perc.get_xlim()[1]
                ),
                y_pos_perc[i],
                t_text,
                va="center",
                ha="left" if t_plot_perc > 0 else "right",
                fontsize=8,
            )

        # Marine text (left)
        m_text = (
            "Inf"
            if np.isinf(m_orig_perc)
            else f"{m_orig_perc:.1f}%"
            if pd.notna(m_orig_perc)
            else "N/A"
        )
        if m_plot_perc != 0 or np.isinf(
            m_orig_perc
        ):  # Check plot value or if original was inf
            ax_perc.text(
                -m_plot_perc
                - (
                    0.01 * ax_perc.get_xlim()[1]
                    if m_plot_perc > 0
                    else -0.03 * ax_perc.get_xlim()[1]
                ),
                y_pos_perc[i],
                m_text,
                va="center",
                ha="right" if m_plot_perc > 0 else "left",
                fontsize=8,
            )

    # Adjust x-limits to ensure text is visible
    current_xlim_perc = ax_perc.get_xlim()
    expansion_factor_perc = 0.25  # Increase to 25%
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


# def plot_sankey_diagram_mpl():


def generate_dual_habitat_protein_counts_summary():
    """
    Generates a CSV file summarizing protein counts per habitat, year, and protein family,
    considering only families present in both terrestrial and marine habitats (dual-habitat).
    The output format will be: ProteinFamily, 2017_Terrestrial, 2017_Marine, 2025_Terrestrial, 2025_Marine.
    """
    print("Generating summary file of dual-habitat protein counts in wide format...")

    # These Series (pf_terrestrial_2017, etc.) are already filtered for dual-habitat families
    # and contain counts for each family within that specific year and habitat.

    s_t_2017 = pf_terrestrial_2017.rename("2017_Terrestrial")
    s_m_2017 = pf_marine_2017.rename("2017_Marine")
    s_t_2025 = pf_terrestrial_2025.rename("2025_Terrestrial")
    s_m_2025 = pf_marine_2025.rename("2025_Marine")

    # Concatenate the series. The index (ProteinFamily) will align them.
    summary_df = pd.concat([s_t_2017, s_m_2017, s_t_2025, s_m_2025], axis=1)

    # Fill NaN with 0 (for families that might be in one series but not another, though unlikely here)
    # and convert counts to integers.
    summary_df = summary_df.fillna(0).astype(int)

    # Reset index to make 'ProteinFamily' a column
    # The original index name is 'Protein families' from the value_counts() source
    summary_df = summary_df.reset_index().rename(
        columns={"Protein families": "ProteinFamily"}
    )

    if summary_df.empty:
        print(
            "No data to summarize for dual-habitat protein counts file. Skipping file generation."
        )
        return

    # Sort for consistency and readability
    summary_df = summary_df.sort_values(by="ProteinFamily").reset_index(drop=True)

    output_filename = "dual_habitat_protein_counts_summary.csv"
    output_filepath = FIGURE_OUTPUT_PATH / output_filename
    summary_df.to_csv(output_filepath, index=False)
    print(f"Generated summary file: {output_filepath}")


# --- Main Execution ---
if __name__ == "__main__":
    # ... (File I/O test and initial checks remain the same) ...
    # test_file_path = FIGURE_OUTPUT_PATH / "_io_test.txt" # Removed I/O Test
    # try:
    #     with open(test_file_path, "w") as f:
    #         f.write("Basic file I/O test successful.")
    #     print(f"Successfully wrote test file to: {test_file_path}")
    #     # test_file_path.unlink()
    # except Exception as e:
    #     print(f"ERROR: Could not write a test file to {FIGURE_OUTPUT_PATH}. Error: {e}")

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
        generate_dual_habitat_protein_counts_summary()  # Generate the summary file
        plot_100_stacked_bar_habitat_vs_year_mpl()
        plot_percentage_increase_by_habitat_mpl()
        plot_percentage_increase_by_protein_family_mpl()
        # plot_sankey_diagrams_mpl() # Sankey plots are commented out
        plot_heatmap_percentage_increase_mpl()
        plot_diverging_bar_charts_protein_family_mpl()
        print(
            f"All Matplotlib plots and summary file attempted to be saved to {FIGURE_OUTPUT_PATH}"
        )
