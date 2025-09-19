import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# --- Configuration ---
TOP_N_GO_TERMS = 15

# --- Path Definitions ---
try:
    BASE_PATH = Path(__file__).resolve().parent.parent
except NameError:  # Fallback for interactive environments
    BASE_PATH = (Path.cwd()) if Path.cwd().name == "toxprot25" else Path.cwd() / ".."

DATA_PATH = BASE_PATH / "data" / "processed"
FIGURE_OUTPUT_PATH = BASE_PATH / "figures" / "go_terms"
FIGURE_OUTPUT_PATH.mkdir(parents=True, exist_ok=True)

# --- Data Loading ---
try:
    toxprot_2017_path = DATA_PATH / "toxprot_2017.csv"
    toxprot_2025_path = DATA_PATH / "toxprot_2025.csv"

    df_2017 = pd.read_csv(toxprot_2017_path)
    df_2025 = pd.read_csv(toxprot_2025_path)

except FileNotFoundError as e:
    print(f"Error loading data: {e}. Please ensure CSV files are in {DATA_PATH}")
    exit()


# --- GO Term Analysis Functions ---
def get_go_counts(df, go_col="Gene Ontology (GO)"):
    """Extract and count GO terms from a DataFrame."""
    # Dropna, split by <; >, flatten, and count
    go_series = df[go_col].dropna().astype(str)
    all_terms = []
    for entry in go_series:
        terms = [t.strip() for t in entry.split(";") if t.strip()]
        all_terms.extend(terms)
    return pd.Series(all_terms).value_counts(), all_terms, go_series


def calculate_percentage_change(old_val, new_val):
    """Safely calculate percentage change, handling division by zero."""
    if old_val == 0:
        return np.inf if new_val > 0 else 0.0
    return ((new_val - old_val) / old_val) * 100


def get_go_counts_for_plotting():
    """Get GO term counts for plotting without printing summary statistics."""
    go_counts_2025, _, _ = get_go_counts(df_2025)
    go_counts_2017, _, _ = get_go_counts(df_2017)
    return go_counts_2017, go_counts_2025


def plot_go_terms_stacked_bar():
    """Create a stacked bar plot showing GO term counts with 2017 as base and increase stacked on top."""

    # Get GO term counts
    go_counts_2017, go_counts_2025 = get_go_counts_for_plotting()

    # Get top N GO terms by 2025 count (current usage)
    top_go_terms = go_counts_2025.head(TOP_N_GO_TERMS).index.tolist()

    if not top_go_terms:
        print("No GO terms found for plotting.")
        return

    # Create summary DataFrame
    go_summary = pd.DataFrame(index=top_go_terms)
    go_summary["2017_count"] = go_counts_2017.reindex(top_go_terms, fill_value=0)
    go_summary["2025_count"] = go_counts_2025.reindex(top_go_terms, fill_value=0)
    go_summary["absolute_change"] = go_summary["2025_count"] - go_summary["2017_count"]
    go_summary["percentage_change"] = [
        calculate_percentage_change(old, new)
        for old, new in zip(go_summary["2017_count"], go_summary["2025_count"])
    ]

    # Sort by 2025 count (descending) to show terms with highest usage at top
    go_summary = go_summary.sort_values(by="2025_count", ascending=False)

    # Create the plot
    title = f"Top {len(go_summary)} GO Terms"
    y_labels = go_summary.index

    values_2017 = go_summary["2017_count"]
    values_2025 = go_summary["2025_count"]
    # Calculate the increase (will be 0 if count decreased)
    values_increase = np.maximum(0, values_2025 - values_2017)
    # For stacking, if there was a decrease, we show the full 2025 count as the base
    values_base = np.minimum(values_2017, values_2025)

    fig, ax = plt.subplots(figsize=(12, max(8, len(y_labels) * 0.5)))
    y_pos = np.arange(len(y_labels))

    # Define colors with better differentiation
    color_base = "steelblue"  # Base (existing) portion - darker blue
    color_increase = "#87CEEB"  # Sky blue for increase - lighter, more distinguishable

    # Create horizontal stacked bars
    # Base portion (2017 count)
    ax.barh(y_pos, values_base, color=color_base, label="2017", alpha=0.8)

    # Increase portion (stacked on top of base)
    ax.barh(
        y_pos,
        values_increase,
        left=values_base,
        color=color_increase,
        label="2025",
        alpha=0.8,
    )

    # Formatting
    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [term[:60] + "..." if len(term) > 60 else term for term in y_labels],
        fontsize=10,
    )
    ax.invert_yaxis()  # Highest change at the top
    ax.set_xlabel("Number of GO Term Assignments")
    ax.set_title(title, pad=20, fontsize=14)
    ax.legend(loc="lower right")

    # Add text annotations showing the counts and changes
    for i, go_term in enumerate(y_labels):
        count_2017 = values_2017[go_term]
        count_2025 = values_2025[go_term]

        # Position text after the end of the stacked bar (2025 count)
        text_x_position = count_2025 + 0.01 * ax.get_xlim()[1]

        # Create annotation text with "New" indicator for terms not in 2017
        if count_2017 == 0:
            annotation_text = f"{count_2017}→{count_2025} (New)"
        else:
            annotation_text = f"{count_2017}→{count_2025}"

        # Annotation showing the progression
        ax.text(
            text_x_position,
            y_pos[i],
            annotation_text,
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold",
        )

    # Adjust x-limits to ensure text is visible
    current_xlim = ax.get_xlim()
    expansion_factor = 0.25
    ax.set_xlim(
        current_xlim[0],
        current_xlim[1] + abs(current_xlim[1]) * expansion_factor,
    )

    plt.tight_layout()

    # Save the plot
    filepath = FIGURE_OUTPUT_PATH / "go_terms_stacked_comparison_2017_vs_2025.png"
    plt.savefig(filepath, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Save summary data
    summary_filepath = FIGURE_OUTPUT_PATH / "go_terms_summary.csv"
    go_summary_export = go_summary.copy()
    go_summary_export.index.name = "GO_Term"
    go_summary_export.to_csv(summary_filepath)


def generate_go_counts_comparison_table():
    """Generate a detailed comparison table of GO term counts."""

    go_counts_2017, go_counts_2025 = get_go_counts(df_2017), get_go_counts(df_2025)
    go_counts_2017, go_counts_2025 = (
        go_counts_2017[0],
        go_counts_2025[0],
    )  # Extract just the counts

    # Get all unique GO terms
    all_terms = set(go_counts_2017.index) | set(go_counts_2025.index)

    # Create comparison DataFrame
    comparison_df = pd.DataFrame(index=sorted(all_terms))
    comparison_df["2017_count"] = go_counts_2017.reindex(
        comparison_df.index, fill_value=0
    )
    comparison_df["2025_count"] = go_counts_2025.reindex(
        comparison_df.index, fill_value=0
    )
    comparison_df["absolute_change"] = (
        comparison_df["2025_count"] - comparison_df["2017_count"]
    )
    comparison_df["percentage_change"] = [
        calculate_percentage_change(old, new)
        for old, new in zip(comparison_df["2017_count"], comparison_df["2025_count"])
    ]

    # Add total counts for percentage calculations
    total_2017 = comparison_df["2017_count"].sum()
    total_2025 = comparison_df["2025_count"].sum()

    comparison_df["2017_percentage"] = (
        (comparison_df["2017_count"] / total_2017 * 100) if total_2017 > 0 else 0
    )
    comparison_df["2025_percentage"] = (
        (comparison_df["2025_count"] / total_2025 * 100) if total_2025 > 0 else 0
    )

    # Sort by 2025 count (descending)
    comparison_df = comparison_df.sort_values(by="2025_count", ascending=False)

    # Save full comparison
    comparison_filepath = FIGURE_OUTPUT_PATH / "go_terms_full_comparison.csv"
    comparison_df.index.name = "GO_Term"
    comparison_df.to_csv(comparison_filepath)

    return comparison_df


# --- Main Execution ---
if __name__ == "__main__":
    # Generate plots and data files
    try:
        # Generate the main stacked bar plot
        plot_go_terms_stacked_bar()

        # Generate detailed comparison table
        generate_go_counts_comparison_table()

    except Exception as e:
        print(f"Error during analysis: {e}")
        raise
