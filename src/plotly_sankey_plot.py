import pandas as pd
import plotly.graph_objects as go
from pathlib import Path

# Define file paths
BASE_PATH = Path(__file__).resolve().parent.parent
DATA_PATH = BASE_PATH / "data" / "processed"
FIGURE_PATH = BASE_PATH / "figures" / "taxa"
FIGURE_PATH.mkdir(parents=True, exist_ok=True)

# Load data
try:
    toxprot_2017_df = pd.read_csv(DATA_PATH / "toxprot_2017.csv")
    toxprot_2025_df = pd.read_csv(DATA_PATH / "toxprot_2025.csv")
except FileNotFoundError as e:
    print(f"Error loading data: {e}")
    print(
        f"Please ensure the files toxprot_2017.csv and toxprot_2025.csv exist in {DATA_PATH}"
    )
    exit()


def get_protein_family_by_habitat(df, habitat_type, split_char):
    """Extract protein family counts for a specific habitat type"""
    habitat_data = df[df["Habitat"] == habitat_type].copy()
    # Ensure 'Protein families' is string type before splitting
    habitat_data["Protein families"] = (
        habitat_data["Protein families"].astype(str).str.split(split_char).str[0]
    )
    return habitat_data["Protein families"].value_counts().sort_index()


# Get data for both years and habitats
terrestrial_2017 = get_protein_family_by_habitat(toxprot_2017_df, "terrestrial", ".")
marine_2017 = get_protein_family_by_habitat(toxprot_2017_df, "marine", ".")
terrestrial_2025 = get_protein_family_by_habitat(toxprot_2025_df, "terrestrial", ",")
marine_2025 = get_protein_family_by_habitat(toxprot_2025_df, "marine", ",")

# Find protein families present in both habitats across both years (or either for broader scope)
# For Sankey, we often want to see flows, so common families are a good start.
# If you want all families from selected top_n, the logic below needs adjustment.
common_families_2017 = set(terrestrial_2017.index) & set(marine_2017.index)
common_families_2025 = set(terrestrial_2025.index) & set(marine_2025.index)
# Or, consider all families present in either terrestrial or marine for a given year, then take top N
all_families_2017 = set(terrestrial_2017.index) | set(marine_2017.index)
all_families_2025 = set(terrestrial_2025.index) | set(marine_2025.index)

# Union of common families across both years to have a consistent set of top families if desired
# both_habitats_overall = common_families_2017 | common_families_2025
# Using all families from both years and then picking top N might be more representative
# For simplicity, let's use the original notebook's logic for top_families for now:
# (set(terrestrial_2017.index) & set(marine_2017.index)) | (set(terrestrial_2025.index) & set(marine_2025.index))

combined_top_families_set = (set(terrestrial_2017.index) & set(marine_2017.index)) | (
    set(terrestrial_2025.index) & set(marine_2025.index)
)

top_n = 15
# Ensure top_families are sorted for consistent node ordering if used across plots
top_families = sorted(list(combined_top_families_set))[:top_n]


def prepare_plotly_sankey_data(terrestrial_data, marine_data, current_top_families):
    """Prepare data for Plotly Sankey diagram."""
    labels = ["Terrestrial", "Marine"] + current_top_families
    source_indices = []
    target_indices = []
    values = []

    node_map = {label: i for i, label in enumerate(labels)}

    # Terrestrial flows
    for family in current_top_families:
        if family in terrestrial_data.index:
            source_indices.append(node_map["Terrestrial"])
            target_indices.append(node_map[family])
            values.append(terrestrial_data[family])

    # Marine flows
    for family in current_top_families:
        if family in marine_data.index:
            source_indices.append(node_map["Marine"])
            target_indices.append(node_map[family])
            values.append(marine_data[family])

    return dict(
        source=source_indices, target=target_indices, value=values, label=labels
    )


def create_and_save_sankey(sankey_data, year_label, output_filename):
    """Creates a Sankey diagram using Plotly and saves it."""
    fig = go.Figure(
        data=[
            go.Sankey(
                node=dict(
                    pad=15,
                    thickness=20,
                    line=dict(color="black", width=0.5),
                    label=sankey_data["label"],
                    # Add colors if desired, e.g., based on node type
                    # color=["blue", "green"] + ["gray"] * len(top_families)
                ),
                link=dict(
                    source=sankey_data["source"],
                    target=sankey_data["target"],
                    value=sankey_data["value"],
                ),
            )
        ]
    )

    fig.update_layout(
        title_text=f"Protein Family Distribution by Habitat ({year_label})",
        font_size=10,
    )
    fig.write_html(FIGURE_PATH / output_filename)
    print(f"Saved Sankey diagram to {FIGURE_PATH / output_filename}")


# Generate and save Sankey for 2017
# We might want different top families per year, or a combined set.
# Using the globally defined `top_families` for consistency:
sankey_data_2017 = prepare_plotly_sankey_data(
    terrestrial_2017, marine_2017, top_families
)
if sankey_data_2017["source"]:  # Check if there's any data to plot
    create_and_save_sankey(
        sankey_data_2017, "2017", "habitat_protein_families_plotly_sankey_2017.html"
    )
else:
    print("No data to plot for 2017 with the current top_families selection.")


# Generate and save Sankey for 2025
# Using the globally defined `top_families` for consistency:
sankey_data_2025 = prepare_plotly_sankey_data(
    terrestrial_2025, marine_2025, top_families
)
if sankey_data_2025["source"]:  # Check if there's any data to plot
    create_and_save_sankey(
        sankey_data_2025, "2025", "habitat_protein_families_plotly_sankey_2025.html"
    )
else:
    print("No data to plot for 2025 with the current top_families selection.")

print("Script completed.")
