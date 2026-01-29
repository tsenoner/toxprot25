#!/usr/bin/env python3
"""
Protein Evidence Sankey Plot - Final Version
============================================

Creates a Sankey diagram showing the flow of protein evidence categories
from ToxProt 2017 to ToxProt 2025 with proper node positioning.

Author: Tobias Senoner
"""

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go


def normalize_pe_category(category):
    """
    Normalize protein existence categories to standard names.

    Args:
        category (str): Raw category from dataset

    Returns:
        str: Normalized category name
    """
    if pd.isna(category):
        return "Unknown"

    category = str(category).strip()

    # Map numbered categories to readable names
    if category.startswith("1:") or "protein level" in category.lower():
        return "Evidence at protein level"
    elif category.startswith("2:") or "transcript level" in category.lower():
        return "Evidence at transcript level"
    elif category.startswith("3:") or "homology" in category.lower():
        return "Inferred from homology"
    elif category.startswith("4:") or "predicted" in category.lower():
        return "Predicted"
    elif category.startswith("5:") or "uncertain" in category.lower():
        return "Uncertain"
    else:
        return category


def hex_to_rgba(hex_color, alpha=0.4):
    """Convert hex color to rgba with transparency."""
    hex_color = hex_color.lstrip("#")
    if len(hex_color) != 6:
        return f"rgba(200,200,200,{alpha})"
    try:
        r, g, b = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
        return f"rgba({r},{g},{b},{alpha})"
    except ValueError:
        return f"rgba(200,200,200,{alpha})"


def calculate_node_positions(categories_2017, categories_2025, flow_counter, has_new_node=False):
    """
    Calculate proper Y positions for nodes to avoid overlap.

    Args:
        categories_2017: List of categories for 2017 side
        categories_2025: List of categories for 2025 side
        flow_counter: Counter of flows between categories
        has_new_node: Boolean indicating if there's a "New" node in the middle

    Returns:
        tuple: (y_positions_2017, y_positions_2025, y_position_new)
    """
    # Calculate total flow for each node to determine relative sizes
    node_flows_2017 = {cat: 0 for cat in categories_2017}
    node_flows_2025 = {cat: 0 for cat in categories_2025}
    new_node_flow = 0

    for (src, tgt), count in flow_counter.items():
        if src == "New":
            new_node_flow += count
        elif src in node_flows_2017:
            node_flows_2017[src] += count
        if tgt in node_flows_2025:
            node_flows_2025[tgt] += count

    # Calculate positions for 2017 side - leave space at bottom for New node
    def calculate_2017_positions(categories, node_flows):
        """Calculate y positions for 2017 side, leaving bottom space."""
        if not categories:
            return []

        # Calculate total flow
        total_flow = sum(node_flows[cat] for cat in categories)
        if total_flow == 0:
            # If no flow, distribute evenly
            n = len(categories)
            return [0.05 + i * 0.5 / max(n - 1, 1) for i in range(n)]

        # Use only top portion of space to leave room for New node at bottom
        margin_top = 0.05  # 5% margin at top
        margin_bottom = 0.35  # 35% margin at bottom for New node
        usable_height = 1.0 - margin_top - margin_bottom

        # Minimum spacing between nodes
        min_spacing = 0.015  # 1.5% minimum gap
        total_spacing = min_spacing * (len(categories) - 1)

        # Available height for nodes after accounting for spacing
        available_for_nodes = usable_height - total_spacing

        positions = []
        current_y = margin_top

        for _i, cat in enumerate(categories):
            # Calculate this node's height based on its flow
            node_height = (node_flows[cat] / total_flow) * available_for_nodes

            # Position is at the center of the node
            y_position = current_y + node_height / 2

            # Ensure position is within valid range
            y_position = max(0.001, min(0.999, y_position))
            positions.append(y_position)

            # Move to next node position
            current_y += node_height + min_spacing

        return positions

    # Calculate positions for 2025 side - use full height
    def calculate_2025_positions(categories, node_flows):
        """Calculate y positions for 2025 side using full height."""
        if not categories:
            return []

        # Calculate total flow
        total_flow = sum(node_flows[cat] for cat in categories)
        if total_flow == 0:
            # If no flow, distribute evenly
            n = len(categories)
            return [0.05 + i * 0.9 / max(n - 1, 1) for i in range(n)]

        margin_top = 0.01  # 1% margin at top
        margin_bottom = 0.05  # 5% margin at bottom
        usable_height = 1.0 - margin_top - margin_bottom

        # Minimum spacing between nodes
        min_spacing = 0.02  # 2% minimum gap
        total_spacing = min_spacing * (len(categories) - 1)

        # Available height for nodes after accounting for spacing
        available_for_nodes = usable_height - total_spacing

        positions = []
        current_y = margin_top

        for _i, cat in enumerate(categories):
            # Calculate this node's height based on its flow
            node_height = (node_flows[cat] / total_flow) * available_for_nodes

            # Position is at the center of the node
            y_position = current_y + node_height / 2

            # Ensure position is within valid range
            y_position = max(0.001, min(0.999, y_position))
            positions.append(y_position)

            # Move to next node position
            current_y += node_height + min_spacing

        return positions

    y_positions_2017 = calculate_2017_positions(categories_2017, node_flows_2017)
    y_positions_2025 = calculate_2025_positions(categories_2025, node_flows_2025)

    # Position for the "New" node - at the bottom, with some clearance
    y_position_new = 0.95 if has_new_node else None

    return y_positions_2017, y_positions_2025, y_position_new


def create_protein_evidence_sankey(toxprot_2017_path, toxprot_2025_path, output_path=None):
    """
    Create a Sankey diagram showing protein evidence category transitions.

    Args:
        toxprot_2017_path (str): Path to 2017 ToxProt data
        toxprot_2025_path (str): Path to 2025 ToxProt data
        output_path (str, optional): Path to save the plot

    Returns:
        go.Figure: Plotly figure object
    """

    # Load data
    print("Loading data...")
    df_2017 = pd.read_csv(toxprot_2017_path)
    df_2025 = pd.read_csv(toxprot_2025_path)

    # Determine protein ID column
    id_col = "Entry" if "Entry" in df_2017.columns else df_2017.columns[0]

    print(f"Using '{id_col}' as protein ID column")
    print(f"2017 dataset: {len(df_2017)} proteins")
    print(f"2025 dataset: {len(df_2025)} proteins")

    # Normalize protein existence categories
    df_2017["PE_normalized"] = df_2017["Protein existence"].apply(normalize_pe_category)
    df_2025["PE_normalized"] = df_2025["Protein existence"].apply(normalize_pe_category)

    # Map protein ID to normalized PE for each year
    pe_2017 = df_2017.set_index(id_col)["PE_normalized"]
    pe_2025 = df_2025.set_index(id_col)["PE_normalized"]

    print(f"2017 categories: {pe_2017.value_counts().to_dict()}")
    print(f"2025 categories: {pe_2025.value_counts().to_dict()}")

    # Define canonical PE categories and colors
    pe_categories = [
        "Evidence at protein level",
        "Evidence at transcript level",
        "Inferred from homology",
        "Predicted",
        "Uncertain",
    ]

    # Define colors - using a more distinct palette
    pe_colors = [
        "#2E86AB",  # Blue - Evidence at protein level
        "#A23B72",  # Purple - Evidence at transcript level
        "#F18F01",  # Orange - Inferred from homology
        "#C73E1D",  # Red - Predicted
        "#592941",  # Dark purple - Uncertain
    ]

    special_categories = ["Removed", "New"]
    special_colors = ["#808080", "#A0A0A0"]  # Gray shades for removed/new

    # Find all categories present in the data
    all_categories_2017 = set(pe_2017.dropna().unique())
    all_categories_2025 = set(pe_2025.dropna().unique())
    all_categories = all_categories_2017.union(all_categories_2025)

    # Order categories: canonical first, then others
    ordered_categories = []
    for cat in pe_categories:
        if cat in all_categories:
            ordered_categories.append(cat)

    # Add any extra categories
    for cat in sorted(all_categories):
        if cat not in ordered_categories:
            ordered_categories.append(cat)

    print(f"All categories: {ordered_categories}")

    # Build color mapping
    color_map = {}
    for _i, cat in enumerate(ordered_categories):
        if cat in pe_categories:
            idx = pe_categories.index(cat)
            color_map[cat] = pe_colors[idx]
        elif cat in special_categories:
            idx = special_categories.index(cat)
            color_map[cat] = special_colors[idx]
        else:
            color_map[cat] = "#CCCCCC"  # Default gray

    # Calculate flows
    print("Calculating protein flows...")
    flow_counter = Counter()

    # Proteins present in both years
    common_ids = set(pe_2017.index) & set(pe_2025.index)
    print(f"Common proteins: {len(common_ids)}")

    for pid in common_ids:
        src = pe_2017[pid]
        tgt = pe_2025[pid]
        if pd.notna(src) and pd.notna(tgt):
            flow_counter[(src, tgt)] += 1

    # Proteins only in 2017 (removed in 2025)
    removed_ids = set(pe_2017.index) - set(pe_2025.index)
    print(f"Removed proteins: {len(removed_ids)}")

    for pid in removed_ids:
        src = pe_2017[pid]
        if pd.notna(src):
            flow_counter[(src, "Removed")] += 1

    # Proteins only in 2025 (new in 2025)
    new_ids = set(pe_2025.index) - set(pe_2017.index)
    print(f"New proteins: {len(new_ids)}")

    for pid in new_ids:
        tgt = pe_2025[pid]
        if pd.notna(tgt):
            flow_counter[("New", tgt)] += 1

    # Determine which categories appear on each side
    flow_categories_2017 = set()  # Source categories (2017 side)
    flow_categories_2025 = set()  # Target categories (2025 side)
    has_new_node = False

    for (src, tgt), _count in flow_counter.items():
        if src == "New":
            has_new_node = True
        else:
            flow_categories_2017.add(src)
        flow_categories_2025.add(tgt)

    # Create separate category lists for each side in a logical order
    categories_2017 = []
    categories_2025 = []

    # Add regular categories that appear on each side
    for cat in ordered_categories:
        if cat in flow_categories_2017:
            categories_2017.append(cat)
        if cat in flow_categories_2025:
            categories_2025.append(cat)

    # Add "Removed" to 2025 side if present
    if "Removed" in flow_categories_2025:
        categories_2025.append("Removed")

    # Update color mapping for any new categories
    all_unique_categories = set(categories_2017 + categories_2025)
    if has_new_node:
        all_unique_categories.add("New")
    for cat in all_unique_categories:
        if cat not in color_map:
            if cat in special_categories:
                idx = special_categories.index(cat)
                color_map[cat] = special_colors[idx]
            else:
                color_map[cat] = "#CCCCCC"  # Default gray

    print(f"Total flows: {len(flow_counter)}")
    print(f"2017 categories: {categories_2017}")
    print(f"2025 categories: {categories_2025}")
    print(f"Has 'New' node: {has_new_node}")

    # Calculate proper Y positions
    y_positions_2017, y_positions_2025, y_position_new = calculate_node_positions(
        categories_2017, categories_2025, flow_counter, has_new_node
    )

    # Calculate node counts and percentages for labels
    node_counts_2017 = {cat: 0 for cat in categories_2017}
    node_counts_2025 = {cat: 0 for cat in categories_2025}
    new_node_count = 0

    for (src, tgt), count in flow_counter.items():
        if src == "New":
            new_node_count += count
        elif src in node_counts_2017:
            node_counts_2017[src] += count
        if tgt in node_counts_2025:
            node_counts_2025[tgt] += count

    # Calculate totals for percentage calculations
    # For 2017: only count proteins that were in 2017
    total_2017 = sum(node_counts_2017.values())
    # For 2025: count all proteins in 2025 (including new ones)
    total_2025 = sum(node_counts_2025.values())

    # Build Sankey diagram data
    # Create node labels with protein counts and percentages
    labels_2017 = []
    for cat in categories_2017:
        count = node_counts_2017[cat]
        percentage = (count / total_2017 * 100) if total_2017 > 0 else 0
        if percentage < 1.0:
            labels_2017.append(f"{cat} {count:,} proteins ({percentage:.1f}%)")
        else:
            labels_2017.append(f"{cat}<br>{count:,} proteins ({percentage:.1f}%)")

    labels_2025 = []
    for cat in categories_2025:
        count = node_counts_2025[cat]
        percentage = (count / total_2025 * 100) if total_2025 > 0 else 0
        if percentage < 1.0:
            labels_2025.append(f"{cat} {count:,} proteins ({percentage:.1f}%)")
        else:
            labels_2025.append(f"{cat}<br>{count:,} proteins ({percentage:.1f}%)")

    # Build all labels list
    all_labels = labels_2017 + labels_2025
    if has_new_node:
        all_labels.append(f"New in 2025<br>{new_node_count:,} proteins")

    # Create node colors
    colors_2017 = [color_map[cat] for cat in categories_2017]
    colors_2025 = [color_map[cat] for cat in categories_2025]
    node_colors = colors_2017 + colors_2025
    if has_new_node:
        node_colors.append(color_map["New"])

    # Create node index mappings
    cat_to_2017_idx = {cat: i for i, cat in enumerate(categories_2017)}
    cat_to_2025_idx = {cat: i + len(categories_2017) for i, cat in enumerate(categories_2025)}
    if has_new_node:
        new_node_idx = len(categories_2017) + len(categories_2025)

    # Build links
    sources, targets, values, link_colors = [], [], [], []

    for (src, tgt), count in flow_counter.items():
        if src == "New" and has_new_node:
            sources.append(new_node_idx)
            targets.append(cat_to_2025_idx[tgt])
            values.append(count)
            link_colors.append(hex_to_rgba(color_map["New"], 0.5))
        elif src in cat_to_2017_idx and tgt in cat_to_2025_idx:
            sources.append(cat_to_2017_idx[src])
            targets.append(cat_to_2025_idx[tgt])
            values.append(count)
            link_colors.append(hex_to_rgba(color_map[src], 0.5))

    print(f"Created {len(sources)} links")

    # X positions: 2017 on left, 2025 on right, New in middle
    x_positions_2017 = [0.001] * len(categories_2017)
    x_positions_2025 = [0.999] * len(categories_2025)
    all_x_positions = x_positions_2017 + x_positions_2025
    if has_new_node:
        all_x_positions.append(0.35)  # Middle position for "New" node

    # Combine all Y positions
    all_y_positions = y_positions_2017 + y_positions_2025
    if has_new_node:
        all_y_positions.append(y_position_new)

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="snap",  # Use snap for some automatic adjustment
                node=dict(
                    pad=15,  # Padding between nodes
                    thickness=25,  # Node thickness
                    line=dict(color="black", width=1),
                    label=all_labels,
                    color=node_colors,
                    x=all_x_positions,
                    y=all_y_positions,
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=link_colors,
                ),
                textfont=dict(size=24, color="black"),  # Significantly increased font size
            )
        ]
    )

    # Add year labels as annotations above each column
    fig.add_annotation(
        x=0.001,
        y=1.08,
        text="<b>2017</b>",
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(size=24, color="black"),
        xanchor="center",
    )

    fig.add_annotation(
        x=0.999,
        y=1.08,
        text="<b>2025</b>",
        showarrow=False,
        xref="paper",
        yref="paper",
        font=dict(size=24, color="black"),
        xanchor="center",
    )

    # Update layout
    fig.update_layout(
        title={
            "text": "Protein Evidence Category Transitions: ToxProt Database Evolution",
            "x": 0.5,
            "xanchor": "center",
            "font": {"size": 32},
        },
        font=dict(size=18),  # Significantly increased base font size
        height=900,  # Increased height for better spacing
        width=1400,  # Increased width for better readability with larger text
        margin=dict(t=140, b=40, l=40, r=40),  # More top margin for year labels
        showlegend=False,
    )

    # Save if output path provided
    if output_path:
        print(f"Saving plot to {output_path}")
        fig.write_html(output_path)
        # Also save as PNG if possible
        try:
            png_path = output_path.replace(".html", ".png")
            fig.write_image(png_path, width=1400, height=900, scale=2)
            print(f"Also saved PNG to {png_path}")
        except Exception as e:
            print(f"Could not save PNG: {e}")

    return fig


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(description="Create protein evidence Sankey plot")

    parser.add_argument(
        "--data-2017",
        default="data/processed/toxprot_2017.csv",
        help="Path to 2017 ToxProt data",
    )
    parser.add_argument(
        "--data-2025",
        default="data/processed/toxprot_2025.csv",
        help="Path to 2025 ToxProt data",
    )
    parser.add_argument(
        "--output",
        default="figures/protein_evidence",
        help="Output directory for the plots",
    )
    parser.add_argument("--show", action="store_true", help="Display the plot in browser")

    args = parser.parse_args()

    # Convert to absolute paths
    script_dir = Path(__file__).parent.parent
    data_2017_path = script_dir / args.data_2017
    data_2025_path = script_dir / args.data_2025
    output_dir = script_dir / args.output
    output_path = output_dir / "protein_evidence_sankey.html"

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create the plot
    fig = create_protein_evidence_sankey(str(data_2017_path), str(data_2025_path), str(output_path))

    # Show plot if requested
    if args.show:
        fig.show()

    print("Done!")


if __name__ == "__main__":
    main()
