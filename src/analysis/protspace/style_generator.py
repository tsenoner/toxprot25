"""Generate style.json for ProtSpace visualization.

Creates a style.json file dynamically from colors.py, ensuring color consistency
with the protein family analysis figures.
"""

import json
from pathlib import Path

import pandas as pd

from ..analyze_protein_families import get_reference_families
from ..colors import CATEGORICAL_PALETTE, NAN_COLOR, OTHER_COLOR
from .config import TOP_N


def hex_to_rgb(hex_color: str) -> str:
    """Convert hex color to rgb() format for ProtSpace.

    Args:
        hex_color: Color in hex format (e.g., "#1f77b4")

    Returns:
        Color in rgb() format (e.g., "rgb(31, 119, 180)")
    """
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgb({r}, {g}, {b})"


def hex_to_rgba(hex_color: str, alpha: float = 0.5) -> str:
    """Convert hex color to rgba() format for ProtSpace.

    Args:
        hex_color: Color in hex format (e.g., "#c7c7c7")
        alpha: Alpha transparency value (0-1)

    Returns:
        Color in rgba() format (e.g., "rgba(199, 199, 199, 0.5)")
    """
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgba({r}, {g}, {b}, {alpha})"


def generate_protein_family_colors(
    reference_families: list[str],
) -> dict[str, str]:
    """Generate color mapping for protein families.

    Colors are assigned based on the order in reference_families,
    using CATEGORICAL_PALETTE from colors.py.

    Args:
        reference_families: Ordered list of top N family names

    Returns:
        Dictionary mapping family names to rgba() colors with alpha=1.0 for opaque colors
    """
    colors = {}

    # Assign colors to reference families in order (opaque with alpha=1.0)
    for i, family in enumerate(reference_families):
        colors[family] = hex_to_rgba(CATEGORICAL_PALETTE[i % len(CATEGORICAL_PALETTE)], alpha=1.0)

    # Add special categories (semi-transparent)
    colors["Other"] = hex_to_rgba(OTHER_COLOR, 0.5)
    colors["NaN"] = hex_to_rgba(NAN_COLOR, 0.5)

    return colors


def generate_phylum_colors() -> dict[str, str]:
    """Generate color mapping for phylum categories.

    Returns:
        Dictionary mapping phylum names to rgba() colors with alpha=1.0
    """
    # Standard phylum colors (ColorBrewer Set1-like) with alpha=1.0
    return {
        "Annelida": "rgba(228, 26, 28, 1.0)",
        "Arthropoda": "rgba(55, 126, 184, 1.0)",
        "Chordata": "rgba(77, 175, 74, 1.0)",
        "Cnidaria": "rgba(152, 78, 163, 1.0)",
        "Echinodermata": "rgba(255, 127, 0, 1.0)",
        "Mollusca": "rgba(255, 255, 51, 1.0)",
        "Nematoda": "rgba(166, 86, 40, 1.0)",
        "Nemertea": "rgba(247, 129, 191, 1.0)",
        "Porifera": "rgba(153, 153, 153, 1.0)",
    }


def generate_boolean_colors() -> dict[str, dict[str, str]]:
    """Generate color mappings for boolean features.

    Returns:
        Dictionary with color mappings for has_fragment and has_signal_peptide
    """
    return {
        "has_fragment": {
            "yes": "rgba(193, 62, 62, 0.8)",
            "no": "rgba(62, 193, 62, 0.8)",
        },
        "has_signal_peptide": {
            "yes": "rgba(62, 62, 193, 0.8)",
            "no": "rgba(193, 193, 62, 0.8)",
        },
    }


def generate_style_json(
    df_reference: pd.DataFrame,
    output_path: Path,
    top_n: int = TOP_N,
    verbose: bool = True,
) -> Path:
    """Generate style.json file for ProtSpace visualization.

    Uses the reference DataFrame (typically 2025 venom_tissue) to determine
    the top N protein families, ensuring color consistency with other analyses.

    Args:
        df_reference: Reference DataFrame for determining top families
        output_path: Path to save style.json
        top_n: Number of top families to include
        verbose: Print progress messages

    Returns:
        Path to generated style.json
    """
    # Get reference families (same order as protein family analysis)
    reference_families = get_reference_families(df_reference, top_n=top_n)

    if verbose:
        print(f"Top {top_n} protein families:")
        for i, fam in enumerate(reference_families, 1):
            print(f"  {i:2d}. {fam}")

    # Generate color mappings
    protein_family_colors = generate_protein_family_colors(reference_families)
    phylum_colors = generate_phylum_colors()
    boolean_colors = generate_boolean_colors()

    # Build style.json structure
    style = {
        "Protein families": {"colors": protein_family_colors},
        "Phylum": {"colors": phylum_colors},
        **{k: {"colors": v} for k, v in boolean_colors.items()},
    }

    # Write to file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(style, f, indent=2)

    if verbose:
        print(f"\nGenerated style.json: {output_path}")

    return output_path


def get_reference_families_for_style(
    processed_csv: Path,
    top_n: int = TOP_N,
    definition: str = "venom_tissue",
) -> list[str]:
    """Load reference families from processed CSV.

    Helper function to get reference families for style generation,
    applying the standard venom_tissue filter.

    Args:
        processed_csv: Path to processed CSV file
        top_n: Number of top families
        definition: ToxProt definition filter

    Returns:
        List of top N family names
    """
    df = pd.read_csv(processed_csv)

    # Apply venom_tissue filter (matches protein family analysis default)
    if "ToxProt definition" in df.columns and definition == "venom_tissue":
        df = df[df["ToxProt definition"].isin(["venom_tissue", "both"])]

    return get_reference_families(df, top_n=top_n)
