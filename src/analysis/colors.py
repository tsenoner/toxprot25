"""
Centralized color configuration for ToxProt figures.

This module provides consistent colors across all analysis figures
to ensure visual coherence in the manuscript.
"""

# =============================================================================
# Year Colors - For temporal comparisons (2005, 2015, 2025)
# =============================================================================
# Sunset gradient (bright gold → tangerine → maroon)
YEAR_COLORS = {
    2005: "#ffb703",  # Bright gold - earliest
    2015: "#fb5607",  # Tangerine - middle
    2025: "#800020",  # Maroon - latest
}


# =============================================================================
# ToxProt Definition Colors
# =============================================================================
DEFINITION_COLORS = {
    "venom_tissue": "#4C78A8",  # Blue - venom tissue annotation
    "kw_toxin": "#F58518",  # Orange - toxin keyword
    "both": "#54A24B",  # Green - both criteria
}

# Shorthand aliases
VENOM_COLOR = DEFINITION_COLORS["venom_tissue"]
KEYWORD_COLOR = DEFINITION_COLORS["kw_toxin"]
BOTH_COLOR = DEFINITION_COLORS["both"]


# =============================================================================
# Taxa Order Colors - For top 5 orders + Others
# =============================================================================
TAXA_ORDER_COLORS = {
    "Squamata": "#1f77b4",  # Blue - Snakes & Lizards
    "Araneae": "#ff7f0e",  # Orange - Spiders
    "Neogastropoda": "#2ca02c",  # Green - Cone Snails
    "Scorpiones": "#d62728",  # Red - Scorpions
    "Hymenoptera": "#9467bd",  # Purple - Bees & Wasps
    "Others": "#c7c7c7",  # Light gray - Other orders
}

# List version for stackplots (maintains order)
TAXA_ORDER_COLORS_LIST = [
    "#1f77b4",  # Squamata
    "#ff7f0e",  # Araneae
    "#2ca02c",  # Neogastropoda
    "#d62728",  # Scorpiones
    "#9467bd",  # Hymenoptera
    "#c7c7c7",  # Others
]


# =============================================================================
# Habitat Colors
# =============================================================================
HABITAT_COLORS = {
    "terrestrial": "#2e7d32",  # Forest green
    "marine": "#1565c0",  # Ocean blue
}

# Detailed habitat colors (for freshwater/estuarine if present)
HABITAT_COLORS_DETAILED = {
    "terrestrial": "#2e7d32",  # Forest green
    "marine": "#1565c0",  # Ocean blue
    "freshwater": "#00838f",  # Teal
    "estuarine": "#5e35b1",  # Purple
}


# =============================================================================
# Protein Evidence Colors
# =============================================================================
PROTEIN_EVIDENCE_COLORS = {
    "PE1": "#2ecc71",  # Green - Protein level
    "PE2": "#3498db",  # Blue - Transcript level
    "PE3": "#9b59b6",  # Purple - Homology
    "PE4": "#f39c12",  # Yellow - Predicted
    "PE5": "#e74c3c",  # Red - Uncertain
    "Removed": "#808080",  # Gray
    "New": "#A0A0A0",  # Light gray
}


# =============================================================================
# General Palette for categorical data (e.g., protein families)
# =============================================================================
# tab20 without gray tones (indices 14, 15)
CATEGORICAL_PALETTE = [
    "#1f77b4",
    "#aec7e8",
    "#ff7f0e",
    "#ffbb78",
    "#2ca02c",
    "#98df8a",
    "#d62728",
    "#ff9896",
    "#9467bd",
    "#c5b0d5",
    "#8c564b",
    "#c49c94",
    "#e377c2",
    "#f7b6d2",
    "#bcbd22",
    "#dbdb8d",
    "#17becf",
    "#9edae5",
]


def get_categorical_color(index: int) -> str:
    """Get color from categorical palette by index (cycles if needed)."""
    return CATEGORICAL_PALETTE[index % len(CATEGORICAL_PALETTE)]


# =============================================================================
# Special colors
# =============================================================================
OTHER_COLOR = "#c7c7c7"  # Light gray for "Other" categories
NAN_COLOR = "#e0e0e0"  # Lighter gray for NaN/missing
HIGHLIGHT_COLOR = "#e74c3c"  # Red for highlighting
