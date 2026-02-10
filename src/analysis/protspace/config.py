"""Configuration constants for ProtSpace analysis.

This module centralizes all configuration for the protspace analysis pipeline,
ensuring consistency across all commands and matching the protein family analysis.
"""

# Number of top protein families to track (matches protein family analysis)
TOP_N = 10

# Colab exchange subdirectory
COLAB_SUBDIR = "colab"

# UMAP parameters
DEFAULT_N_NEIGHBORS = 50
DEFAULT_MIN_DIST = 0.5

# Default year for analysis
DEFAULT_YEAR = "2025"

# Variant configurations
# Each variant defines:
#   - name: variant identifier
#   - description: human-readable description
#   - uses_mature: whether to use mature (SP-cleaved) sequences
#   - exclude_other: exclude proteins not in top N families
#   - exclude_nan: exclude proteins without family annotation
#   - exclude_fragments: exclude fragment sequences
VARIANT_CONFIGS = {
    "full": {
        "name": "full",
        "description": "All proteins, full sequences (with signal peptides)",
        "uses_mature": False,
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
    },
    "mature": {
        "name": "mature",
        "description": "All proteins, mature sequences (UniProt SP cleavage)",
        "uses_mature": True,
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
    },
    "mature_clean": {
        "name": "mature_clean",
        "description": "All proteins, mature sequences, no fragments",
        "uses_mature": True,
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": True,
    },
}


# File naming patterns
def get_fasta_filename(year: str, is_mature: bool) -> str:
    """Get FASTA filename for a given year and sequence type.

    Args:
        year: Dataset year
        is_mature: Whether this is mature sequences (signal peptide removed)

    Returns:
        FASTA filename string
    """
    if not is_mature:
        return f"toxprot_{year}_full.fasta"
    return f"toxprot_{year}_mature.fasta"


def get_h5_base_filename(year: str, is_mature: bool) -> str:
    """Get base H5 filename (from Colab embeddings).

    Args:
        year: Dataset year
        is_mature: Whether this is mature sequences

    Returns:
        H5 filename string
    """
    if not is_mature:
        return f"toxprot_{year}_full.h5"
    return f"toxprot_{year}_mature.h5"


def get_h5_variant_filename(year: str, variant: str) -> str:
    """Get H5 filename for a specific variant."""
    return f"toxprot_{year}_{variant}.h5"


def get_metadata_filename(year: str, variant: str) -> str:
    """Get metadata CSV filename for a specific variant."""
    return f"metadata_{year}_{variant}.csv"


def get_protspace_output_filename(year: str, variant: str) -> str:
    """Get protspace intermediate (unstyled) parquetbundle filename."""
    return f"protspace_{year}_{variant}_tmp.parquetbundle"


def get_protspace_styled_filename(year: str, variant: str) -> str:
    """Get protspace final (styled) parquetbundle filename."""
    return f"protspace_{year}_{variant}.parquetbundle"


def get_plot_filename(year: str, variant: str) -> str:
    """Get plot filename for a specific variant."""
    return f"{year}_{variant}.png"
