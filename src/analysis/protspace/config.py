"""Configuration constants for ProtSpace analysis.

This module centralizes all configuration for the protspace analysis pipeline,
ensuring consistency across all commands and matching the protein family analysis.
"""

# Number of top protein families to track (matches protein family analysis)
TOP_N = 10

# Colab exchange subdirectory
COLAB_SUBDIR = "colab"

# Intermediates subdirectory (metadata CSVs and variant H5 files)
INTERMEDIATES_SUBDIR = "intermediates"

# UMAP parameters
DEFAULT_N_NEIGHBORS = 50
DEFAULT_MIN_DIST = 0.5

# Default year for analysis
DEFAULT_YEAR = "2025"

# Sequence types for FASTA/H5 base files
SEQ_TYPES = ("full", "mature", "active")

# Variant configurations
# Each variant defines:
#   - name: variant identifier
#   - description: human-readable description
#   - seq_type: sequence type ("full", "mature", or "active")
#     - full: complete precursor sequences (with signal peptides)
#     - mature: signal peptide removed
#     - active: signal peptide + propeptide removed
#   - exclude_other: exclude proteins not in top N families
#   - exclude_nan: exclude proteins without family annotation
#   - exclude_fragments: exclude fragment sequences
VARIANT_CONFIGS = {
    "full": {
        "name": "full",
        "description": "All proteins, full sequences (with signal peptides)",
        "seq_type": "full",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
    },
    "mature": {
        "name": "mature",
        "description": "All proteins, mature sequences (UniProt SP cleavage)",
        "seq_type": "mature",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
    },
    "mature_clean": {
        "name": "mature_clean",
        "description": "All proteins, mature sequences, no fragments",
        "seq_type": "mature",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": True,
    },
    "active": {
        "name": "active",
        "description": "All proteins, active sequences (SP + propeptide cleaved)",
        "seq_type": "active",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
    },
    "active_clean": {
        "name": "active_clean",
        "description": "All proteins, active sequences, no fragments",
        "seq_type": "active",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": True,
    },
}


# File naming patterns
def get_fasta_filename(year: str, seq_type: str) -> str:
    """Get FASTA filename for a given year and sequence type.

    Args:
        year: Dataset year
        seq_type: Sequence type ("full", "mature", or "active")

    Returns:
        FASTA filename string
    """
    return f"toxprot_{year}_{seq_type}.fasta"


def get_h5_base_filename(year: str, seq_type: str) -> str:
    """Get base H5 filename (from Colab embeddings).

    Args:
        year: Dataset year
        seq_type: Sequence type ("full", "mature", or "active")

    Returns:
        H5 filename string
    """
    return f"toxprot_{year}_{seq_type}.h5"


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
