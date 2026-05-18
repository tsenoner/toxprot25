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

# Curated demo bundles subdirectory
DEMO_SUBDIR = "demo"

# UMAP parameters
DEFAULT_N_NEIGHBORS = 50
DEFAULT_MIN_DIST = 0.5

# Default year for analysis
DEFAULT_YEAR = "2025"

# Embedder model name embedded in H5 input via `:model` syntax in protspace v4
PROTSPACE_MODEL_NAME = "prot_t5"

# Sequence types for FASTA/H5 base files
SEQ_TYPES = ("full", "mature", "active")

# Annotation columns shipped in the curated demo bundles, in display order.
# Keys are CSV-source column names, values are the canonical snake_case
# annotation names recognised by the ProtSpace web UI grouping table
# (packages/core/src/components/control-bar/annotation-categories.ts).
# Matching names get auto-grouped under UniProt / Taxonomy; unmatched names
# fall into "Other".
DEMO_ANNOTATION_COLUMNS: dict[str, str] = {
    "identifier": "identifier",
    "Protein families": "protein_families",
    "Length": "length",
    "Phylum": "phylum",
    "Class": "class",
    "Order": "order",
    "Family": "family",
    "Genus": "genus",
    "Species": "species",
    "Habitat": "habitat",
    "Habitat_Detailed": "habitat_detailed",
    "Source tissues": "source_tissues",
    "PTM Keywords": "ptm_keywords",
    "PTM Summary": "ptm_summary",
    "Protein existence": "protein_existence",
    "has_propeptide": "has_propeptide",
    "has_fragment": "has_fragment",
}

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
#   - bundle_kind: "analysis" (top-N + Other, lean schema, manuscript palette)
#                  or "demo" (full family names, expanded schema, no styling)
VARIANT_CONFIGS = {
    "full": {
        "name": "full",
        "description": "All proteins, full sequences (with signal peptides)",
        "seq_type": "full",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
        "bundle_kind": "analysis",
    },
    "mature": {
        "name": "mature",
        "description": "All proteins, mature sequences (UniProt SP cleavage)",
        "seq_type": "mature",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
        "bundle_kind": "analysis",
    },
    "mature_clean": {
        "name": "mature_clean",
        "description": "All proteins, mature sequences, no fragments",
        "seq_type": "mature",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": True,
        "bundle_kind": "analysis",
    },
    "active": {
        "name": "active",
        "description": "All proteins, active sequences (SP + propeptide cleaved)",
        "seq_type": "active",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
        "bundle_kind": "analysis",
    },
    "active_clean": {
        "name": "active_clean",
        "description": "All proteins, active sequences, no fragments",
        "seq_type": "active",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": True,
        "bundle_kind": "analysis",
    },
    "demo_mature": {
        "name": "demo_mature",
        "description": "Demo bundle: mature sequences (SP removed), fragments included",
        "seq_type": "mature",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": False,
        "bundle_kind": "demo",
    },
    "demo_mature_clean": {
        "name": "demo_mature_clean",
        "description": "Demo bundle: mature sequences (SP removed), no fragments",
        "seq_type": "mature",
        "exclude_other": False,
        "exclude_nan": False,
        "exclude_fragments": True,
        "bundle_kind": "demo",
    },
}


ANALYSIS_VARIANTS = tuple(
    name for name, cfg in VARIANT_CONFIGS.items() if cfg["bundle_kind"] == "analysis"
)
DEMO_VARIANTS = tuple(
    name for name, cfg in VARIANT_CONFIGS.items() if cfg["bundle_kind"] == "demo"
)


# File naming patterns
def get_fasta_filename(year: str, seq_type: str) -> str:
    """Get FASTA filename for a given year and sequence type."""
    return f"toxprot_{year}_{seq_type}.fasta"


def get_h5_base_filename(year: str, seq_type: str) -> str:
    """Get base H5 filename (from Colab embeddings)."""
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
    """Get protspace final parquetbundle filename.

    For analysis variants this is the styled bundle; for demo variants there
    is no separate styled file (the demo bundles ship unstyled).
    """
    return f"protspace_{year}_{variant}.parquetbundle"
