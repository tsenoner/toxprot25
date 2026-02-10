"""Prepare metadata and H5 file variants for ProtSpace visualization.

This module handles:
1. Creating metadata CSV files for each variant (all, top10_full, top10_mature, top10_mature_clean)
2. Filtering H5 embedding files to match metadata variants
"""

from pathlib import Path

import h5py
import pandas as pd

from ..analyze_protein_families import get_reference_families, normalize_family_name
from .config import (
    COLAB_SUBDIR,
    TOP_N,
    VARIANT_CONFIGS,
    get_h5_base_filename,
    get_h5_variant_filename,
    get_metadata_filename,
)


def process_protein_families(
    df: pd.DataFrame,
    reference_families: list[str],
    column: str = "Protein families",
) -> pd.DataFrame:
    """Process protein families column to use top N + Other + NaN.

    Args:
        df: Input DataFrame
        reference_families: List of top N family names to keep
        column: Column name containing family information

    Returns:
        DataFrame with processed Protein families column
    """
    df = df.copy()

    def categorize(value):
        if pd.isna(value):
            return value  # Keep NaN as is
        normalized = normalize_family_name(value)
        if normalized in reference_families:
            return normalized
        return "Other"

    df[column] = df[column].apply(categorize)
    return df


def create_metadata_csv(
    df: pd.DataFrame,
    variant_config: dict,
    reference_families: list[str],
    output_path: Path,
    verbose: bool = True,
) -> int:
    """Create metadata CSV for a specific variant.

    Args:
        df: Base DataFrame with all data
        variant_config: Configuration for this variant
        reference_families: List of top N family names
        output_path: Path to save metadata CSV
        verbose: Print progress messages

    Returns:
        Number of entries in the metadata file
    """
    df_variant = df.copy()

    # Apply filters based on variant configuration
    if variant_config["exclude_nan"]:
        df_variant = df_variant[df_variant["Protein families"].notna()]

    if variant_config["exclude_other"]:
        df_variant = df_variant[df_variant["Protein families"] != "Other"]

    if variant_config["exclude_fragments"]:
        df_variant = df_variant[df_variant["has_fragment"] == "no"]

    # Select columns for output
    columns = [
        "identifier",
        "Protein families",
        "Phylum",
        "has_fragment",
        "has_signal_peptide",
    ]

    # Remove has_fragment column if fragments are excluded (redundant)
    if variant_config["exclude_fragments"]:
        columns.remove("has_fragment")

    # Remove has_signal_peptide for mature variants (signal peptide already cleaved)
    if variant_config["uses_mature"] and "has_signal_peptide" in columns:
        columns.remove("has_signal_peptide")

    df_variant[columns].to_csv(output_path, index=False)

    if verbose:
        print(f"  {variant_config['name']}: {len(df_variant)} entries -> {output_path.name}")

    return len(df_variant)


def filter_h5_by_metadata(
    h5_input: Path,
    metadata_csv: Path,
    h5_output: Path,
    verbose: bool = True,
) -> tuple[int, int]:
    """Filter H5 file to include only identifiers present in metadata.

    Args:
        h5_input: Path to input H5 file (base embeddings)
        metadata_csv: Path to metadata CSV with identifiers to keep
        h5_output: Path to output filtered H5 file
        verbose: Print progress messages

    Returns:
        Tuple of (entries_kept, total_in_h5)
    """
    if not h5_input.exists():
        raise FileNotFoundError(f"H5 file not found: {h5_input}")

    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_csv}")

    # Read identifiers from metadata
    df = pd.read_csv(metadata_csv)
    identifiers_to_keep = set(df["identifier"].tolist())

    # Filter H5 file
    with h5py.File(h5_input, "r") as input_file:
        h5_keys = set(input_file.keys())
        proteins_to_keep = identifiers_to_keep.intersection(h5_keys)

        h5_output.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(h5_output, "w") as output_file:
            for protein_id in proteins_to_keep:
                input_file.copy(protein_id, output_file)

    if verbose:
        print(f"    H5: {len(proteins_to_keep)}/{len(h5_keys)} embeddings -> {h5_output.name}")

    return len(proteins_to_keep), len(h5_keys)


def prepare_all_variants(
    processed_csv: Path,
    interim_tsv: Path,
    protspace_dir: Path,
    year: str = "2025",
    top_n: int = TOP_N,
    definition: str = "venom_tissue",
    verbose: bool = True,
) -> dict[str, dict]:
    """Prepare metadata and H5 files for all variants.

    Args:
        processed_csv: Path to processed CSV (e.g., toxprot_2025.csv)
        interim_tsv: Path to interim TSV with signal peptide info
        protspace_dir: Directory for protspace files
        year: Dataset year
        top_n: Number of top families to track
        definition: ToxProt definition filter
        verbose: Print progress messages

    Returns:
        Dictionary with results for each variant
    """
    protspace_dir.mkdir(parents=True, exist_ok=True)

    # Load processed CSV
    if verbose:
        print(f"Loading {processed_csv}...")
    df = pd.read_csv(processed_csv)

    # Apply ToxProt definition filter
    if "ToxProt definition" in df.columns and definition == "venom_tissue":
        df = df[df["ToxProt definition"].isin(["venom_tissue", "both"])]

    if verbose:
        print(f"Loaded {len(df)} entries (definition: {definition})")

    # Get reference families (same order as protein family analysis)
    reference_families = get_reference_families(df, top_n=top_n)

    if verbose:
        print(f"\nTop {top_n} protein families:")
        for i, fam in enumerate(reference_families, 1):
            print(f"  {i:2d}. {fam}")

    # Load signal peptide information from interim TSV
    if verbose:
        print(f"\nLoading signal peptide info from {interim_tsv}...")
    df_interim = pd.read_csv(
        interim_tsv,
        sep="\t",
        usecols=["Entry", "Signal peptide (range)"],
    )

    # Merge signal peptide information
    df = df.merge(df_interim, on="Entry", how="left")

    # Rename Entry to identifier (required by protspace)
    df = df.rename(columns={"Entry": "identifier"})

    # Add helper columns
    df["has_signal_peptide"] = df["Signal peptide (range)"].notna().map(
        {True: "yes", False: "no"}
    )
    df["has_fragment"] = (df["Fragment"].astype(str) == "fragment").map(
        {True: "yes", False: "no"}
    )

    # Process protein families to top N + Other
    df = process_protein_families(df, reference_families)

    # Check for base H5 files (in colab/ subdirectory)
    colab_dir = protspace_dir / COLAB_SUBDIR
    h5_full = colab_dir / get_h5_base_filename(year, is_mature=False)
    h5_mature = colab_dir / get_h5_base_filename(year, is_mature=True)

    # Map sequence type to H5 file
    h5_files = {
        "full": h5_full,
        "mature": h5_mature,
    }

    # Check which H5 files are needed based on variant configs
    h5_missing = []
    if not h5_full.exists():
        h5_missing.append(h5_full)
    if not h5_mature.exists():
        h5_missing.append(h5_mature)

    if h5_missing:
        _print_h5_missing_error(h5_missing, year)
        raise FileNotFoundError("Required H5 embedding files not found")

    # Create variants
    if verbose:
        print("\nCreating metadata and H5 variants...")

    results = {}
    for variant_name, config in VARIANT_CONFIGS.items():
        if verbose:
            print(f"\n{config['description']}:")

        # Determine which base H5 to use
        if config["uses_mature"]:
            h5_key = "mature"
        else:
            h5_key = "full"

        h5_base = h5_files.get(h5_key)

        # Skip variant if H5 not available
        if h5_base is None or not h5_base.exists():
            if verbose:
                print(f"  Skipping: H5 file not found ({h5_base})")
            continue

        # Create metadata
        metadata_path = protspace_dir / get_metadata_filename(year, variant_name)
        n_entries = create_metadata_csv(
            df, config, reference_families, metadata_path, verbose=verbose
        )

        # Filter H5
        h5_variant = protspace_dir / get_h5_variant_filename(year, variant_name)
        n_kept, n_total = filter_h5_by_metadata(
            h5_base, metadata_path, h5_variant, verbose=verbose
        )

        results[variant_name] = {
            "metadata": metadata_path,
            "h5": h5_variant,
            "n_entries": n_entries,
            "n_embeddings": n_kept,
        }

    # Print summary
    if verbose:
        print("\n" + "=" * 60)
        print("Summary:")
        print("=" * 60)
        for variant_name, info in results.items():
            print(f"  {variant_name}: {info['n_entries']} entries, {info['n_embeddings']} embeddings")

    return results


def _print_h5_missing_error(missing_files: list[Path], year: str) -> None:
    """Print helpful error message when H5 files are missing."""
    print("\n" + "=" * 60)
    print("ERROR: Embedding files not found:")
    print("=" * 60)
    for f in missing_files:
        print(f"  - {f}")
    print("\nTo generate embeddings:")
    print("1. Run: toxprot analysis protspace generate-fasta")
    print("2. Upload FASTA files to Google Colab")
    print(
        "3. Open: https://colab.research.google.com/github/tsenoner/protspace/"
        "blob/master/colab/ProtSpace_Embeddings.ipynb"
    )
    print(f"4. Download H5 files to {missing_files[0].parent}/")
    print("5. Re-run this command")
