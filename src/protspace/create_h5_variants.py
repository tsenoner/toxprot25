"""
Create H5 file variants for different protspace visualizations.

Generates four .h5 files from base embeddings:
1. All data (uses full sequences .h5)
2. Top 15 only (uses full sequences .h5, filtered by metadata)
3. Top 15 mature (uses mature sequences .h5, filtered by metadata)
4. Top 15 mature no fragments (uses mature sequences .h5, filtered by metadata)
"""

from pathlib import Path

import h5py
import pandas as pd


def filter_h5_by_metadata(
    h5_input_path: Path, metadata_csv: Path, h5_output_path: Path, variant_name: str
):
    """
    Filter H5 file to only include identifiers from metadata.

    Args:
        h5_input_path: Path to input H5 file
        metadata_csv: Path to metadata CSV with identifiers to keep
        h5_output_path: Path to output filtered H5 file
        variant_name: Name of variant for logging
    """
    if not h5_input_path.exists():
        print(f"  ✗ H5 file not found: {h5_input_path}")
        return False

    if not metadata_csv.exists():
        print(f"  ✗ Metadata file not found: {metadata_csv}")
        return False

    # Read metadata to get identifiers
    df = pd.read_csv(metadata_csv)
    identifiers_to_keep = set(df["identifier"].tolist())

    # Open input H5 and filter
    with h5py.File(h5_input_path, "r") as input_file:
        h5_keys = set(input_file.keys())
        proteins_to_keep = identifiers_to_keep.intersection(h5_keys)

        # Create output H5 with filtered proteins
        h5_output_path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(h5_output_path, "w") as output_file:
            for protein_id in proteins_to_keep:
                input_file.copy(protein_id, output_file)

    kept = len(proteins_to_keep)
    total_in_h5 = len(h5_keys)
    total_in_metadata = len(identifiers_to_keep)

    print(f"  ✓ {variant_name}: {kept}/{total_in_metadata} proteins from metadata")
    print(f"    (H5 had {total_in_h5} total, kept {kept})")

    return True


def create_h5_variants(base_dir: Path, year: str = "2025"):
    """
    Create all H5 variants for a given year.

    Args:
        base_dir: Base directory of the project
        year: Dataset year
    """
    protspace_dir = base_dir / "data" / "processed" / "protspace"

    # Base H5 files
    # Note: toxprot_2025_full.h5 contains embeddings from full sequences (with signal peptides)
    # Note: toxprot_2025_mature.h5 contains embeddings from mature sequences (signal peptides removed)
    h5_full = protspace_dir / f"toxprot_{year}_full.h5"
    h5_mature = protspace_dir / f"toxprot_{year}_mature.h5"

    print(f"Creating H5 variants for {year} dataset...")
    print("=" * 60)
    print("Base H5 files:")
    print(f"  Full sequences: {h5_full.name}")
    print(f"  Mature sequences: {h5_mature.name}")
    print()

    # Variant 1: All data (using full sequences)
    print("1. All data variant (using full sequences):")
    metadata_all = protspace_dir / f"metadata_{year}_all.csv"
    h5_all = protspace_dir / f"toxprot_{year}_all.h5"
    filter_h5_by_metadata(h5_full, metadata_all, h5_all, "All data")

    # Variant 2: Top 15 only (using full sequences)
    print("\n2. Top 15 full sequences variant:")
    metadata_top15 = protspace_dir / f"metadata_{year}_top15.csv"
    h5_top15 = protspace_dir / f"toxprot_{year}_top15.h5"
    filter_h5_by_metadata(h5_full, metadata_top15, h5_top15, "Top 15 full")

    # Variant 3: Top 15 mature (using mature sequences - signal peptides removed)
    print("\n3. Top 15 mature sequences variant (using mature sequence embeddings):")
    metadata_top15_mature = protspace_dir / f"metadata_{year}_top15_mature.csv"
    h5_top15_mature = protspace_dir / f"toxprot_{year}_top15_mature.h5"
    filter_h5_by_metadata(h5_mature, metadata_top15_mature, h5_top15_mature, "Top 15 mature")

    # Variant 4: Top 15 mature no fragments (using mature sequences, excluding fragments)
    print("\n4. Top 15 mature no fragments variant (using mature sequence embeddings):")
    metadata_top15_mature_no_frag = protspace_dir / f"metadata_{year}_top15_mature_no_fragments.csv"
    h5_top15_mature_no_frag = protspace_dir / f"toxprot_{year}_top15_mature_no_fragments.h5"
    filter_h5_by_metadata(
        h5_mature,
        metadata_top15_mature_no_frag,
        h5_top15_mature_no_frag,
        "Top 15 mature no fragments",
    )

    print("\n" + "=" * 60)
    print("✓ H5 variant creation complete!")


def main():
    """Main execution function."""
    base_dir = Path(__file__).resolve().parent.parent.parent
    create_h5_variants(base_dir, year="2025")


if __name__ == "__main__":
    main()
