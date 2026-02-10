"""Generate FASTA files for ProtT5 embedding generation.

Creates FASTA files:
1. Full sequences (with signal peptides) - for all and top10_full variants
2. Mature sequences (signal peptides removed based on UniProt) - for top10_mature variants
"""

from pathlib import Path

import pandas as pd

from .config import COLAB_SUBDIR, get_fasta_filename


def remove_signal_peptide(sequence: str, signal_range: str | None) -> str:
    """Remove signal peptide from sequence based on UniProt annotation range.

    Args:
        sequence: Full protein sequence
        signal_range: Range string like "1-22" indicating signal peptide positions
            (1-based, inclusive). None or invalid format returns original sequence.

    Returns:
        Mature sequence with signal peptide removed, or original if no valid range.
    """
    if not isinstance(signal_range, str) or "-" not in signal_range:
        return sequence

    try:
        # Signal peptide range is 1-based, e.g., "1-22" means positions 1 to 22
        end_pos = int(signal_range.split("-")[1])

        if 0 < end_pos < len(sequence):
            # Return sequence after signal peptide (0-based indexing)
            return sequence[end_pos:]
        return sequence
    except (ValueError, IndexError):
        return sequence


def generate_fasta_files(
    interim_tsv: Path,
    output_dir: Path,
    year: str = "2025",
    processed_csv: Path | None = None,
    definition: str = "venom_tissue",
    verbose: bool = True,
) -> tuple[Path, Path]:
    """Generate FASTA files for embedding generation.

    Args:
        interim_tsv: Path to interim TSV with sequences and signal peptide info
        output_dir: Base protspace directory (FASTA files go to colab/ subdirectory)
        year: Dataset year
        processed_csv: Path to processed CSV with ToxProt definition column.
            If provided, filters entries by definition.
        definition: Definition filter to apply ("venom_tissue" filters to
            venom_tissue + both). Only used if processed_csv is provided.
        verbose: Print progress messages

    Returns:
        Tuple of (full_fasta_path, mature_fasta_path)
    """
    # FASTA files go to colab/ subdirectory for exchange with Google Colab
    colab_dir = output_dir / COLAB_SUBDIR
    colab_dir.mkdir(parents=True, exist_ok=True)

    # Load valid entries from processed CSV if provided
    valid_entries: set[str] | None = None
    if processed_csv is not None and processed_csv.exists():
        if verbose:
            print(f"Loading definition filter from {processed_csv}...")
        df_processed = pd.read_csv(
            processed_csv, usecols=["Entry", "ToxProt definition"]
        )
        if "ToxProt definition" in df_processed.columns and definition == "venom_tissue":
            df_processed = df_processed[
                df_processed["ToxProt definition"].isin(["venom_tissue", "both"])
            ]
        valid_entries = set(df_processed["Entry"])
        if verbose:
            print(f"  Filtered to {len(valid_entries)} entries with definition: {definition}")

    # Read interim TSV
    if verbose:
        print(f"Reading {interim_tsv}...")

    df = pd.read_csv(
        interim_tsv,
        sep="\t",
        usecols=["Entry", "Sequence", "Signal peptide (range)"],
    )

    total_entries = len(df)

    # Filter to valid entries if specified
    if valid_entries is not None:
        df = df[df["Entry"].isin(valid_entries)]

    if verbose:
        if valid_entries is not None:
            print(f"Loaded {len(df)} entries (filtered from {total_entries})")
        else:
            print(f"Loaded {len(df)} entries")

    # Output files (in colab/ subdirectory)
    fasta_full = colab_dir / get_fasta_filename(year, is_mature=False)
    fasta_mature = colab_dir / get_fasta_filename(year, is_mature=True)

    # Generate full sequences FASTA
    if verbose:
        print(f"\nGenerating full sequences FASTA: {fasta_full}")

    full_count = 0
    with open(fasta_full, "w") as f:
        for _, row in df.iterrows():
            entry = row["Entry"]
            sequence = row["Sequence"]

            if isinstance(sequence, str) and sequence:
                f.write(f">{entry}\n")
                f.write(f"{sequence}\n")
                full_count += 1

    if verbose:
        print(f"  Wrote {full_count} sequences")

    # Generate mature sequences FASTA (with signal peptides removed)
    if verbose:
        print(f"\nGenerating mature sequences FASTA: {fasta_mature}")

    mature_modified = 0
    mature_unchanged = 0

    with open(fasta_mature, "w") as f:
        for _, row in df.iterrows():
            entry = row["Entry"]
            sequence = row["Sequence"]
            signal_range = row.get("Signal peptide (range)")

            if isinstance(sequence, str) and sequence:
                mature_sequence = remove_signal_peptide(sequence, signal_range)

                if mature_sequence != sequence:
                    mature_modified += 1
                else:
                    mature_unchanged += 1

                f.write(f">{entry}\n")
                f.write(f"{mature_sequence}\n")

    if verbose:
        print(f"  Wrote {mature_modified + mature_unchanged} sequences")
        print(f"    - {mature_modified} with signal peptide removed")
        print(f"    - {mature_unchanged} unchanged (no signal peptide annotation)")

    return fasta_full, fasta_mature


def print_next_steps(
    fasta_full: Path,
    fasta_mature: Path,
    year: str = "2025",
) -> None:
    """Print instructions for generating embeddings via Colab.

    Args:
        fasta_full: Path to full sequences FASTA
        fasta_mature: Path to mature sequences FASTA (UniProt-based)
        year: Dataset year
    """
    colab_dir = fasta_full.parent
    print("\n" + "=" * 60)
    print("Next Steps: Generate embeddings via Google Colab")
    print("=" * 60)
    print(f"\n1. Upload FASTA files from {colab_dir}/ to Google Colab:")
    print(f"   - {fasta_full.name}")
    print(f"   - {fasta_mature.name}")
    print(
        "\n2. Open: https://colab.research.google.com/github/tsenoner/protspace/"
        "blob/master/colab/ProtSpace_Embeddings.ipynb"
    )
    print(f"\n3. Download H5 files to {colab_dir}/:")
    print(f"   - toxprot_{year}_full.h5")
    print(f"   - toxprot_{year}_mature.h5")
    print("\n4. Run: toxprot analysis protspace prepare")
