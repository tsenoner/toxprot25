#!/usr/bin/env python3
"""
Generate FASTA files for ProtT5 embedding generation.

Creates two FASTA files:
1. Full sequences (with signal peptides) - for variants 1 & 2
2. Mature sequences (signal peptides removed) - for variants 3 & 4
"""

from pathlib import Path

import pandas as pd


def remove_signal_peptide(sequence: str, signal_range: str) -> str:
    """
    Remove signal peptide from sequence based on range.

    Args:
        sequence: Full protein sequence
        signal_range: Range string like "1-22" indicating signal peptide positions

    Returns:
        Mature sequence with signal peptide removed
    """
    if not isinstance(signal_range, str) or "-" not in signal_range:
        # No valid signal peptide range, return original sequence
        return sequence

    try:
        # Signal peptide range is 1-based, e.g., "1-22" means positions 1 to 22
        end_pos = int(signal_range.split("-")[1])

        if end_pos > 0 and end_pos < len(sequence):
            # Return sequence after signal peptide (0-based indexing)
            return sequence[end_pos:]
        else:
            return sequence
    except (ValueError, IndexError):
        # If parsing fails, return original sequence
        return sequence


def generate_fasta_files(interim_tsv: Path, output_dir: Path, year: str = "2025"):
    """
    Generate FASTA files for embedding generation.

    Args:
        interim_tsv: Path to interim TSV with sequences and signal peptide info
        output_dir: Directory to save FASTA files
        year: Dataset year
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read interim TSV
    print(f"Reading {interim_tsv}...")
    df = pd.read_csv(
        interim_tsv, sep="\t", usecols=["Entry", "Sequence", "Signal peptide (range)"]
    )
    print(f"Loaded {len(df)} entries")

    # Output files
    fasta_full = output_dir / f"toxprot_{year}_full.fasta"
    fasta_mature = output_dir / f"toxprot_{year}_mature.fasta"

    # Generate full sequences FASTA
    print(f"\nGenerating full sequences FASTA: {fasta_full}")
    with open(fasta_full, "w") as f:
        for _, row in df.iterrows():
            entry = row["Entry"]
            sequence = row["Sequence"]

            if isinstance(sequence, str) and sequence:
                f.write(f">{entry}\n")
                f.write(f"{sequence}\n")

    print(f"  ✓ Wrote {len(df)} sequences")

    # Generate mature sequences FASTA (with signal peptides removed)
    print(f"\nGenerating mature sequences FASTA: {fasta_mature}")
    mature_count = 0
    unchanged_count = 0

    with open(fasta_mature, "w") as f:
        for _, row in df.iterrows():
            entry = row["Entry"]
            sequence = row["Sequence"]
            signal_range = row.get("Signal peptide (range)")

            if isinstance(sequence, str) and sequence:
                # Remove signal peptide if present
                mature_sequence = remove_signal_peptide(sequence, signal_range)

                if mature_sequence != sequence:
                    mature_count += 1
                else:
                    unchanged_count += 1

                f.write(f">{entry}\n")
                f.write(f"{mature_sequence}\n")

    print(f"  ✓ Wrote {len(df)} sequences")
    print(f"    - {mature_count} sequences with signal peptide removed")
    print(f"    - {unchanged_count} sequences unchanged (no signal peptide)")

    # Summary
    print("\n" + "=" * 60)
    print("Summary:")
    print("=" * 60)
    print(f"Full sequences FASTA: {fasta_full}")
    print(f"Mature sequences FASTA: {fasta_mature}")
    print("\nNext steps:")
    print(f"1. Generate ProtT5 embeddings for {fasta_full.name}")
    print(f"   → Save as toxprot_{year}_full.h5")
    print(f"2. Generate ProtT5 embeddings for {fasta_mature.name}")
    print(f"   → Save as toxprot_{year}_mature.h5")


def main():
    """Main execution function."""
    base_dir = Path(__file__).resolve().parent.parent.parent
    interim_tsv = base_dir / "data" / "interim" / "toxprot_2025.tsv"
    output_dir = base_dir / "data" / "processed" / "protspace"

    generate_fasta_files(interim_tsv, output_dir, year="2025")

    print("\n✓ FASTA generation complete!")


if __name__ == "__main__":
    main()
