"""Generate FASTA files for ProtT5 embedding generation.

Creates FASTA files:
1. Full sequences (with signal peptides) - for full variants
2. Mature sequences (signal peptides removed) - for mature variants
3. Active sequences (signal peptides + propeptides removed) - for active variants
"""

from pathlib import Path

import pandas as pd

from .config import COLAB_SUBDIR, get_fasta_filename


def _parse_ranges(range_str: str | None) -> list[tuple[int, int]]:
    """Parse semicolon-separated ranges into list of (start, end) tuples.

    Args:
        range_str: Range string like "1-22" or "23-45; 67-89" (1-based inclusive).
            None or invalid format returns empty list.

    Returns:
        List of (start, end) tuples (1-based inclusive).
    """
    if not isinstance(range_str, str) or not range_str:
        return []
    ranges = []
    for part in range_str.split(";"):
        part = part.strip()
        if "-" not in part:
            continue
        try:
            start_str, end_str = part.split("-")
            ranges.append((int(start_str), int(end_str)))
        except (ValueError, IndexError):
            continue
    return ranges


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


def remove_regions(sequence: str, ranges: list[tuple[int, int]]) -> str:
    """Remove specified 1-based inclusive ranges from sequence.

    Args:
        sequence: Full protein sequence
        ranges: List of (start, end) tuples (1-based inclusive) to remove.

    Returns:
        Sequence with specified regions removed. Returns original if no
        valid ranges or result would be empty.
    """
    if not ranges:
        return sequence

    keep = [True] * len(sequence)
    for start, end in ranges:
        # Convert 1-based inclusive to 0-based
        for i in range(start - 1, min(end, len(sequence))):
            keep[i] = False

    result = "".join(c for c, k in zip(sequence, keep) if k)
    return result if result else sequence


def _extract_regions(sequence: str, ranges: list[tuple[int, int]]) -> str:
    """Extract specified 1-based inclusive ranges from sequence and concatenate.

    Args:
        sequence: Full protein sequence
        ranges: List of (start, end) tuples (1-based inclusive) to extract.

    Returns:
        Concatenated subsequences from the specified ranges.
        Returns original sequence if no valid ranges or result would be empty.
    """
    if not ranges:
        return sequence

    parts = []
    for start, end in sorted(ranges):
        # Convert 1-based inclusive to 0-based
        s = max(start - 1, 0)
        e = min(end, len(sequence))
        if s < e:
            parts.append(sequence[s:e])

    result = "".join(parts)
    return result if result else sequence


def _resolve_chain_ranges(ranges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Resolve overlapping chain ranges by removing parent ranges that contain sub-chains.

    UniProt sometimes lists a parent chain (e.g., 23-1642) alongside its processed
    sub-chains (e.g., 23-649, 733-984, 1264-1642). We keep only the most specific
    (non-containing) ranges.

    If no containment relationships exist, returns all ranges unchanged.
    """
    if len(ranges) <= 1:
        return ranges

    # Check if any range contains another
    resolved = []
    for i, (s1, e1) in enumerate(ranges):
        is_parent = False
        for j, (s2, e2) in enumerate(ranges):
            if i != j and s1 <= s2 and e1 >= e2 and (s1 != s2 or e1 != e2):
                # range i contains range j — i is a parent
                is_parent = True
                break
        if not is_parent:
            resolved.append((s1, e1))

    return resolved if resolved else ranges


def remove_signal_and_propeptide(
    sequence: str,
    signal_range: str | None,
    propeptide_range: str | None,
    chain_range: str | None = None,
    peptide_range: str | None = None,
) -> str:
    """Get active peptide sequence using the best available annotations.

    Priority:
    1. PEPTIDE annotations (most specific active peptide products)
    2. CHAIN annotations (mature protein chains, with parent chains removed)
    3. Subtraction: remove SP + propeptide positions from full sequence
    4. Remove SP only
    5. Return original sequence

    Args:
        sequence: Full protein sequence
        signal_range: Signal peptide range string (e.g., "1-22")
        propeptide_range: Propeptide range string (e.g., "23-45; 67-89")
        chain_range: Chain range string (e.g., "46-66; 90-120")
        peptide_range: Peptide range string (e.g., "51-71")

    Returns:
        Active sequence based on best available annotation.
    """
    # Priority 1: Use PEPTIDE annotations if present
    peptide_ranges = _parse_ranges(peptide_range)
    if peptide_ranges:
        return _extract_regions(sequence, peptide_ranges)

    # Priority 2: Use CHAIN annotations if present
    chain_ranges = _parse_ranges(chain_range)
    if chain_ranges:
        resolved = _resolve_chain_ranges(chain_ranges)
        return _extract_regions(sequence, resolved)

    # Priority 3: Subtraction — remove SP + propeptide positions
    all_ranges = _parse_ranges(signal_range) + _parse_ranges(propeptide_range)
    if all_ranges:
        return remove_regions(sequence, all_ranges)

    return sequence


def generate_fasta_files(
    interim_tsv: Path,
    output_dir: Path,
    year: str = "2025",
    processed_csv: Path | None = None,
    definition: str = "venom_tissue",
    verbose: bool = True,
) -> tuple[Path, Path, Path]:
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
        Tuple of (full_fasta_path, mature_fasta_path, active_fasta_path)
    """
    # FASTA files go to colab/ subdirectory for exchange with Google Colab
    colab_dir = output_dir / COLAB_SUBDIR
    colab_dir.mkdir(parents=True, exist_ok=True)

    # Load valid entries from processed CSV if provided
    valid_entries: set[str] | None = None
    if processed_csv is not None and processed_csv.exists():
        if verbose:
            print(f"Loading definition filter from {processed_csv}...")
        df_processed = pd.read_csv(processed_csv, usecols=["Entry", "ToxProt definition"])
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

    # Determine available columns
    with open(interim_tsv) as f:
        header = f.readline().strip().split("\t")

    usecols = ["Entry", "Sequence", "Signal peptide (range)"]
    has_propeptide_range = "Propeptide (range)" in header
    if has_propeptide_range:
        usecols.append("Propeptide (range)")
    has_chain_range = "Chain (range)" in header
    if has_chain_range:
        usecols.append("Chain (range)")
    has_peptide_range = "Peptide (range)" in header
    if has_peptide_range:
        usecols.append("Peptide (range)")

    df = pd.read_csv(interim_tsv, sep="\t", usecols=usecols)

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
    fasta_full = colab_dir / get_fasta_filename(year, seq_type="full")
    fasta_mature = colab_dir / get_fasta_filename(year, seq_type="mature")
    fasta_active = colab_dir / get_fasta_filename(year, seq_type="active")

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

    # Generate active sequences FASTA (with signal peptides + propeptides removed)
    if verbose:
        print(f"\nGenerating active sequences FASTA: {fasta_active}")

    active_sp_only = 0
    active_pp_only = 0
    active_both = 0
    active_unchanged = 0

    active_used_peptide = 0
    active_used_chain = 0

    with open(fasta_active, "w") as f:
        for _, row in df.iterrows():
            entry = row["Entry"]
            sequence = row["Sequence"]
            signal_range = row.get("Signal peptide (range)")
            propeptide_range = row.get("Propeptide (range)") if has_propeptide_range else None
            chain_range = row.get("Chain (range)") if has_chain_range else None
            peptide_range = row.get("Peptide (range)") if has_peptide_range else None

            if isinstance(sequence, str) and sequence:
                active_sequence = remove_signal_and_propeptide(
                    sequence, signal_range, propeptide_range,
                    chain_range, peptide_range,
                )

                has_sp = bool(_parse_ranges(signal_range))
                has_pp = bool(_parse_ranges(propeptide_range))
                has_pep = bool(_parse_ranges(peptide_range))
                has_ch = bool(_parse_ranges(chain_range))

                if active_sequence != sequence:
                    if has_pep:
                        active_used_peptide += 1
                    elif has_ch:
                        active_used_chain += 1
                    elif has_sp and has_pp:
                        active_both += 1
                    elif has_sp:
                        active_sp_only += 1
                    elif has_pp:
                        active_pp_only += 1
                else:
                    active_unchanged += 1

                f.write(f">{entry}\n")
                f.write(f"{active_sequence}\n")

    if verbose:
        total_modified = active_used_peptide + active_used_chain + active_sp_only + active_pp_only + active_both
        print(f"  Wrote {total_modified + active_unchanged} sequences")
        print(f"    - {active_used_peptide} from PEPTIDE annotations")
        print(f"    - {active_used_chain} from CHAIN annotations")
        print(f"    - {active_both} with SP + propeptide removed (subtraction)")
        print(f"    - {active_sp_only} with SP only removed")
        print(f"    - {active_pp_only} with propeptide only removed")
        print(f"    - {active_unchanged} unchanged (no annotations)")

    return fasta_full, fasta_mature, fasta_active


def print_next_steps(
    fasta_full: Path,
    fasta_mature: Path,
    fasta_active: Path,
    year: str = "2025",
) -> None:
    """Print instructions for generating embeddings via Colab.

    Args:
        fasta_full: Path to full sequences FASTA
        fasta_mature: Path to mature sequences FASTA (UniProt-based)
        fasta_active: Path to active sequences FASTA (SP + propeptide removed)
        year: Dataset year
    """
    colab_dir = fasta_full.parent
    print("\n" + "=" * 60)
    print("Next Steps: Generate embeddings via Google Colab")
    print("=" * 60)
    print(f"\n1. Upload FASTA files from {colab_dir}/ to Google Colab:")
    print(f"   - {fasta_full.name}")
    print(f"   - {fasta_mature.name}")
    print(f"   - {fasta_active.name}")
    print(
        "\n2. Open: https://colab.research.google.com/github/tsenoner/protspace/"
        "blob/master/colab/ProtSpace_Embeddings.ipynb"
    )
    print(f"\n3. Download H5 files to {colab_dir}/:")
    print(f"   - toxprot_{year}_full.h5")
    print(f"   - toxprot_{year}_mature.h5")
    print(f"   - toxprot_{year}_active.h5")
    print("\n4. Run: toxprot analysis protspace prepare")
