"""
Prepare metadata variants for protspace visualization.

Generates four metadata files for the 2025 dataset:
1. All data (top 15 protein families + NaN + Other)
2. Top 15 only (excluding NaN and Other, full sequences)
3. Top 15 mature (excluding NaN/Other, mature sequences only)
4. Top 15 mature no fragments (excluding NaN/Other/fragments, mature sequences)
"""

import pandas as pd
from pathlib import Path


def process_top_n_families(
    df: pd.DataFrame, column: str = "Protein families", top_n: int = 15
) -> pd.DataFrame:
    """
    Process protein families to keep top N actual families, label rest as 'Other'.

    Note: NaN values are preserved as NaN and excluded from the top N count.
    "Other" is the label for all families outside the top N.

    Args:
        df: Input DataFrame
        column: Column name to process
        top_n: Number of top families to keep (excluding NaN and Other)

    Returns:
        DataFrame with processed column
    """
    df = df.copy()

    # Get value counts excluding NaN - these are the actual protein families
    counts = df[column].value_counts()

    # Get top N actual families (this naturally excludes NaN since value_counts does)
    top_values = set(counts.nlargest(top_n).index)

    # Create processed column: keep top N, NaN stays NaN, rest becomes "Other"
    def categorize(x):
        if pd.isna(x):
            return x  # Keep NaN as is
        elif x in top_values:
            return x
        else:
            return "Other"

    df[column] = df[column].apply(categorize)

    return df


def create_metadata_variants(
    input_csv: Path, interim_tsv: Path, output_dir: Path, year: str = "2025"
):
    """
    Create four metadata variants from base CSV file.

    Args:
        input_csv: Path to processed CSV (e.g., toxprot_2025.csv)
        interim_tsv: Path to interim TSV with signal peptide info
        output_dir: Directory to save output files
        year: Dataset year (for naming)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read processed CSV
    print(f"Reading {input_csv}...")
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} entries")

    # Read interim TSV to get signal peptide information
    print(f"Reading signal peptide info from {interim_tsv}...")
    df_interim = pd.read_csv(
        interim_tsv, sep="\t", usecols=["Entry", "Signal peptide (range)"]
    )

    # Merge signal peptide information
    df = df.merge(df_interim, on="Entry", how="left")

    # Rename Entry to identifier (required by protspace)
    df = df.rename(columns={"Entry": "identifier"})

    # Add 'has_signal_peptide' column (yes/no)
    df["has_signal_peptide"] = (
        df["Signal peptide (range)"].notna().map({True: "yes", False: "no"})
    )

    # Add 'has_fragment' column (yes/no)
    df["has_fragment"] = (df["Fragment"].astype(str) == "fragment").map(
        {True: "yes", False: "no"}
    )

    # Process protein families to get top 15
    df = process_top_n_families(df, "Protein families", top_n=15)

    # Variant 1: All data (including NaN and Other)
    df_all = df.copy()
    output_all = output_dir / f"metadata_{year}_all.csv"

    # Keep only specified columns
    columns_all = [
        "identifier",
        "Protein families",
        "Phylum",
        "has_fragment",
        "has_signal_peptide",
    ]
    df_all[columns_all].to_csv(output_all, index=False)
    print(f"\n1. All data: {len(df_all)} entries → {output_all}")

    # Variant 2: Top 15 only (exclude NaN and Other)
    df_top15 = df[
        df["Protein families"].notna() & (df["Protein families"] != "Other")
    ].copy()
    output_top15 = output_dir / f"metadata_{year}_top15.csv"

    # Keep columns (still has has_fragment and has_signal_peptide as filtering hasn't been applied yet)
    columns_top15 = [
        "identifier",
        "Protein families",
        "Phylum",
        "has_fragment",
        "has_signal_peptide",
    ]
    df_top15[columns_top15].to_csv(output_top15, index=False)
    print(f"2. Top 15 (full sequences): {len(df_top15)} entries → {output_top15}")

    # Variant 3: Top 15 mature (exclude NaN/Other, use mature sequences - signal peptides cut off)
    df_top15_mature = df[
        df["Protein families"].notna() & (df["Protein families"] != "Other")
    ].copy()
    output_top15_mature = output_dir / f"metadata_{year}_top15_mature.csv"

    # Keep has_fragment and has_signal_peptide for reference
    columns_top15_mature = [
        "identifier",
        "Protein families",
        "Phylum",
        "has_fragment",
        "has_signal_peptide",
    ]
    df_top15_mature[columns_top15_mature].to_csv(output_top15_mature, index=False)
    print(
        f"3. Top 15 mature (cut off signal peptides): {len(df_top15_mature)} entries → {output_top15_mature}"
    )

    # Variant 4: Top 15 mature no fragments (exclude NaN/Other/fragments, use mature sequences)
    df_top15_mature_no_frag = df[
        df["Protein families"].notna()
        & (df["Protein families"] != "Other")
        & (df["has_fragment"] == "no")
    ].copy()
    output_top15_mature_no_frag = (
        output_dir / f"metadata_{year}_top15_mature_no_fragments.csv"
    )

    # Remove has_fragment since filtered, keep has_signal_peptide for reference
    columns_top15_mature_no_frag = [
        "identifier",
        "Protein families",
        "Phylum",
        "has_signal_peptide",
    ]
    df_top15_mature_no_frag[columns_top15_mature_no_frag].to_csv(
        output_top15_mature_no_frag, index=False
    )
    print(
        f"4. Top 15 mature no fragments: {len(df_top15_mature_no_frag)} entries → {output_top15_mature_no_frag}"
    )

    # Print summary statistics
    print("\n" + "=" * 60)
    print("Summary:")
    print("=" * 60)
    print(f"Total entries: {len(df)}")
    print(f"  → All data: {len(df_all)}")
    print(f"  → Top 15 only: {len(df_top15)} ({len(df_top15) / len(df) * 100:.1f}%)")
    print(
        f"  → Top 15 mature: {len(df_top15_mature)} ({len(df_top15_mature) / len(df) * 100:.1f}%)"
    )
    print(
        f"  → Top 15 mature no fragments: {len(df_top15_mature_no_frag)} ({len(df_top15_mature_no_frag) / len(df) * 100:.1f}%)"
    )

    # Show top 15 protein families (excluding Other and NaN)
    print("\nTop 15 Protein Families (excluding 'Other' and NaN):")
    all_counts = df_all["Protein families"].value_counts()

    # Filter out "Other" and get top 15 actual families
    top15_families = all_counts[all_counts.index != "Other"].head(15)

    for i, (family, count) in enumerate(top15_families.items(), 1):
        print(f"  {i:2d}. {family}: {count}")

    # Also show Other and NaN counts for context
    other_count = all_counts.get("Other", 0)
    nan_count = df_all["Protein families"].isna().sum()
    print("\n  Additional:")
    print(f"     Other: {other_count}")
    print(f"     NaN: {nan_count}")


def main():
    """Main execution function."""
    # Define paths
    base_dir = Path(__file__).resolve().parent.parent.parent
    input_csv = base_dir / "data" / "processed" / "toxprot_2025.csv"
    interim_tsv = base_dir / "data" / "interim" / "toxprot_2025.tsv"
    output_dir = base_dir / "data" / "processed" / "protspace"

    # Create metadata variants
    create_metadata_variants(input_csv, interim_tsv, output_dir, year="2025")

    print("\n✓ Metadata preparation complete!")


if __name__ == "__main__":
    main()
