#!/usr/bin/env python3
"""
Clean and process ToxProt data from TSV files.

This script:
1. Updates protein family names for consistency
2. Creates FASTA files (with signal peptide removal)
3. Adds taxonomic information using taxopy
4. Adds habitat classification (marine/terrestrial)
5. Exports cleaned CSV files
"""

import argparse
import json
from pathlib import Path

import pandas as pd
import taxopy


def update_protfams(df):
    """Update and standardize protein family names."""
    # Split on common delimiters and take first part
    df["Protein families"] = df["Protein families"].str.split(r"[.,;]").str[0]

    # Map of family name corrections
    family_corrections = {
        "I1 superfamily": "Conotoxin I1 superfamily",
        "O1 superfamily": "Conotoxin O1 superfamily",
        "O2 superfamily": "Conotoxin O2 superfamily",
        "E superfamily": "Conotoxin E superfamily",
        "F superfamily": "Conotoxin F superfamily",
        "Conotoxin M family": "Conotoxin M superfamily",
        "Conotoxin B2 family": "Conotoxin B2 superfamily",
        "Conotoxin O1 family": "Conotoxin O1 superfamily",
        "Conotoxin O2 family": "Conotoxin O2 superfamily",
        "Bradykinin- potentiating peptide family": "Bradykinin-potentiating peptide family",
    }

    # Apply all corrections at once
    df["Protein families"] = df["Protein families"].replace(family_corrections)

    return df


def create_fasta_file(
    df: pd.DataFrame,
    entry_col: str,
    sequence_col: str,
    fasta_output_path: Path,
    signal_peptide_range_column: str = None,
):
    """Create a FASTA file from a DataFrame, optionally removing signal peptides."""
    fasta_output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(fasta_output_path, "w") as f_out:
        for _, row in df.iterrows():
            entry = row[entry_col]
            original_sequence = row[sequence_col]

            # Ensure sequence is a valid non-empty string
            if not isinstance(original_sequence, str) or not original_sequence:
                continue

            sequence_to_write = original_sequence

            # If a signal peptide range column is provided, attempt to remove the signal peptide
            if signal_peptide_range_column:
                signal_range_str = row.get(signal_peptide_range_column)

                # Process if range string is a valid string and contains a hyphen (e.g., "1-22")
                if isinstance(signal_range_str, str) and "-" in signal_range_str:
                    try:
                        # Signal peptide range (e.g., "1-22") is 1-based; slice at end_pos_1based for mature protein.
                        end_pos_1based = int(signal_range_str.split("-")[1])

                        if end_pos_1based > 0:
                            # Slice from end_pos_1based. Python handles end_pos_1based >= len correctly (empty string).
                            sequence_to_write = original_sequence[end_pos_1based:]
                        # If end_pos_1based is not positive, original_sequence is kept (invalid range for cut).
                    except (ValueError, IndexError):
                        # If parsing fails (e.g. "1-foo", "1-", "-"), keep the original sequence.
                        pass

            # Write to FASTA only if the (potentially modified) sequence is not empty
            if sequence_to_write:
                f_out.write(f">{entry}\n")
                f_out.write(f"{sequence_to_write}\n")


def process_toxprot_tsv(tsv_input_path: Path, update_protfams_func, create_fasta_func):
    """
    Process a ToxProt TSV file:
    - Reads the TSV
    - Renames columns for consistency
    - Updates protein families
    - Merges Gene Ontology columns
    - Writes a cleaned CSV (without sequence and original GO columns, except 'Gene Ontology (molecular function)')
    - Writes a FASTA file (removing signal peptide if present)
    - Displays summary info

    Args:
        tsv_input_path (Path): Path to the input TSV file.
        update_protfams_func (callable): Function to update protein families.
        create_fasta_func (callable): Function to create FASTA file.
    """
    # Infer output paths
    base = tsv_input_path.with_suffix("")
    csv_output_path = base.with_suffix(".csv")
    fasta_output_path = base.with_suffix(".fasta")

    # Read columns from the TSV file
    # Try to infer columns present in the file
    with open(tsv_input_path, "r") as f:
        header = f.readline().strip().split("\t")
    # Required columns
    required_cols = [
        "Entry",
        "Organism",
        "Organism (ID)",
        "Protein families",
        "Length",
        "Fragment",
        "Toxic dose",
        "Post-translational modification",
        "PTM Summary",
        "Sequence",
        "Signal peptide (range)",
        "Protein existence",
    ]
    # Optional GO columns
    go_cols = [
        "Gene Ontology (GO)",
        "Gene Ontology (biological process)",
        "Gene Ontology (cellular component)",
        "Gene Ontology (molecular function)",
    ]
    # Only use columns that exist in the file
    usecols = [col for col in required_cols + go_cols if col in header]

    df = pd.read_csv(tsv_input_path, sep="\t", usecols=usecols)

    # Rename 'Post-translational modification' to 'PTM' if present
    if "Post-translational modification" in df.columns:
        df = df.rename(columns={"Post-translational modification": "PTM"})

    # Update protein families
    df = update_protfams_func(df)

    # Merge GO columns if present
    go_merge_cols = [
        col
        for col in [
            "Gene Ontology (biological process)",
            "Gene Ontology (cellular component)",
            "Gene Ontology (molecular function)",
        ]
        if col in df.columns
    ]
    if go_merge_cols:

        def merge_go_terms(row):
            terms = []
            for col in go_merge_cols:
                val = row.get(col)
                if pd.notnull(val):
                    val = str(val).strip()
                    if val:
                        terms.append(val)
            return "; ".join(terms) if terms else pd.NA

        df["Gene Ontology (GO)"] = df.apply(merge_go_terms, axis=1)

    # Create FASTA file if 'Sequence' and 'Entry' columns exist
    if "Sequence" in df.columns and "Entry" in df.columns:
        signal_peptide_col = (
            "Signal peptide (range)" if "Signal peptide (range)" in df.columns else None
        )
        create_fasta_func(
            df, "Entry", "Sequence", fasta_output_path, signal_peptide_col
        )

    # Prepare columns for CSV output: remove 'Sequence', 'Signal peptide (range)', and original GO columns,
    # but KEEP 'Gene Ontology (molecular function)' in the CSV
    drop_cols = ["Sequence", "Signal peptide (range)"] + [
        col for col in go_merge_cols if col != "Gene Ontology (molecular function)"
    ]
    columns_for_csv = [col for col in df.columns if col not in drop_cols]

    # Save CSV
    df.to_csv(csv_output_path, index=False, columns=columns_for_csv)

    return df


def setup_db_paths():
    """Setup and return the database paths."""
    home_dir = Path.home() / ".cache"
    db_dir = home_dir / "taxopy_db"
    db_dir.mkdir(parents=True, exist_ok=True)
    nodes_file = db_dir / "nodes.dmp"
    names_file = db_dir / "names.dmp"
    merged_file = db_dir / "merged.dmp"

    return db_dir, nodes_file, names_file, merged_file


def initialize_taxdb():
    """Initialize and return the taxonomy database."""
    # Get the database paths
    db_dir, nodes_file, names_file, merged_file = setup_db_paths()

    if nodes_file.exists() and names_file.exists():
        taxdb = taxopy.TaxDb(
            nodes_dmp=str(nodes_file),
            names_dmp=str(names_file),
            merged_dmp=str(merged_file),
        )
    else:
        taxdb = taxopy.TaxDb(taxdb_dir=str(db_dir), keep_files=True)

    return taxdb


def get_taxonomy_info(taxon_id, taxdb):
    """Get order, family, genus, species info for a taxon ID."""
    # Get the Taxon object
    taxon = taxopy.Taxon(taxon_id, taxdb)

    # Get the rank information
    ranks = taxon.rank_name_dictionary

    return {
        "taxon_name": taxon.name,
        "phylum": ranks.get("phylum", ""),
        "class": ranks.get("class", ""),
        "order": ranks.get("order", ""),
        "family": ranks.get("family", ""),
        "genus": ranks.get("genus", ""),
        "species": ranks.get("species", ""),
    }


# Mapping of deprecated/merged NCBI taxonomy IDs to current IDs
LEGACY_TAXID_MAP = {
    184771: 666126,  # Oxyopes kitabensis -> Oxyopes takobius
}


def build_taxonomy_cache(df, taxdb):
    """Build a cache of taxonomy information for all unique organism IDs."""
    taxonomy_cache = {}
    for taxon_id in df["Organism (ID)"].unique():
        if pd.notna(taxon_id):
            # Map legacy taxonomy IDs to current ones
            lookup_id = LEGACY_TAXID_MAP.get(int(taxon_id), int(taxon_id))
            taxonomy_cache[taxon_id] = get_taxonomy_info(lookup_id, taxdb)

    return taxonomy_cache


def add_taxonomy_columns(df, taxonomy_cache):
    """Add taxonomy columns to the dataframe."""

    # Create a mapping function that extracts all taxonomy info at once
    def get_taxonomy_info_from_cache(taxon_id):
        cache_entry = taxonomy_cache.get(taxon_id, {})
        return pd.Series(
            {
                "Scientific_Name": cache_entry.get("taxon_name", ""),
                "Phylum": cache_entry.get("phylum", ""),
                "Class": cache_entry.get("class", ""),
                "Order": cache_entry.get("order", ""),
                "Family": cache_entry.get("family", ""),
                "Genus": cache_entry.get("genus", ""),
                "Species": cache_entry.get("species", ""),
            }
        )

    # Apply the mapping function once to get all columns
    taxonomy_df = df["Organism (ID)"].apply(get_taxonomy_info_from_cache)

    # Concatenate the new columns with the original dataframe
    return pd.concat([df, taxonomy_df], axis=1)


def process_dataframe_with_taxonomy(input_path, output_path, taxdb):
    """Process the dataframe: load, add taxonomy, remove Organism column, save."""
    # Load the dataframe
    df = pd.read_csv(input_path)

    # Build the taxonomy cache
    taxonomy_cache = build_taxonomy_cache(df, taxdb)

    # Add taxonomy columns
    df = add_taxonomy_columns(df, taxonomy_cache)

    # Remove the Organism column if it exists
    if "Organism" in df.columns:
        df = df.drop(columns=["Organism"])

    # Save the updated dataframe
    df.to_csv(output_path, index=False)

    return df


def determine_habitat(row, habitat_mapping):
    """Determine habitat based on order and genus."""
    order = row["Order"]
    genus = row.get("Genus", "")  # Get genus if available, otherwise empty string

    # Check if order is in clear_orders
    if order in habitat_mapping["clear_orders"]["terrestrial"]:
        return "terrestrial"
    elif order in habitat_mapping["clear_orders"]["marine"]:
        return "marine"

    # Check if order is in ambiguous_orders
    if order in habitat_mapping["ambiguous_orders"]:
        # Check if genus is in the terrestrial list for this order
        if genus in habitat_mapping["ambiguous_orders"][order].get("terrestrial", {}):
            return "terrestrial"
        # Check if genus is in the marine list for this order
        elif genus in habitat_mapping["ambiguous_orders"][order].get("marine", {}):
            return "marine"

    # If we can't determine, return 'unknown'
    return "unknown"


def determine_habitat_detailed(row, habitat_mapping):
    """Determine detailed habitat (terrestrial/freshwater/estuarine/marine) based on order and genus."""
    order = row["Order"]
    genus = row.get("Genus", "")  # Get genus if available, otherwise empty string

    # Check if order is in clear_orders
    for habitat_type in ["terrestrial", "freshwater", "estuarine", "marine"]:
        if order in habitat_mapping["clear_orders"].get(habitat_type, {}):
            return habitat_type

    # Check if order is in ambiguous_orders
    if order in habitat_mapping["ambiguous_orders"]:
        # Check each habitat type for this order
        for habitat_type in ["terrestrial", "freshwater", "estuarine", "marine"]:
            if genus in habitat_mapping["ambiguous_orders"][order].get(
                habitat_type, {}
            ):
                return habitat_type

    # If we can't determine, return 'unknown'
    return "unknown"


def add_habitat_classification(csv_path, habitat_mapping_path, habitat_detailed_path):
    """Add habitat classification to the CSV file."""
    # Load the marine/terrestrial mapping
    with open(habitat_mapping_path, "r") as f:
        habitat_mapping = json.load(f)

    # Load the detailed habitat mapping
    with open(habitat_detailed_path, "r") as f:
        habitat_detailed_mapping = json.load(f)

    # Load the CSV
    df = pd.read_csv(csv_path)

    # Add habitat column (simple: marine/terrestrial)
    df["Habitat"] = df.apply(
        lambda row: determine_habitat(row, habitat_mapping), axis=1
    )

    # Add detailed habitat column (terrestrial/freshwater/estuarine/marine)
    df["Habitat_Detailed"] = df.apply(
        lambda row: determine_habitat_detailed(row, habitat_detailed_mapping), axis=1
    )

    # Save updated CSV
    df.to_csv(csv_path, index=False)

    return df


def main():
    """Main function to process ToxProt data."""
    parser = argparse.ArgumentParser(
        description="Clean and process ToxProt TSV data files."
    )
    parser.add_argument(
        "--year",
        choices=["2005", "2015", "2017", "2025", "all"],
        default="all",
        help="Which year(s) to process (default: all)",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data"),
        help="Base data directory (default: data)",
    )

    args = parser.parse_args()

    # Define paths
    data_path = args.data_dir
    interim_path = data_path / "interim"
    processed_path = data_path / "processed"
    habitat_mapping_path = data_path / "raw" / "marine_terrestrial.json"
    habitat_detailed_path = data_path / "raw" / "habitat_detailed.json"

    # Ensure directories exist
    processed_path.mkdir(parents=True, exist_ok=True)

    # Determine which years to process
    all_years = ["2005", "2015", "2017", "2025"]
    if args.year == "all":
        years_to_process = all_years
    else:
        years_to_process = [args.year]

    # Initialize taxonomy database
    print("Initializing taxonomy database...")
    taxdb = initialize_taxdb()

    for year in years_to_process:
        print(f"\nProcessing ToxProt {year}...")

        # Step 1: Process TSV to CSV and FASTA
        tsv_path = interim_path / f"toxprot_{year}.tsv"
        csv_path = interim_path / f"toxprot_{year}.csv"
        fasta_path = interim_path / f"toxprot_{year}.fasta"

        print("  → Processing TSV and creating FASTA...")
        process_toxprot_tsv(tsv_path, update_protfams, create_fasta_file)
        print(f"    ✓ CSV: {csv_path}")
        print(f"    ✓ FASTA: {fasta_path}")

        # Step 2: Add taxonomic information
        interim_csv_path = interim_path / f"toxprot_{year}.csv"
        processed_csv_path = processed_path / f"toxprot_{year}.csv"

        print("  → Adding taxonomy and habitat classification...")
        df = process_dataframe_with_taxonomy(
            interim_csv_path, processed_csv_path, taxdb
        )

        # Step 3: Add habitat classification
        df = add_habitat_classification(
            processed_csv_path, habitat_mapping_path, habitat_detailed_path
        )

        # Print summary statistics
        print(f"    ✓ Processed: {len(df)} entries")
        print(f"    ✓ Unique families: {df['Protein families'].nunique()}")
        habitat_counts = {k: int(v) for k, v in df["Habitat"].value_counts().items()}
        print(f"    ✓ Habitat: {habitat_counts}")
        habitat_detailed_counts = {
            k: int(v) for k, v in df["Habitat_Detailed"].value_counts().items()
        }
        print(f"    ✓ Habitat (detailed): {habitat_detailed_counts}")
        print(f"    ✓ Output: {processed_csv_path}")

    print(f"\n{'=' * 60}")
    print("✓ All processing complete!")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
