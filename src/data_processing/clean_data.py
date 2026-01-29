#!/usr/bin/env python3
"""
Clean and process ToxProt data from TSV files.

This script:
1. Updates protein family names for consistency
2. Extracts PTM summaries from JSON (with resolved modification names)
3. Creates FASTA files (with signal peptide removal)
4. Adds taxonomic information using taxopy
5. Adds habitat classification (marine/terrestrial)
6. Exports cleaned CSV files
"""

import argparse
import json
import logging
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
import taxopy

from .parse_sprot_dat import PTMVocabulary

# Module logger - configuration handled by caller or main()
logger = logging.getLogger(__name__)

# Global PTM vocabulary instance, lazy-initialized
_ptm_vocab = None


def _get_ptm_vocab():
    """Get or initialize the global PTM vocabulary."""
    global _ptm_vocab
    if _ptm_vocab is None:
        _ptm_vocab = PTMVocabulary("data/raw")
        _ptm_vocab.ensure_ptmlist_available()
        _ptm_vocab.load_ptmlist()
    return _ptm_vocab


def extract_ptm_summary(ptm_json: str) -> str:
    """Extract PTM type counts from PTM_Features JSON.

    For "Modified residue" entries, resolves the modification name to its
    biological keyword (e.g., "4-hydroxyproline" -> "Hydroxylation").

    Args:
        ptm_json: JSON string with PTM data.

    Returns:
        Summary string like 'Disulfide bond:4; Hydroxylation:2; Amidation:1'.
    """
    if not ptm_json or pd.isna(ptm_json):
        return ""
    try:
        data = json.loads(ptm_json)
        counts = Counter()

        for ptm_type, features in data.items():
            if not isinstance(features, list):
                continue

            if ptm_type == "Modified residue":
                # Resolve each modification to its keyword
                vocab = _get_ptm_vocab()
                for feat in features:
                    note = feat.get("note", "")
                    if note:
                        keyword = vocab.resolve_mod_res_keyword(note)
                        counts[keyword] += 1
            else:
                counts[ptm_type] += len(features)

        if not counts:
            return ""

        # Sort by count descending, then alphabetically
        sorted_items = sorted(counts.items(), key=lambda x: (-x[1], x[0]))
        return "; ".join(f"{k}:{v}" for k, v in sorted_items)
    except (json.JSONDecodeError, TypeError):
        return ""


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
    with open(tsv_input_path) as f:
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
        "PTM_Features",
        "PTM Keywords",
        "Sequence",
        "Signal peptide (range)",
        "Protein existence",
        "ToxProt definition",
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

    # Update protein families
    df = update_protfams_func(df)

    # Extract PTM summaries from PTM_Features JSON and place next to PTM Keywords
    if "PTM_Features" in df.columns:
        ptm_summary = df["PTM_Features"].apply(extract_ptm_summary)
        # Insert PTM Summary right after PTM Keywords
        if "PTM Keywords" in df.columns:
            ptm_kw_idx = df.columns.get_loc("PTM Keywords")
            df.insert(ptm_kw_idx + 1, "PTM Summary", ptm_summary)
        else:
            df["PTM Summary"] = ptm_summary

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

    # Prepare columns for CSV output: remove columns not needed in final CSV
    # Keep 'Gene Ontology (molecular function)' from GO columns
    drop_cols = [
        "Sequence",
        "Signal peptide (range)",
        "PTM_Features",  # Full JSON not needed in final CSV
    ] + [col for col in go_merge_cols if col != "Gene Ontology (molecular function)"]
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
    with open(habitat_mapping_path) as f:
        habitat_mapping = json.load(f)

    # Load the detailed habitat mapping
    with open(habitat_detailed_path) as f:
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


def get_year_from_filename(filepath: Path) -> str:
    """Extract year from filename like toxprot_2005.tsv -> 2005"""
    import re
    match = re.search(r"toxprot_(\d{4})", filepath.stem)
    return match.group(1) if match else filepath.stem


def main():
    """Main function to process ToxProt data."""
    # Configure logging for standalone execution
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
        ],
    )

    parser = argparse.ArgumentParser(
        description="Clean and process ToxProt TSV data files.",
        epilog="""Examples:
  # Process all files in default input directory
  python clean_data.py

  # Process specific years
  python clean_data.py --years 2005 2010 2015 2020 2025

  # Process files from custom directory
  python clean_data.py --input-dir data/interim/toxprot_parsed""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("data/interim/toxprot_parsed"),
        help="Directory containing toxprot_*.tsv files (default: data/interim/toxprot_parsed)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed/toxprot"),
        help="Directory for output files (default: data/processed/toxprot)",
    )
    parser.add_argument(
        "--years",
        type=int,
        nargs="+",
        help="Specific years to process (default: all found in input-dir)",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("data"),
        help="Base data directory for habitat mappings (default: data)",
    )

    args = parser.parse_args()

    # Define paths for habitat mappings
    habitat_mapping_path = args.data_dir / "raw" / "marine_terrestrial.json"
    habitat_detailed_path = args.data_dir / "raw" / "habitat_detailed.json"

    # Ensure output directory exists
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Find input files
    if args.years:
        input_files = [args.input_dir / f"toxprot_{year}.tsv" for year in args.years]
        input_files = [f for f in input_files if f.exists()]
    else:
        input_files = sorted(args.input_dir.glob("toxprot_*.tsv"))

    if not input_files:
        logger.error(f"No toxprot_*.tsv files found in {args.input_dir}")
        return

    logger.info("=" * 60)
    logger.info("ToxProt Data Cleaning Pipeline")
    logger.info("=" * 60)
    logger.info(f"Input directory: {args.input_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Files to process: {len(input_files)}")
    logger.info("=" * 60)

    # Initialize taxonomy database
    logger.info("Initializing taxonomy database...")
    taxdb = initialize_taxdb()

    successful = 0
    for tsv_path in input_files:
        year = get_year_from_filename(tsv_path)
        logger.info(f"[{year}] Processing {tsv_path.name}...")

        # Step 1: Process TSV to CSV and FASTA (in input directory as intermediate)
        logger.info("  Processing TSV and creating FASTA...")
        process_toxprot_tsv(tsv_path, update_protfams, create_fasta_file)

        # Intermediate files created next to input
        interim_csv_path = tsv_path.with_suffix(".csv")
        interim_fasta_path = tsv_path.with_suffix(".fasta")

        # Final output paths
        processed_csv_path = args.output_dir / f"toxprot_{year}.csv"
        processed_fasta_path = args.output_dir / f"toxprot_{year}.fasta"

        # Step 2: Add taxonomic information
        logger.info("  Adding taxonomy...")
        df = process_dataframe_with_taxonomy(
            interim_csv_path, processed_csv_path, taxdb
        )

        # Step 3: Add habitat classification
        logger.info("  Adding habitat classification...")
        df = add_habitat_classification(
            processed_csv_path, habitat_mapping_path, habitat_detailed_path
        )

        # Move FASTA to output directory
        interim_fasta_path.rename(processed_fasta_path)

        # Clean up intermediate CSV
        interim_csv_path.unlink()

        # Log summary statistics
        logger.info(f"    Entries: {len(df)}")
        logger.info(f"    Families: {df['Protein families'].nunique()}")
        habitat_counts = df["Habitat"].value_counts().to_dict()
        logger.info(f"    Habitat: {habitat_counts}")

        successful += 1

    logger.info("=" * 60)
    logger.info(f"Processing complete! {successful}/{len(input_files)} files processed")
    logger.info(f"  Output: {args.output_dir}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
