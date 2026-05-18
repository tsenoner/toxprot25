"""Prepare metadata and H5 file variants for ProtSpace visualization.

This module handles:
1. Creating metadata CSV files for each variant. Analysis variants ship a lean
   schema with top-N + Other family grouping; demo variants ship an expanded
   schema with full family names and additional biological columns.
2. Filtering H5 embedding files to match metadata variants.
"""

from pathlib import Path

import h5py
import pandas as pd

from ..analyze_protein_families import get_reference_families, normalize_family_name
from .config import (
    COLAB_SUBDIR,
    DEMO_ANNOTATION_COLUMNS,
    INTERMEDIATES_SUBDIR,
    SEQ_TYPES,
    TOP_N,
    VARIANT_CONFIGS,
    get_h5_base_filename,
    get_h5_variant_filename,
    get_metadata_filename,
)

# Analysis-variant annotation columns (lean schema with top-N + Other grouping).
ANALYSIS_BASE_COLUMNS = [
    "identifier",
    "Protein families",
    "Phylum",
    "has_fragment",
    "has_signal_peptide",
    "has_propeptide",
]


def process_protein_families(
    df: pd.DataFrame,
    reference_families: list[str],
    column: str = "Protein families",
) -> pd.DataFrame:
    """Collapse protein-family values to the top-N reference set + "Other".

    NaN values are preserved as-is (the protspace UI surfaces them as a
    distinct category).
    """
    df = df.copy()

    def categorize(value):
        if pd.isna(value):
            return value
        normalized = normalize_family_name(value)
        if normalized in reference_families:
            return normalized
        return "Other"

    df[column] = df[column].apply(categorize)
    return df


def normalize_family_names_full(df: pd.DataFrame, column: str = "Protein families") -> pd.DataFrame:
    """Normalize family names without collapsing the long tail to "Other"."""
    df = df.copy()
    df[column] = df[column].apply(lambda v: v if pd.isna(v) else normalize_family_name(v))
    return df


def _select_columns(variant_config: dict) -> list[str]:
    """Return the annotation source-column list for a variant.

    The returned names are CSV-source names. For demo variants the caller
    renames them to canonical UI names via ``DEMO_ANNOTATION_COLUMNS``.

    Drops columns that become meaningless after filtering:
    - mature/active variants: drop ``has_signal_peptide`` (SP already cleaved)
    - active variants: also drop ``has_propeptide`` (propeptide already cleaved)
    - variants with ``exclude_fragments``: drop ``has_fragment`` (now constant)
    """
    if variant_config["bundle_kind"] == "demo":
        columns = list(DEMO_ANNOTATION_COLUMNS.keys())
    else:
        columns = list(ANALYSIS_BASE_COLUMNS)

    seq_type = variant_config["seq_type"]

    if variant_config["exclude_fragments"] and "has_fragment" in columns:
        columns.remove("has_fragment")
    if seq_type in ("mature", "active") and "has_signal_peptide" in columns:
        columns.remove("has_signal_peptide")
    if seq_type == "active" and "has_propeptide" in columns:
        columns.remove("has_propeptide")

    return columns


def create_metadata_csv(
    df: pd.DataFrame,
    variant_config: dict,
    reference_families: list[str],
    output_path: Path,
    verbose: bool = True,
) -> int:
    """Create the metadata CSV for a single variant.

    For analysis variants the families are collapsed to top-N + Other; for
    demo variants the full curated family names are kept.
    """
    df_variant = df.copy()

    if variant_config["exclude_nan"]:
        df_variant = df_variant[df_variant["Protein families"].notna()]
    if variant_config["exclude_other"]:
        df_variant = df_variant[df_variant["Protein families"] != "Other"]
    if variant_config["exclude_fragments"]:
        df_variant = df_variant[df_variant["has_fragment"] == "no"]

    if variant_config["bundle_kind"] == "demo":
        df_variant = normalize_family_names_full(df_variant)
    else:
        df_variant = process_protein_families(df_variant, reference_families)

    columns = _select_columns(variant_config)
    missing = [c for c in columns if c not in df_variant.columns]
    if missing:
        raise KeyError(
            f"Variant {variant_config['name']} expects columns {missing} "
            f"but they are not present in the input frame."
        )

    df_out = df_variant[columns]
    if variant_config["bundle_kind"] == "demo":
        df_out = df_out.rename(columns={k: DEMO_ANNOTATION_COLUMNS[k] for k in columns})

    df_out.to_csv(output_path, index=False)

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

    Also writes a ``model_name`` attribute (``prot_t5``) to the output file
    because protspace v4 requires it.
    """
    if not h5_input.exists():
        raise FileNotFoundError(f"H5 file not found: {h5_input}")
    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_csv}")

    df = pd.read_csv(metadata_csv)
    identifiers_to_keep = set(df["identifier"].tolist())

    with h5py.File(h5_input, "r") as input_file:
        h5_keys = set(input_file.keys())
        proteins_to_keep = identifiers_to_keep.intersection(h5_keys)

        h5_output.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(h5_output, "w") as output_file:
            for protein_id in proteins_to_keep:
                input_file.copy(protein_id, output_file)
            # protspace v4 requires a model_name attribute on the H5 file
            # so it can name projections (e.g. "ProtT5 — UMAP 2").
            output_file.attrs["model_name"] = "prot_t5"

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
    variants: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, dict]:
    """Prepare metadata and H5 files for the requested variants.

    Args:
        processed_csv: Path to the processed CSV (e.g., toxprot_2025.csv).
        interim_tsv: Path to interim TSV with SP/propeptide annotations.
        protspace_dir: Root directory for protspace files.
        year: Dataset year.
        top_n: Number of top families for analysis variants.
        definition: ToxProt definition filter.
        variants: Subset of variant names to prepare (default: all).
        verbose: Print progress messages.
    """
    protspace_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print(f"Loading {processed_csv}...")
    df = pd.read_csv(processed_csv)

    if "ToxProt definition" in df.columns and definition == "venom_tissue":
        df = df[df["ToxProt definition"].isin(["venom_tissue", "both"])]

    if verbose:
        print(f"Loaded {len(df)} entries (definition: {definition})")

    reference_families = get_reference_families(df, top_n=top_n)
    if verbose:
        print(f"\nTop {top_n} protein families:")
        for i, fam in enumerate(reference_families, 1):
            print(f"  {i:2d}. {fam}")

    # Merge in signal-peptide and propeptide range info from the interim TSV.
    if verbose:
        print(f"\nLoading signal peptide/propeptide info from {interim_tsv}...")
    with open(interim_tsv) as f:
        tsv_header = f.readline().strip().split("\t")
    interim_cols = ["Entry", "Signal peptide (range)"]
    if "Propeptide (range)" in tsv_header:
        interim_cols.append("Propeptide (range)")
    df_interim = pd.read_csv(interim_tsv, sep="\t", usecols=interim_cols)
    df = df.merge(df_interim, on="Entry", how="left")

    # protspace conventionally calls the identifier column `identifier`.
    df = df.rename(columns={"Entry": "identifier"})

    df["has_signal_peptide"] = df["Signal peptide (range)"].notna().map({True: "yes", False: "no"})
    if "Propeptide (range)" in df.columns:
        df["has_propeptide"] = df["Propeptide (range)"].notna().map({True: "yes", False: "no"})
    else:
        df["has_propeptide"] = "no"
    df["has_fragment"] = (df["Fragment"].astype(str) == "fragment").map({True: "yes", False: "no"})

    # Resolve which base H5 files we need.
    colab_dir = protspace_dir / COLAB_SUBDIR
    h5_files = {st: colab_dir / get_h5_base_filename(year, st) for st in SEQ_TYPES}

    selected = list(variants) if variants else list(VARIANT_CONFIGS.keys())
    needed_seq_types = {VARIANT_CONFIGS[v]["seq_type"] for v in selected if v in VARIANT_CONFIGS}
    h5_missing = [h5_files[st] for st in needed_seq_types if not h5_files[st].exists()]
    if h5_missing:
        _print_h5_missing_error(h5_missing, year)
        raise FileNotFoundError("Required H5 embedding files not found")

    intermediates_dir = protspace_dir / INTERMEDIATES_SUBDIR
    intermediates_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("\nCreating metadata and H5 variants...")

    results = {}
    for variant_name in selected:
        if variant_name not in VARIANT_CONFIGS:
            if verbose:
                print(f"\nUnknown variant: {variant_name}")
            continue

        config = VARIANT_CONFIGS[variant_name]
        if verbose:
            print(f"\n{config['description']}:")

        h5_base = h5_files.get(config["seq_type"])
        if h5_base is None or not h5_base.exists():
            if verbose:
                print(f"  Skipping: H5 file not found ({h5_base})")
            continue

        metadata_path = intermediates_dir / get_metadata_filename(year, variant_name)
        n_entries = create_metadata_csv(
            df, config, reference_families, metadata_path, verbose=verbose
        )

        h5_variant = intermediates_dir / get_h5_variant_filename(year, variant_name)
        n_kept, _ = filter_h5_by_metadata(h5_base, metadata_path, h5_variant, verbose=verbose)

        results[variant_name] = {
            "metadata": metadata_path,
            "h5": h5_variant,
            "n_entries": n_entries,
            "n_embeddings": n_kept,
        }

    if verbose:
        print("\n" + "=" * 60)
        print("Summary:")
        print("=" * 60)
        for variant_name, info in results.items():
            print(
                f"  {variant_name}: {info['n_entries']} entries, "
                f"{info['n_embeddings']} embeddings"
            )

    return results


def _print_h5_missing_error(missing_files: list[Path], year: str) -> None:
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
