"""Analyze clustering quality using silhouette scores.

Compares clustering quality across ProtSpace variants using the
silhouette score, which measures how well-separated clusters are.
Score ranges from -1 to 1, where higher is better.

Supports single-file parquetbundle format (delimiter-separated parts).
"""

import io
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.metrics import silhouette_score

from .config import (
    VARIANT_CONFIGS,
    get_protspace_styled_filename,
)

# Delimiter used in parquetbundle files to separate parquet parts
PARQUET_BUNDLE_DELIMITER = b"---PARQUET_DELIMITER---"


def load_protspace_data(data_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load UMAP embeddings and protein family labels from protspace data.

    Supports single-file parquetbundle format (delimiter-separated parts).

    Args:
        data_path: Path to protspace styled parquetbundle file

    Returns:
        Tuple of (embeddings array, labels array)
    """
    if data_path.suffix == ".parquetbundle":
        return _load_from_parquetbundle(data_path)
    elif data_path.suffix == ".json":
        return _load_from_json(data_path)
    else:
        raise ValueError(f"Unsupported file format: {data_path.suffix}")


def _load_from_parquetbundle(bundle_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load data from single-file parquetbundle format.

    The parquetbundle contains multiple parquet parts separated by delimiter:
    - Part 0: Annotations (protein_id, Protein families, etc.)
    - Part 1: Projection metadata
    - Part 2: Coordinates (identifier, x, y, z)
    - Part 3: Settings (optional)

    Args:
        bundle_path: Path to parquetbundle file

    Returns:
        Tuple of (embeddings array, labels array)
    """
    # Read and split the bundle file
    with open(bundle_path, "rb") as f:
        content = f.read()

    parts = content.split(PARQUET_BUNDLE_DELIMITER)

    if len(parts) < 3:
        raise ValueError(
            f"Invalid parquetbundle format: expected at least 3 parts, got {len(parts)}"
        )

    # Part 0: Annotations
    df_ann = pq.ParquetFile(io.BytesIO(parts[0])).read().to_pandas()

    # Part 2: Coordinates
    df_coords = pq.ParquetFile(io.BytesIO(parts[2])).read().to_pandas()

    # Extract embeddings (x, y coordinates)
    embeddings = df_coords[["x", "y"]].values

    # Get protein families from annotations
    # Merge on identifier/protein_id
    id_col_ann = "protein_id" if "protein_id" in df_ann.columns else "identifier"
    id_col_coords = "identifier" if "identifier" in df_coords.columns else "protein_id"

    df_merged = df_coords[[id_col_coords]].merge(
        df_ann[[id_col_ann, "Protein families"]],
        left_on=id_col_coords,
        right_on=id_col_ann,
        how="left",
    )
    labels = df_merged["Protein families"].fillna("Unknown").values

    return embeddings, labels


def _load_from_json(json_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load data from legacy JSON format.

    Args:
        json_path: Path to protspace styled JSON file

    Returns:
        Tuple of (embeddings array, labels array)
    """
    with open(json_path) as f:
        data = json.load(f)

    embeddings = []
    labels = []

    projections_list = data.get("projections", [])
    protein_data = data.get("protein_data", {})

    # Find UMAP projection
    umap_projection = None
    for proj in projections_list:
        if proj.get("name") == "UMAP_2":
            umap_projection = proj
            break

    if umap_projection is None:
        raise ValueError("UMAP projection 'UMAP_2' not found in JSON file")

    # Extract coordinates and labels
    for item in umap_projection.get("data", []):
        protein_id = item.get("identifier")
        coords_dict = item.get("coordinates", {})
        coords = [coords_dict.get("x", 0), coords_dict.get("y", 0)]
        embeddings.append(coords)

        features = protein_data.get(protein_id, {}).get("features", {})
        family = features.get("Protein families", "Unknown")
        labels.append(family)

    return np.array(embeddings), np.array(labels)


def calculate_silhouette_score(
    embeddings: np.ndarray,
    labels: np.ndarray,
    exclude_labels: list[str] | None = None,
) -> float:
    """Calculate silhouette score for clustering.

    Args:
        embeddings: 2D array of coordinates
        labels: 1D array of cluster labels
        exclude_labels: Labels to exclude from calculation

    Returns:
        Silhouette score, or NaN if calculation not possible
    """
    if exclude_labels is None:
        exclude_labels = ["Unknown", "nan", "Other", "NaN"]

    # Filter out excluded labels
    valid_mask = ~np.isin(labels, exclude_labels)

    if valid_mask.sum() < 2:
        return np.nan

    embeddings_filtered = embeddings[valid_mask]
    labels_filtered = labels[valid_mask]

    # Need at least 2 unique labels
    unique_labels = np.unique(labels_filtered)
    if len(unique_labels) < 2:
        return np.nan

    try:
        return silhouette_score(embeddings_filtered, labels_filtered)
    except Exception:
        return np.nan


def analyze_variant(
    data_path: Path,
    variant_name: str,
    verbose: bool = True,
) -> dict:
    """Analyze clustering quality for a single variant.

    Args:
        data_path: Path to protspace styled parquetbundle or JSON
        variant_name: Name of the variant
        verbose: Print progress messages

    Returns:
        Dictionary with analysis results
    """
    if not data_path.exists():
        if verbose:
            print(f"  Data file not found: {data_path}")
        return {
            "variant": variant_name,
            "n_proteins": 0,
            "n_families": 0,
            "silhouette_score": np.nan,
        }

    embeddings, labels = load_protspace_data(data_path)

    n_proteins = len(labels)
    n_families = len(np.unique(labels))
    score = calculate_silhouette_score(embeddings, labels)

    if verbose:
        print(f"  Proteins: {n_proteins}")
        print(f"  Families: {n_families}")
        score_str = f"{score:.3f}" if not np.isnan(score) else "N/A"
        print(f"  Silhouette score: {score_str}")

    return {
        "variant": variant_name,
        "n_proteins": n_proteins,
        "n_families": n_families,
        "silhouette_score": score,
    }


def analyze_all_variants(
    protspace_dir: Path,
    year: str = "2025",
    variants: list[str] | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Analyze clustering quality for all variants.

    Args:
        protspace_dir: Directory containing styled parquetbundle files
        year: Dataset year
        variants: List of variant names to analyze (None = all variants)
        verbose: Print progress messages

    Returns:
        DataFrame with analysis results
    """
    # Default: analyze all variants
    if variants is None:
        variants = list(VARIANT_CONFIGS.keys())

    if verbose:
        print("=" * 60)
        print(f"Silhouette Score Analysis - {year}")
        print("=" * 60)

    results = []
    for variant_name in variants:
        if variant_name not in VARIANT_CONFIGS:
            if verbose:
                print(f"\nUnknown variant: {variant_name}")
            continue

        config = VARIANT_CONFIGS[variant_name]
        if verbose:
            print(f"\n{config['description']}:")

        data_path = protspace_dir / get_protspace_styled_filename(year, variant_name)
        result = analyze_variant(data_path, variant_name, verbose=verbose)
        result["description"] = config["description"]
        results.append(result)

    return pd.DataFrame(results)


def plot_silhouette_comparison(
    df: pd.DataFrame,
    output_path: Path,
    verbose: bool = True,
) -> bool:
    """Create bar chart visualization of silhouette scores.

    Args:
        df: DataFrame with analysis results
        output_path: Path to save plot
        verbose: Print progress messages

    Returns:
        True if successful
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Filter out NaN scores and prepare display names
    df_plot = df[df["silhouette_score"].notna()].copy()
    df_plot = df_plot.iloc[::-1]  # Reverse for horizontal bar chart

    # Create display labels
    df_plot["display_name"] = df_plot["variant"].str.replace("_", "\n")

    # Create horizontal bar chart
    ax.barh(
        df_plot["display_name"],
        df_plot["silhouette_score"],
        color="#b1b3c6",
        zorder=3,
    )

    # Styling
    ax.set_xlabel("Silhouette Score", fontsize=18)
    ax.set_title(
        "Clustering Quality Comparison\n(Higher score = better separation)",
        fontsize=18,
        fontweight="bold",
    )
    ax.set_xlim(0, 1)
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.grid(axis="x", alpha=0.5, linestyle="--", linewidth=2, zorder=0)

    # Add score labels
    for i, (_, row) in enumerate(df_plot.iterrows()):
        score = row["silhouette_score"]
        ax.text(score + 0.02, i, f"{score:.3f}", va="center", fontsize=18)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

    if verbose:
        print(f"\nPlot saved to: {output_path}")

    return True


def run_silhouette_analysis(
    protspace_dir: Path,
    figures_dir: Path,
    year: str = "2025",
    variants: list[str] | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Run complete silhouette analysis and generate outputs.

    Args:
        protspace_dir: Directory containing styled JSON files
        figures_dir: Directory to save outputs
        year: Dataset year
        variants: List of variant names to analyze
        verbose: Print progress messages

    Returns:
        DataFrame with analysis results
    """
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Analyze variants
    df_results = analyze_all_variants(protspace_dir, year=year, variants=variants, verbose=verbose)

    # Save CSV
    csv_path = figures_dir / "silhouette_comparison.csv"
    df_results.to_csv(csv_path, index=False)
    if verbose:
        print(f"\nResults saved to: {csv_path}")

    # Generate plot
    plot_path = figures_dir / "silhouette_comparison.png"
    plot_silhouette_comparison(df_results, plot_path, verbose=verbose)

    # Print summary
    if verbose:
        print("\n" + "=" * 60)
        print("Summary")
        print("=" * 60)
        for _, row in df_results.iterrows():
            score_str = (
                f"{row['silhouette_score']:.3f}" if not np.isnan(row["silhouette_score"]) else "N/A"
            )
            print(f"  {row['variant']}: {score_str} ({row['n_proteins']} proteins)")

    return df_results
