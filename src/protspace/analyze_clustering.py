#!/usr/bin/env python3
"""
Analyze clustering quality using silhouette scores.

Compares clustering quality across three 2025 dataset variants:
1. Top 15 families (full sequences)
2. Top 15 families (mature sequences)
3. Top 15 families (mature sequences, no fragments)

The silhouette score measures how well-separated clusters are.
Score ranges from -1 to 1, where higher is better.
"""

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score


def load_protspace_json(json_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Load UMAP embeddings and protein family labels from protspace JSON.

    Args:
        json_path: Path to protspace JSON file

    Returns:
        Tuple of (embeddings, labels)
    """
    with open(json_path) as f:
        data = json.load(f)

    embeddings = []
    labels = []

    # Get projections and protein data
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

    # Extract data from UMAP projection
    umap_data = umap_projection.get("data", [])

    for item in umap_data:
        protein_id = item.get("identifier")
        coords_dict = item.get("coordinates", {})
        coords = [coords_dict.get("x", 0), coords_dict.get("y", 0)]
        embeddings.append(coords)

        # Get protein family label
        features = protein_data.get(protein_id, {}).get("features", {})
        family = features.get("Protein families", "Unknown")
        labels.append(family)

    return np.array(embeddings), np.array(labels)


def calculate_silhouette_score(embeddings: np.ndarray, labels: np.ndarray) -> float:
    """
    Calculate silhouette score for clustering.

    Args:
        embeddings: 2D array of coordinates
        labels: 1D array of cluster labels

    Returns:
        Silhouette score
    """
    # Filter out invalid labels
    valid_mask = (labels != "Unknown") & (labels != "nan") & (labels != "Other")

    if valid_mask.sum() < 2:
        return np.nan

    embeddings_filtered = embeddings[valid_mask]
    labels_filtered = labels[valid_mask]

    # Need at least 2 unique labels
    unique_labels = np.unique(labels_filtered)
    if len(unique_labels) < 2:
        return np.nan

    try:
        score = silhouette_score(embeddings_filtered, labels_filtered)
        return score
    except Exception as e:
        print(f"    Warning: Could not calculate silhouette score: {e}")
        return np.nan


def analyze_variants(base_dir: Path, year: str = "2025") -> pd.DataFrame:
    """
    Analyze clustering quality for dataset variants.

    Args:
        base_dir: Base directory of the project
        year: Dataset year

    Returns:
        DataFrame with analysis results
    """
    protspace_dir = base_dir / "data" / "processed" / "protspace"

    variants = [
        {
            "name": "Top 15",
            "description": "Top 15 families, full sequences",
            "json": protspace_dir / f"protspace_{year}_top15_style.json",
        },
        {
            "name": "Top 15\nmature",
            "description": "Top 15 families, mature sequences",
            "json": protspace_dir / f"protspace_{year}_top15_mature_style.json",
        },
        {
            "name": "Top 15\nmature\nno fragments",
            "description": "Top 15 families, mature, no fragments",
            "json": protspace_dir
            / f"protspace_{year}_top15_mature_no_fragments_style.json",
        },
    ]

    results = []

    print("=" * 60)
    print(f"Silhouette Score Analysis - {year} Dataset")
    print("=" * 60)

    for variant in variants:
        print(f"\nAnalyzing: {variant['name']}")
        print(f"  {variant['description']}")

        if not variant["json"].exists():
            print(f"  ✗ JSON file not found: {variant['json']}")
            results.append(
                {
                    "Variant": variant["name"],
                    "Description": variant["description"],
                    "N_proteins": 0,
                    "N_families": 0,
                    "Silhouette_score": np.nan,
                }
            )
            continue

        # Load embeddings and labels
        embeddings, labels = load_protspace_json(variant["json"])

        n_proteins = len(labels)
        n_families = len(np.unique(labels))

        # Calculate silhouette score
        score = calculate_silhouette_score(embeddings, labels)

        print(f"  Proteins: {n_proteins}")
        print(f"  Families: {n_families}")
        print(
            f"  Silhouette score: {score:.2f}"
            if not np.isnan(score)
            else "  Silhouette score: N/A"
        )

        results.append(
            {
                "Variant": variant["name"],
                "Description": variant["description"],
                "N_proteins": n_proteins,
                "N_families": n_families,
                "Silhouette_score": score,
            }
        )

    return pd.DataFrame(results)


def plot_silhouette_comparison(df: pd.DataFrame, output_path: Path):
    """
    Create visualization of silhouette scores across variants.

    Args:
        df: DataFrame with analysis results
        output_path: Path to save plot
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Filter out NaN scores and reverse order
    df_plot = df[df["Silhouette_score"].notna()].copy()
    df_plot = df_plot.iloc[::-1]  # Reverse order

    # Create bar plot with single color
    ax.barh(
        df_plot["Variant"],
        df_plot["Silhouette_score"],
        color="#b1b3c6",  # Monocolor blue
        zorder=3,  # Put bars in front of grid
    )

    # Customize plot
    ax.set_xlabel("Silhouette Score", fontsize=18)
    ax.set_title(
        "Clustering Quality Comparison\n(Higher score = better separation)",
        fontsize=18,
        fontweight="bold",
    )
    ax.set_xlim(0, 1)

    # Increase tick label font size
    ax.tick_params(axis="both", which="major", labelsize=16)

    # Add dashed grid behind bars
    ax.grid(axis="x", alpha=0.5, linestyle="--", linewidth=2, zorder=0)

    # Add score labels on bars (rounded to 2 decimal places)
    for i, (_idx, row) in enumerate(df_plot.iterrows()):
        score = row["Silhouette_score"]
        ax.text(score + 0.02, i, f"{score:.2f}", va="center", fontsize=18)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"\n✓ Plot saved to: {output_path}")
    plt.close()


def main():
    """Main execution function."""
    base_dir = Path(__file__).resolve().parent.parent.parent
    figures_dir = base_dir / "figures" / "protspace"

    year = "2025"

    # Analyze variants
    df_results = analyze_variants(base_dir, year)

    # Save results to CSV
    output_csv = figures_dir / "silhouette_comparison.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df_results.to_csv(output_csv, index=False)
    print(f"\n✓ Results saved to: {output_csv}")

    # Create visualization
    output_plot = figures_dir / "silhouette_comparison.png"
    plot_silhouette_comparison(df_results, output_plot)

    # Print summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(df_results.to_string(index=False))
    print("\n✓ Analysis complete!")


if __name__ == "__main__":
    main()
