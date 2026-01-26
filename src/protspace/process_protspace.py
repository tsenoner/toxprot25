#!/usr/bin/env python3
"""
Process protspace for ToxProt datasets.

Generates JSON and styled files for all 2025 data variants:
1. All data (including NaN and Other)
2. Top 15 families (full sequences)
3. Top 15 families (mature sequences)
4. Top 15 families (mature sequences, no fragments)
"""

import shutil
import subprocess
from pathlib import Path


def run_command(cmd: str) -> bool:
    """
    Run a shell command and return success status.

    Args:
        cmd: Command to execute

    Returns:
        True if successful, False otherwise
    """
    print(f"  Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ✗ Error: {result.stderr}")
        return False
    return True


def process_variant(
    h5_path: Path,
    metadata_path: Path,
    output_json: Path,
    styled_json: Path,
    style_file: Path,
    methods: str = "umap2",
    n_neighbors: int = 50,
    min_dist: float = 0.5,
    variant_name: str = "",
) -> bool:
    """
    Process a single variant: generate JSON and apply styling.

    Args:
        h5_path: Path to H5 file
        metadata_path: Path to metadata CSV
        output_json: Path for output JSON
        styled_json: Path for styled output JSON
        style_file: Path to style JSON
        methods: Methods for dimensionality reduction
        n_neighbors: UMAP n_neighbors parameter
        min_dist: UMAP min_dist parameter
        variant_name: Name of variant for logging

    Returns:
        True if successful, False otherwise
    """
    print(f"\nProcessing: {variant_name}")
    print("-" * 60)

    # Check if input files exist
    if not h5_path.exists():
        print(f"  ✗ H5 file not found: {h5_path}")
        return False
    if not metadata_path.exists():
        print(f"  ✗ Metadata file not found: {metadata_path}")
        return False

    print("  ✓ Input files found")

    # Generate protspace JSON
    cmd_json = (
        f"protspace-local -i {h5_path} -f {metadata_path} -o {output_json} "
        f"--non-binary --bundled false "
        f"-m {methods} --n_neighbors {n_neighbors} --min_dist {min_dist}"
    )

    if not run_command(cmd_json):
        return False

    # Workaround for protspace bug: always check for bug-path output and move
    bug_subdir = output_json.parent / output_json.stem  # e.g. protspace_2025_all/
    bug_file = bug_subdir / "selected_features_projections.json"

    if bug_file.exists():
        shutil.move(str(bug_file), str(output_json))
        # Remove the bug subdirectory after moving file, if possible
        try:
            shutil.rmtree(str(bug_subdir))
        except Exception:
            pass
        print(f"  ✓ JSON moved from protspace bug path: {output_json.name}")

    if not output_json.exists():
        print(f"  ✗ Failed to create JSON: {output_json}")
        return False

    size_mb = output_json.stat().st_size / (1024 * 1024)
    print(f"  ✓ JSON created: {output_json.name} ({size_mb:.1f} MB)")

    # Apply styling
    if not style_file.exists():
        print(f"  ✗ Style file not found: {style_file}")
        return False

    cmd_style = (
        f"protspace-feature-colors --feature_styles {style_file} "
        f"{output_json} {styled_json}"
    )

    if not run_command(cmd_style):
        return False

    if not styled_json.exists():
        print(f"  ✗ Failed to create styled JSON: {styled_json}")
        return False

    # Remove the original output_json file after styling
    try:
        output_json.unlink()
        print(f"  ✓ Removed intermediate JSON: {output_json.name}")
    except Exception as e:
        print(f"  ⚠ Could not remove intermediate JSON: {output_json.name}: {e}")

    size_mb = styled_json.stat().st_size / (1024 * 1024)
    print(f"  ✓ Styled JSON created: {styled_json.name} ({size_mb:.1f} MB)")

    return True


def main():
    """Main execution function."""
    base_dir = Path(__file__).resolve().parent.parent.parent
    protspace_dir = base_dir / "data" / "processed" / "protspace"

    # Configuration
    year = "2025"
    methods = "umap2"
    n_neighbors = 50
    min_dist = 0.5
    style_file = protspace_dir / "style.json"

    # Define all variants to process
    variants = [
        {
            "name": "All data (including NaN and Other)",
            "h5": protspace_dir / f"toxprot_{year}_all.h5",
            "metadata": protspace_dir / f"metadata_{year}_all.csv",
            "output": protspace_dir / f"protspace_{year}_all.json",
            "styled": protspace_dir / f"protspace_{year}_all_style.json",
        },
        {
            "name": "Top 15 families (full sequences)",
            "h5": protspace_dir / f"toxprot_{year}_top15.h5",
            "metadata": protspace_dir / f"metadata_{year}_top15.csv",
            "output": protspace_dir / f"protspace_{year}_top15.json",
            "styled": protspace_dir / f"protspace_{year}_top15_style.json",
        },
        {
            "name": "Top 15 families (mature sequences)",
            "h5": protspace_dir / f"toxprot_{year}_top15_mature.h5",
            "metadata": protspace_dir / f"metadata_{year}_top15_mature.csv",
            "output": protspace_dir / f"protspace_{year}_top15_mature.json",
            "styled": protspace_dir / f"protspace_{year}_top15_mature_style.json",
        },
        {
            "name": "Top 15 families (mature, no fragments)",
            "h5": protspace_dir / f"toxprot_{year}_top15_mature_no_fragments.h5",
            "metadata": protspace_dir
            / f"metadata_{year}_top15_mature_no_fragments.csv",
            "output": protspace_dir
            / f"protspace_{year}_top15_mature_no_fragments.json",
            "styled": protspace_dir
            / f"protspace_{year}_top15_mature_no_fragments_style.json",
        },
    ]

    print("=" * 60)
    print(f"Protspace Processing - {year} Dataset")
    print("=" * 60)

    # Process all variants
    success_count = 0
    for variant in variants:
        success = process_variant(
            h5_path=variant["h5"],
            metadata_path=variant["metadata"],
            output_json=variant["output"],
            styled_json=variant["styled"],
            style_file=style_file,
            methods=methods,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            variant_name=variant["name"],
        )
        if success:
            success_count += 1

    # Summary
    print("\n" + "=" * 60)
    print("Processing Summary")
    print("=" * 60)
    print(f"Successfully processed: {success_count}/{len(variants)} variants")

    if success_count == len(variants):
        print("✓ All variants processed successfully!")
    else:
        print("⚠ Some variants failed to process")


if __name__ == "__main__":
    main()
