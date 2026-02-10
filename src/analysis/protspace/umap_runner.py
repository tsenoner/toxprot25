"""Run ProtSpace UMAP dimensionality reduction.

Executes protspace-local and applies annotation styling by manually appending
a settings parquet to the parquetbundle file.
"""

import io
import json
import re
import subprocess
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .config import (
    DEFAULT_MIN_DIST,
    DEFAULT_N_NEIGHBORS,
    VARIANT_CONFIGS,
    get_h5_variant_filename,
    get_metadata_filename,
    get_protspace_output_filename,
    get_protspace_styled_filename,
)

# Delimiter used in parquetbundle files to separate parquet parts
PARQUET_BUNDLE_DELIMITER = b"---PARQUET_DELIMITER---"


def run_command(cmd: str, verbose: bool = True) -> tuple[bool, str]:
    """Run a shell command and return success status.

    Args:
        cmd: Command to execute
        verbose: Print command being run

    Returns:
        Tuple of (success, error_message)
    """
    if verbose:
        print(f"  Running: {cmd}")

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        return False, result.stderr
    return True, ""


def process_variant(
    h5_path: Path,
    metadata_path: Path,
    output_bundle: Path,
    styled_bundle: Path,
    style_file: Path,
    n_neighbors: int = DEFAULT_N_NEIGHBORS,
    min_dist: float = DEFAULT_MIN_DIST,
    verbose: bool = True,
) -> bool:
    """Process a single variant: generate UMAP parquetbundle and apply styling.

    Args:
        h5_path: Path to H5 embedding file
        metadata_path: Path to metadata CSV
        output_bundle: Path for intermediate output parquetbundle
        styled_bundle: Path for final styled output parquetbundle
        style_file: Path to style.json
        n_neighbors: UMAP n_neighbors parameter
        min_dist: UMAP min_dist parameter
        verbose: Print progress messages

    Returns:
        True if successful, False otherwise
    """
    # Check input files
    if not h5_path.exists():
        if verbose:
            print(f"  H5 file not found: {h5_path}")
        return False

    if not metadata_path.exists():
        if verbose:
            print(f"  Metadata file not found: {metadata_path}")
        return False

    if not style_file.exists():
        if verbose:
            print(f"  Style file not found: {style_file}")
        return False

    if verbose:
        print("  Input files found")

    # Generate protspace parquetbundle with UMAP (new v3.1.1 syntax)
    cmd_bundle = (
        f"protspace-local -i {h5_path} -a {metadata_path} -o {output_bundle} "
        f"-m umap2 --n_neighbors {n_neighbors} --min_dist {min_dist}"
    )

    success, error = run_command(cmd_bundle, verbose=verbose)
    if not success:
        if verbose:
            print(f"  Error generating parquetbundle: {error}")
        return False

    if not output_bundle.exists():
        if verbose:
            print(f"  Failed to create parquetbundle: {output_bundle}")
        return False

    if verbose:
        print(f"  Parquetbundle created: {output_bundle.name}")

    # Apply styling by appending settings parquet
    success = apply_annotation_styles(output_bundle, style_file, styled_bundle, verbose=verbose)
    if not success:
        return False

    # Remove intermediate parquetbundle
    try:
        output_bundle.unlink()
        if verbose:
            print("  Removed intermediate parquetbundle")
    except Exception:
        pass

    return True


def _rgba_to_hex(rgba_str: str) -> str:
    """Convert 'rgba(r, g, b, a)' to '#RRGGBB'."""
    match = re.match(r"rgba?\((\d+),\s*(\d+),\s*(\d+)", rgba_str)
    if match:
        r, g, b = int(match.group(1)), int(match.group(2)), int(match.group(3))
        return f"#{r:02x}{g:02x}{b:02x}"
    return rgba_str


def _create_settings_parquet(
    style_json: dict,
    annotation_name: str = "Protein families",
) -> bytes:
    """Create a settings parquet from style.json format.

    Args:
        style_json: Style configuration with rgba colors
        annotation_name: Name of the annotation to style

    Returns:
        Bytes of the settings parquet file
    """
    colors = style_json.get(annotation_name, {}).get("colors", {})

    categories = {}
    for i, (name, rgba_color) in enumerate(colors.items()):
        # Skip categories not present in top10 variants
        if name not in ("Other", "NaN"):
            categories[name] = {
                "zOrder": i,
                "color": _rgba_to_hex(rgba_color),
                "shape": "circle",
            }

    settings = {
        annotation_name: {
            "maxVisibleValues": 10,
            "includeShapes": False,
            "shapeSize": 30,
            "sortMode": "size-desc",
            "hiddenValues": [],
            "categories": categories,
            "enableDuplicateStackUI": False,
            "selectedPaletteId": "custom",
        }
    }

    settings_json = json.dumps(settings)
    table = pa.table({"settings_json": [settings_json]})

    buf = io.BytesIO()
    pq.write_table(table, buf)
    return buf.getvalue()


def apply_annotation_styles(
    input_bundle: Path,
    style_file: Path,
    output_bundle: Path,
    verbose: bool = True,
) -> bool:
    """Apply annotation styles (colors, shapes) to a parquetbundle.

    Manually appends a settings parquet to the parquetbundle file,
    as the protspace API has compatibility issues with single-file bundles.

    Args:
        input_bundle: Path to input parquetbundle file
        style_file: Path to style.json with color/shape definitions
        output_bundle: Path for output styled parquetbundle
        verbose: Print progress messages

    Returns:
        True if successful, False otherwise
    """
    try:
        # Load style configuration
        with open(style_file) as f:
            style_json = json.load(f)

        if verbose:
            print(f"  Applying styles from {style_file.name}")

        # Read existing bundle
        with open(input_bundle, "rb") as f:
            content = f.read()

        # Split by delimiter
        parts = content.split(PARQUET_BUNDLE_DELIMITER)

        # Create settings parquet
        settings_parquet = _create_settings_parquet(style_json)

        # Append or replace settings (4th part)
        if len(parts) >= 4:
            parts[3] = settings_parquet
            new_content = PARQUET_BUNDLE_DELIMITER.join(parts)
        else:
            new_content = content + PARQUET_BUNDLE_DELIMITER + settings_parquet

        # Write output
        with open(output_bundle, "wb") as f:
            f.write(new_content)

        if verbose:
            print(f"  Styled parquetbundle created: {output_bundle.name}")

        return True

    except Exception as e:
        if verbose:
            print(f"  Error applying styles: {e}")
        return False


def run_umap_all_variants(
    protspace_dir: Path,
    style_file: Path,
    year: str = "2025",
    n_neighbors: int = DEFAULT_N_NEIGHBORS,
    min_dist: float = DEFAULT_MIN_DIST,
    variants: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, bool]:
    """Run UMAP for all specified variants.

    Args:
        protspace_dir: Directory containing H5 and metadata files
        style_file: Path to style.json
        year: Dataset year
        n_neighbors: UMAP n_neighbors parameter
        min_dist: UMAP min_dist parameter
        variants: List of variant names to process (None = all)
        verbose: Print progress messages

    Returns:
        Dictionary mapping variant names to success status
    """
    if variants is None:
        variants = list(VARIANT_CONFIGS.keys())

    if verbose:
        print("=" * 60)
        print(f"Running UMAP (n_neighbors={n_neighbors}, min_dist={min_dist})")
        print("=" * 60)

    results = {}
    for variant_name in variants:
        if variant_name not in VARIANT_CONFIGS:
            if verbose:
                print(f"\nUnknown variant: {variant_name}")
            results[variant_name] = False
            continue

        config = VARIANT_CONFIGS[variant_name]
        if verbose:
            print(f"\n{config['description']}:")

        h5_path = protspace_dir / get_h5_variant_filename(year, variant_name)
        metadata_path = protspace_dir / get_metadata_filename(year, variant_name)
        output_bundle = protspace_dir / get_protspace_output_filename(year, variant_name)
        styled_bundle = protspace_dir / get_protspace_styled_filename(year, variant_name)

        success = process_variant(
            h5_path=h5_path,
            metadata_path=metadata_path,
            output_bundle=output_bundle,
            styled_bundle=styled_bundle,
            style_file=style_file,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            verbose=verbose,
        )

        results[variant_name] = success

    # Summary
    if verbose:
        success_count = sum(results.values())
        print("\n" + "=" * 60)
        print(f"UMAP Complete: {success_count}/{len(results)} variants successful")
        print("=" * 60)

    return results
