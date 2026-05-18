"""Run ProtSpace projections via the v4 unified CLI.

This module shells out to `protspace prepare` to build the parquetbundle and,
for analysis variants, to `protspace style` to embed the curated palette.
Demo variants ship without styling.
"""

import io
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from .config import (
    DEFAULT_MIN_DIST,
    DEFAULT_N_NEIGHBORS,
    DEMO_SUBDIR,
    INTERMEDIATES_SUBDIR,
    PROTSPACE_MODEL_NAME,
    VARIANT_CONFIGS,
    get_h5_variant_filename,
    get_metadata_filename,
    get_protspace_styled_filename,
)

# protspace `prepare` always writes the bundle as data.parquetbundle.
PROTSPACE_BUNDLE_NAME = "data.parquetbundle"

# Delimiter that separates the parquet parts inside a .parquetbundle.
PARQUET_BUNDLE_DELIMITER = b"---PARQUET_DELIMITER---"

# Annotations the web UI hides from the dropdown (tooltip-only). The first
# non-tooltip annotation column in Part 0 becomes the default selection.
TOOLTIP_ONLY_ANNOTATIONS = ("protein_name", "uniprot_kb_id", "gene_name")


def _reorder_bundle_annotations(bundle_path: Path, primary_annotation: str) -> None:
    """Move ``primary_annotation`` to the very front of the annotation parquet.

    The control-bar dropdown filters out tooltip-only fields and then defaults
    to its ``annotations[0]``, but the scatter-plot component reads the raw
    ``Object.keys(data.annotations)`` (unfiltered) for its initial coloring,
    so a tooltip-only column that happens to come first ends up driving the
    plot. Placing ``primary_annotation`` ahead of the tooltip-only fields
    makes both code paths agree on the same default.
    """
    with open(bundle_path, "rb") as f:
        parts = f.read().split(PARQUET_BUNDLE_DELIMITER)
    if not parts:
        return

    table = pq.read_table(io.BytesIO(parts[0]))
    columns = list(table.column_names)
    if primary_annotation not in columns:
        return

    id_columns = [c for c in columns if c == "protein_id"]
    tooltip = [c for c in columns if c in TOOLTIP_ONLY_ANNOTATIONS and c not in id_columns]
    rest = [
        c
        for c in columns
        if c != primary_annotation and c not in id_columns and c not in tooltip
    ]
    new_order = id_columns + [primary_annotation] + tooltip + rest
    new_table = table.select(new_order)

    buf = io.BytesIO()
    pq.write_table(new_table, buf)
    parts[0] = buf.getvalue()

    with open(bundle_path, "wb") as f:
        f.write(PARQUET_BUNDLE_DELIMITER.join(parts))


def _filter_style_to_variant(
    style_file: Path,
    metadata_path: Path,
) -> dict:
    """Build a variant-specific style dict by intersecting with metadata values.

    `protspace style` (v4) rejects styles that reference values not present in
    the data. Each variant has a different subset of categorical values (e.g.
    `mature_clean` removes fragments so `has_fragment` has only `"no"`), so we
    must trim the shared style.json down to values that actually appear.
    """
    with open(style_file) as f:
        style = json.load(f)
    df = pd.read_csv(metadata_path)

    filtered: dict = {}
    for column, cfg in style.items():
        if column not in df.columns:
            continue
        present = {str(v) for v in df[column].dropna().unique()}
        if df[column].isna().any():
            present.add("nan")
        new_cfg = dict(cfg)
        if "colors" in new_cfg:
            new_cfg["colors"] = {k: v for k, v in new_cfg["colors"].items() if k in present}
            if not new_cfg["colors"]:
                continue
        filtered[column] = new_cfg
    return filtered


def _run(cmd: list[str], verbose: bool = True) -> tuple[bool, str]:
    """Run a shell command via subprocess and return (success, stderr)."""
    if verbose:
        print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return False, result.stderr or result.stdout
    return True, ""


def _output_bundle_path(
    variant_config: dict,
    protspace_dir: Path,
    year: str,
) -> Path:
    """Return the final on-disk parquetbundle path for a variant."""
    filename = get_protspace_styled_filename(year, variant_config["name"])
    if variant_config["bundle_kind"] == "demo":
        return protspace_dir / DEMO_SUBDIR / filename
    return protspace_dir / filename


def process_variant(
    variant_config: dict,
    h5_path: Path,
    metadata_path: Path,
    final_bundle: Path,
    style_file: Path | None,
    n_neighbors: int = DEFAULT_N_NEIGHBORS,
    min_dist: float = DEFAULT_MIN_DIST,
    verbose: bool = True,
) -> bool:
    """Build one variant: prepare → (optionally style) → write to final path.

    Args:
        variant_config: Config dict from ``VARIANT_CONFIGS``.
        h5_path: Filtered H5 embedding file for this variant.
        metadata_path: Annotation CSV for this variant.
        final_bundle: Final parquetbundle output path.
        style_file: Path to style JSON. Required for analysis variants;
            ignored for demo variants (which ship unstyled).
        n_neighbors: UMAP n_neighbors.
        min_dist: UMAP min_dist.
        verbose: Print progress messages.
    """
    if not h5_path.exists():
        if verbose:
            print(f"  H5 file not found: {h5_path}")
        return False
    if not metadata_path.exists():
        if verbose:
            print(f"  Metadata file not found: {metadata_path}")
        return False

    final_bundle.parent.mkdir(parents=True, exist_ok=True)
    bundle_kind = variant_config["bundle_kind"]

    # Demo bundles use a persistent output directory so the protspace `tmp/`
    # cache (UniProt fetches) is reused across re-runs. Analysis bundles use a
    # throw-away temp dir because they don't fetch UniProt annotations.
    if bundle_kind == "demo":
        work_dir_ctx = None
        work_dir = final_bundle.parent / f".{variant_config['name']}_workdir"
        work_dir.mkdir(parents=True, exist_ok=True)
    else:
        work_dir_ctx = tempfile.TemporaryDirectory(prefix="protspace_")
        work_dir = Path(work_dir_ctx.name)

    try:
        cmd_prepare = [
            "protspace",
            "prepare",
            "-i",
            f"{h5_path}:{PROTSPACE_MODEL_NAME}",
            "-a",
            str(metadata_path),
        ]
        if bundle_kind == "demo":
            # Fetch the UniProt annotations the web UI needs for its default
            # entry chip (protein name, mnemonic, etc.). `reviewed` is omitted
            # because every ToxProt entry is Swiss-Prot by construction.
            for ann in ("protein_name", "uniprot_kb_id", "gene_name", "keyword", "ec"):
                cmd_prepare += ["-a", ann]
        cmd_prepare += [
            "-o",
            str(work_dir),
            "-m",
            f"umap2:n_neighbors={n_neighbors};min_dist={min_dist},pca2",
            "--bundled",
            "--no-log",
        ]
        # Keep the protspace tmp/ cache for demo runs to avoid re-fetching
        # UniProt annotations every time.
        cmd_prepare += ["--keep-tmp"] if bundle_kind == "demo" else ["--no-keep-tmp"]
        success, error = _run(cmd_prepare, verbose=verbose)
        if not success:
            if verbose:
                print(f"  Error running `protspace prepare`: {error}")
            return False

        prepared_bundle = work_dir / PROTSPACE_BUNDLE_NAME
        if not prepared_bundle.exists():
            if verbose:
                print(f"  Expected bundle not found: {prepared_bundle}")
            return False

        if bundle_kind == "analysis":
            if style_file is None or not style_file.exists():
                if verbose:
                    print(f"  Style file required for analysis variant: {style_file}")
                return False
            variant_style = _filter_style_to_variant(style_file, metadata_path)
            variant_style_file = work_dir / "style.json"
            with open(variant_style_file, "w") as f:
                json.dump(variant_style, f, indent=2)

            cmd_style = [
                "protspace",
                "style",
                str(prepared_bundle),
                str(final_bundle),
                "--annotation-styles",
                str(variant_style_file),
            ]
            success, error = _run(cmd_style, verbose=verbose)
            if not success:
                if verbose:
                    print(f"  Error running `protspace style`: {error}")
                return False
        else:
            # Demo bundles ship unstyled.
            shutil.copy2(prepared_bundle, final_bundle)
            # Make `protein_families` the default colouring in the web UI.
            _reorder_bundle_annotations(final_bundle, "protein_families")
    finally:
        if work_dir_ctx is not None:
            work_dir_ctx.cleanup()

    if verbose:
        print(f"  Wrote: {final_bundle}")
    return True


def run_umap_all_variants(
    protspace_dir: Path,
    style_file: Path,
    year: str = "2025",
    n_neighbors: int = DEFAULT_N_NEIGHBORS,
    min_dist: float = DEFAULT_MIN_DIST,
    variants: list[str] | None = None,
    verbose: bool = True,
    cleanup_intermediates: bool = True,
) -> dict[str, bool]:
    """Build the requested variants. Returns ``{variant: success}``."""
    if variants is None:
        variants = list(VARIANT_CONFIGS.keys())

    if verbose:
        print("=" * 60)
        print(f"Running UMAP+PCA (n_neighbors={n_neighbors}, min_dist={min_dist})")
        print("=" * 60)

    intermediates_dir = protspace_dir / INTERMEDIATES_SUBDIR

    results: dict[str, bool] = {}
    for variant_name in variants:
        if variant_name not in VARIANT_CONFIGS:
            if verbose:
                print(f"\nUnknown variant: {variant_name}")
            results[variant_name] = False
            continue

        config = VARIANT_CONFIGS[variant_name]
        if verbose:
            print(f"\n{config['description']}:")

        h5_path = intermediates_dir / get_h5_variant_filename(year, variant_name)
        metadata_path = intermediates_dir / get_metadata_filename(year, variant_name)
        final_bundle = _output_bundle_path(config, protspace_dir, year)

        results[variant_name] = process_variant(
            variant_config=config,
            h5_path=h5_path,
            metadata_path=metadata_path,
            final_bundle=final_bundle,
            style_file=style_file if config["bundle_kind"] == "analysis" else None,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            verbose=verbose,
        )

    if cleanup_intermediates and intermediates_dir.exists():
        shutil.rmtree(intermediates_dir)
        if verbose:
            print(f"\n  Cleaned up {INTERMEDIATES_SUBDIR}/ directory")

    if verbose:
        success_count = sum(results.values())
        print("\n" + "=" * 60)
        print(f"Complete: {success_count}/{len(results)} variants successful")
        print("=" * 60)

    return results
