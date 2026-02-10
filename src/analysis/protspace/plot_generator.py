"""Generate ProtSpace visualization plots.

Creates UMAP scatter plots colored by protein families for each variant.
Supports both parquetbundle (v3.1.1+) and legacy JSON formats.
"""

from pathlib import Path

from .config import (
    VARIANT_CONFIGS,
    get_plot_filename,
    get_protspace_styled_filename,
)


def generate_plot(
    styled_bundle: Path,
    output_path: Path,
    annotation: str = "Protein families",
    projection: str = "UMAP_2",
    width: int = 1600,
    height: int = 1200,
    verbose: bool = True,
) -> bool:
    """Generate a plot from styled parquetbundle using ProtSpace.

    Args:
        styled_bundle: Path to styled parquetbundle directory
        output_path: Path for output PNG
        annotation: Annotation to color by
        projection: Projection method name
        width: Plot width in pixels
        height: Plot height in pixels
        verbose: Print progress messages

    Returns:
        True if successful, False otherwise
    """
    if not styled_bundle.exists():
        if verbose:
            print(f"  Styled parquetbundle not found: {styled_bundle}")
        return False

    try:
        # Import here to avoid dependency issues when module loads
        from protspace import ProtSpace

        # Use arrow_dir parameter for parquetbundle format (v3.1.1+)
        protspace = ProtSpace(arrow_dir=str(styled_bundle))

        # ProtSpace adds .png extension, so remove if present
        output_stem = output_path.parent / output_path.stem

        protspace.generate_plot(
            projection=projection,
            annotation=annotation,
            filename=str(output_stem),
            width=width,
            height=height,
            file_format="png",
        )

        # Check if file was created
        expected_output = Path(f"{output_stem}.png")
        if expected_output.exists():
            size_kb = expected_output.stat().st_size / 1024
            if verbose:
                print(f"  Plot created: {expected_output.name} ({size_kb:.1f} KB)")
            return True
        else:
            if verbose:
                print(f"  Plot file not found: {expected_output}")
            return False

    except ImportError:
        if verbose:
            print("  Error: protspace package not installed")
            print("  Install with: pip install protspace")
        return False
    except Exception as e:
        if verbose:
            print(f"  Error creating plot: {e}")
        return False


def generate_all_plots(
    protspace_dir: Path,
    figures_dir: Path,
    year: str = "2025",
    variants: list[str] | None = None,
    verbose: bool = True,
) -> dict[str, bool]:
    """Generate plots for all specified variants.

    Args:
        protspace_dir: Directory containing styled parquetbundle files
        figures_dir: Directory to save plots
        year: Dataset year
        variants: List of variant names to process (None = all)
        verbose: Print progress messages

    Returns:
        Dictionary mapping variant names to success status
    """
    if variants is None:
        variants = list(VARIANT_CONFIGS.keys())

    figures_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("=" * 60)
        print(f"Generating ProtSpace Plots - {year}")
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

        styled_bundle = protspace_dir / get_protspace_styled_filename(year, variant_name)
        output_path = figures_dir / get_plot_filename(year, variant_name)

        success = generate_plot(
            styled_bundle=styled_bundle,
            output_path=output_path,
            verbose=verbose,
        )

        results[variant_name] = success

    # Summary
    if verbose:
        success_count = sum(results.values())
        print("\n" + "=" * 60)
        print(f"Plot Generation: {success_count}/{len(results)} plots created")
        print("=" * 60)
        if success_count > 0:
            print(f"Plots saved to: {figures_dir}")

    return results
