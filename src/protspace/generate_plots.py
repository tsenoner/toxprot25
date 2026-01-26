#!/usr/bin/env python3
"""
Generate protspace visualization plots for 2025 ToxProt dataset.

Creates four plots:
1. All data (including NaN and Other)
2. Top 15 families (full sequences)
3. Top 15 families (mature sequences)
4. Top 15 families (mature sequences, no fragments)

All plots use UMAP embeddings colored by protein families.
"""

from pathlib import Path

from protspace import ProtSpace


def generate_plot(
    styled_json: Path, output_plot: Path, plot_title: str, projection: str = "UMAP_2"
) -> bool:
    """
    Generate a plot from styled JSON file using ProtSpace.

    Args:
        styled_json: Path to styled JSON file
        output_plot: Path for output PNG (without extension)
        plot_title: Title for the plot (used for logging)
        projection: Projection method to use

    Returns:
        True if successful, False otherwise
    """
    print(f"\nGenerating: {plot_title}")
    print("-" * 60)

    if not styled_json.exists():
        print(f"  ✗ Styled JSON not found: {styled_json}")
        return False

    try:
        # Initialize ProtSpace with the JSON file
        protspace = ProtSpace(default_json_file=str(styled_json))

        # Generate plot
        # Remove .png extension if present, as ProtSpace will add it
        output_path_no_ext = output_plot.parent / output_plot.stem

        protspace.generate_plot(
            projection=projection,
            feature="Protein families",
            filename=str(output_path_no_ext),
            width=1600,
            height=1200,
            file_format="png",
        )

        # Check if file was created
        expected_output = Path(f"{output_path_no_ext}.png")
        if expected_output.exists():
            size_kb = expected_output.stat().st_size / 1024
            print(f"  ✓ Plot created: {expected_output.name} ({size_kb:.1f} KB)")
            return True
        else:
            print(f"  ✗ Plot file not found: {expected_output}")
            return False

    except Exception as e:
        print(f"  ✗ Error creating plot: {e}")
        import traceback

        traceback.print_exc()
        return False


def main():
    """Main execution function."""
    base_dir = Path(__file__).resolve().parent.parent.parent
    protspace_dir = base_dir / "data" / "processed" / "protspace"
    figures_dir = base_dir / "figures" / "protspace"

    year = "2025"

    # Define all plots to generate
    plots = [
        {
            "title": "ToxProt 2025 - All Data (Top 15 + Other + NaN)",
            "input": protspace_dir / f"protspace_{year}_all_style.json",
            "output": figures_dir / f"{year}_all_data.png",
        },
        {
            "title": "ToxProt 2025 - Top 15 Families (Full Sequences)",
            "input": protspace_dir / f"protspace_{year}_top15_style.json",
            "output": figures_dir / f"{year}_top15_full_sequences.png",
        },
        {
            "title": "ToxProt 2025 - Top 15 Families (Mature Sequences)",
            "input": protspace_dir / f"protspace_{year}_top15_mature_style.json",
            "output": figures_dir / f"{year}_top15_mature.png",
        },
        {
            "title": "ToxProt 2025 - Top 15 Families (Mature, No Fragments)",
            "input": protspace_dir
            / f"protspace_{year}_top15_mature_no_fragments_style.json",
            "output": figures_dir / f"{year}_top15_mature_no_fragments.png",
        },
    ]

    print("=" * 60)
    print(f"Generating Protspace Plots - {year} Dataset")
    print("=" * 60)

    # Create output directory
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Generate all plots
    success_count = 0
    for plot_config in plots:
        success = generate_plot(
            styled_json=plot_config["input"],
            output_plot=plot_config["output"],
            plot_title=plot_config["title"],
        )
        if success:
            success_count += 1

    # Summary
    print("\n" + "=" * 60)
    print("Plot Generation Summary")
    print("=" * 60)
    print(f"Successfully generated: {success_count}/{len(plots)} plots")

    if success_count == len(plots):
        print(f"✓ All plots saved to: {figures_dir}")
    else:
        print("⚠ Some plots failed to generate")


if __name__ == "__main__":
    main()
