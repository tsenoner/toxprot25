"""CLI commands for ProtSpace analysis."""

from pathlib import Path

import click
import pandas as pd

from .config import (
    DEFAULT_MIN_DIST,
    DEFAULT_N_NEIGHBORS,
    DEFAULT_YEAR,
    TOP_N,
    VARIANT_CONFIGS,
)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def protspace():
    """ProtSpace embedding and visualization commands.

    Generate UMAP visualizations of protein embeddings colored by protein family.

    \b
    Workflow:
    1. generate-fasta  - Create FASTA files for embedding generation
    2. [External]      - Generate H5 embeddings via Google Colab
    3. prepare         - Prepare metadata and filter H5 files
    4. run-umap        - Run ProtSpace UMAP dimensionality reduction
    5. silhouette      - Analyze clustering quality

    Or run all steps (except Colab) with: pipeline
    """
    pass


@protspace.command("generate-fasta")
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/toxprot"),
    show_default=True,
    help="Directory containing processed CSV files.",
)
@click.option(
    "--interim-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/interim/toxprot_parsed"),
    show_default=True,
    help="Directory containing interim TSV files.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("data/processed/protspace"),
    show_default=True,
    help="Directory to save FASTA files.",
)
@click.option(
    "--year",
    type=str,
    default=DEFAULT_YEAR,
    show_default=True,
    help="Dataset year.",
)
def generate_fasta(
    data_dir: Path,
    interim_dir: Path,
    output_dir: Path,
    year: str,
):
    """Generate FASTA files for ProtT5 embedding generation.

    Creates FASTA files:
    - Full sequences (with signal peptides)
    - Mature sequences (signal peptides removed based on UniProt annotations)

    Only includes proteins with ToxProt definition = "venom_tissue" or "both".

    \b
    Examples:
        toxprot analysis protspace generate-fasta
        toxprot analysis protspace generate-fasta --year 2024
    """
    from .fasta_generator import generate_fasta_files, print_next_steps

    processed_csv = data_dir / f"toxprot_{year}.csv"
    interim_tsv = interim_dir / f"toxprot_{year}.tsv"

    if not processed_csv.exists():
        raise click.ClickException(f"Processed CSV not found: {processed_csv}")
    if not interim_tsv.exists():
        raise click.ClickException(f"Interim TSV not found: {interim_tsv}")

    fasta_full, fasta_mature = generate_fasta_files(
        interim_tsv=interim_tsv,
        output_dir=output_dir,
        year=year,
        processed_csv=processed_csv,
    )

    print_next_steps(fasta_full, fasta_mature, year=year)
    click.echo("\nFASTA generation complete.")


@protspace.command("prepare")
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/toxprot"),
    show_default=True,
    help="Directory containing processed CSV files.",
)
@click.option(
    "--interim-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/interim/toxprot_parsed"),
    show_default=True,
    help="Directory containing interim TSV files.",
)
@click.option(
    "--protspace-dir",
    type=click.Path(path_type=Path),
    default=Path("data/processed/protspace"),
    show_default=True,
    help="Directory for protspace files.",
)
@click.option(
    "--year",
    type=str,
    default=DEFAULT_YEAR,
    show_default=True,
    help="Dataset year.",
)
@click.option(
    "--top-n",
    type=int,
    default=TOP_N,
    show_default=True,
    help="Number of top protein families to track.",
)
def prepare(
    data_dir: Path,
    interim_dir: Path,
    protspace_dir: Path,
    year: str,
    top_n: int,
):
    """Prepare metadata and filter H5 files for all variants.

    Requires H5 embedding files from Google Colab:
    - toxprot_{year}_full.h5
    - toxprot_{year}_mature.h5

    Creates metadata CSV and filtered H5 files for 4 variants:
    - all: All data with labels (top N + Other + NaN)
    - top10_full: Top N families, full sequences
    - top10_mature: Top N families, mature sequences
    - top10_mature_clean: Top N families, mature, no fragments

    Also generates style.json from colors.py for color consistency.

    \b
    Examples:
        toxprot analysis protspace prepare
        toxprot analysis protspace prepare --top-n 15
    """
    from .metadata_preparer import prepare_all_variants
    from .style_generator import generate_style_json

    processed_csv = data_dir / f"toxprot_{year}.csv"
    interim_tsv = interim_dir / f"toxprot_{year}.tsv"

    if not processed_csv.exists():
        raise click.ClickException(f"Processed CSV not found: {processed_csv}")
    if not interim_tsv.exists():
        raise click.ClickException(f"Interim TSV not found: {interim_tsv}")

    # Load reference data for style generation
    df_ref = pd.read_csv(processed_csv)
    if "ToxProt definition" in df_ref.columns:
        df_ref = df_ref[df_ref["ToxProt definition"].isin(["venom_tissue", "both"])]

    # Generate style.json
    style_path = protspace_dir / "style.json"
    click.echo("Generating style.json...")
    generate_style_json(df_ref, style_path, top_n=top_n)

    # Prepare all variants
    try:
        prepare_all_variants(
            processed_csv=processed_csv,
            interim_tsv=interim_tsv,
            protspace_dir=protspace_dir,
            year=year,
            top_n=top_n,
        )
    except FileNotFoundError as e:
        raise click.ClickException(str(e)) from e

    click.echo("\nPreparation complete.")


@protspace.command("run-umap")
@click.option(
    "--protspace-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/protspace"),
    show_default=True,
    help="Directory containing H5 and metadata files.",
)
@click.option(
    "--year",
    type=str,
    default=DEFAULT_YEAR,
    show_default=True,
    help="Dataset year.",
)
@click.option(
    "--n-neighbors",
    type=int,
    default=DEFAULT_N_NEIGHBORS,
    show_default=True,
    help="UMAP n_neighbors parameter.",
)
@click.option(
    "--min-dist",
    type=float,
    default=DEFAULT_MIN_DIST,
    show_default=True,
    help="UMAP min_dist parameter.",
)
@click.option(
    "--variant",
    type=click.Choice(list(VARIANT_CONFIGS.keys())),
    multiple=True,
    help="Specific variant(s) to process. Default: all variants.",
)
def run_umap(
    protspace_dir: Path,
    year: str,
    n_neighbors: int,
    min_dist: float,
    variant: tuple[str, ...],
):
    """Run ProtSpace UMAP dimensionality reduction.

    Generates UMAP projections and applies styling from style.json.

    \b
    Examples:
        toxprot analysis protspace run-umap
        toxprot analysis protspace run-umap --n-neighbors 30 --min-dist 0.3
        toxprot analysis protspace run-umap --variant top10_mature
    """
    from .umap_runner import run_umap_all_variants

    style_file = protspace_dir / "style.json"
    if not style_file.exists():
        raise click.ClickException(
            f"Style file not found: {style_file}\nRun 'toxprot analysis protspace prepare' first."
        )

    variants = list(variant) if variant else None

    results = run_umap_all_variants(
        protspace_dir=protspace_dir,
        style_file=style_file,
        year=year,
        n_neighbors=n_neighbors,
        min_dist=min_dist,
        variants=variants,
    )

    failed = [v for v, success in results.items() if not success]
    if failed:
        click.echo(f"\nWarning: Failed variants: {', '.join(failed)}")


@protspace.command("silhouette")
@click.option(
    "--protspace-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/protspace"),
    show_default=True,
    help="Directory containing styled JSON files.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("figures"),
    show_default=True,
    help="Directory to save analysis outputs.",
)
@click.option(
    "--year",
    type=str,
    default=DEFAULT_YEAR,
    show_default=True,
    help="Dataset year.",
)
@click.option(
    "--variant",
    type=click.Choice(list(VARIANT_CONFIGS.keys())),
    multiple=True,
    help="Specific variant(s) to analyze. Default: all variants.",
)
def silhouette(
    protspace_dir: Path,
    output_dir: Path,
    year: str,
    variant: tuple[str, ...],
):
    """Analyze clustering quality using silhouette scores.

    Calculates silhouette scores for each variant and generates
    a comparison plot. Higher scores indicate better cluster separation.

    \b
    Examples:
        toxprot analysis protspace silhouette
        toxprot analysis protspace silhouette --variant top10_mature
    """
    from .clustering_analyzer import run_silhouette_analysis

    variants = list(variant) if variant else None

    run_silhouette_analysis(
        protspace_dir=protspace_dir,
        figures_dir=output_dir,
        year=year,
        variants=variants,
    )

    click.echo("\nSilhouette analysis complete.")


@protspace.command("pipeline")
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/toxprot"),
    show_default=True,
    help="Directory containing processed CSV files.",
)
@click.option(
    "--interim-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/interim/toxprot_parsed"),
    show_default=True,
    help="Directory containing interim TSV files.",
)
@click.option(
    "--protspace-dir",
    type=click.Path(path_type=Path),
    default=Path("data/processed/protspace"),
    show_default=True,
    help="Directory for protspace files.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("figures"),
    show_default=True,
    help="Directory to save output figures.",
)
@click.option(
    "--year",
    type=str,
    default=DEFAULT_YEAR,
    show_default=True,
    help="Dataset year.",
)
@click.option(
    "--top-n",
    type=int,
    default=TOP_N,
    show_default=True,
    help="Number of top protein families to track.",
)
@click.option(
    "--skip-fasta",
    is_flag=True,
    help="Skip FASTA generation (assume files exist).",
)
@click.pass_context
def pipeline(
    ctx,
    data_dir: Path,
    interim_dir: Path,
    protspace_dir: Path,
    output_dir: Path,
    year: str,
    top_n: int,
    skip_fasta: bool,
):
    """Run full ProtSpace analysis pipeline.

    Runs all steps except the external Colab embedding generation:
    1. Generate FASTA files (unless --skip-fasta)
    2. Prepare metadata and filter H5 files
    3. Run UMAP dimensionality reduction
    4. Silhouette analysis

    Note: H5 embedding files must already exist from Google Colab.

    \b
    Examples:
        toxprot analysis protspace pipeline
        toxprot analysis protspace pipeline --skip-fasta
        toxprot analysis protspace pipeline --top-n 15
    """
    click.echo("=" * 60)
    click.echo("ProtSpace Analysis Pipeline")
    click.echo("=" * 60)

    # Step 1: Generate FASTA
    if not skip_fasta:
        click.echo("\n[1/4] Generating FASTA files...")
        ctx.invoke(
            generate_fasta,
            data_dir=data_dir,
            interim_dir=interim_dir,
            output_dir=protspace_dir,
            year=year,
        )
    else:
        click.echo("\n[1/4] Skipping FASTA generation")

    # Step 2: Prepare metadata and H5
    click.echo("\n[2/4] Preparing metadata and H5 variants...")
    ctx.invoke(
        prepare,
        data_dir=data_dir,
        interim_dir=interim_dir,
        protspace_dir=protspace_dir,
        year=year,
        top_n=top_n,
    )

    # Step 3: Run UMAP
    click.echo("\n[3/4] Running UMAP...")
    ctx.invoke(
        run_umap,
        protspace_dir=protspace_dir,
        year=year,
    )

    # Step 4: Silhouette analysis
    click.echo("\n[4/4] Running silhouette analysis...")
    ctx.invoke(
        silhouette,
        protspace_dir=protspace_dir,
        output_dir=output_dir,
        year=year,
    )

    click.echo("\n" + "=" * 60)
    click.echo("Pipeline complete!")
    click.echo("=" * 60)
    click.echo(f"Figures saved to: {output_dir}")
