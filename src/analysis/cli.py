"""CLI commands for ToxProt analysis."""

from pathlib import Path

import click
import pandas as pd


def filter_by_definition(df: pd.DataFrame, definition: str) -> pd.DataFrame:
    """Filter DataFrame by ToxProt definition column.

    Args:
        df: DataFrame with a 'ToxProt definition' column.
        definition: One of 'all', 'venom_tissue', 'kw_toxin', 'both_only'.

    Returns:
        Filtered DataFrame.
    """
    if definition == "all":
        return df
    if "ToxProt definition" not in df.columns:
        click.echo("Warning: 'ToxProt definition' column not found, returning all data.")
        return df

    if definition == "venom_tissue":
        return df[df["ToxProt definition"].isin(["venom_tissue", "both"])]
    elif definition == "kw_toxin":
        return df[df["ToxProt definition"].isin(["kw_toxin", "both"])]
    elif definition == "both_only":
        return df[df["ToxProt definition"] == "both"]
    else:
        return df


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--definition",
    "-d",
    type=click.Choice(["all", "venom_tissue", "kw_toxin", "both_only"]),
    default="venom_tissue",
    show_default=True,
    help="Filter entries by ToxProt definition.",
)
@click.pass_context
def analysis(ctx, definition):
    """Analysis and visualization commands."""
    ctx.ensure_object(dict)
    ctx.obj["definition"] = definition


@analysis.command()
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/toxprot"),
    show_default=True,
    help="Directory containing processed CSV files.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("figures/taxa"),
    show_default=True,
    help="Directory to save output figures.",
)
@click.pass_context
def taxa(ctx, data_dir, output_dir):
    """Run taxonomic distribution analysis.

    Generates trend lines and newcomer analyses for taxa orders across
    ToxProt yearly datasets.

    \b
    Examples:
        toxprot analysis taxa
        toxprot analysis -d all taxa
        toxprot analysis taxa --output-dir figures/taxa_filtered
    """
    from .analyze_taxa import (
        YEARS,
        load_datasets,
        plot_taxa_newcomers,
        plot_top_taxa_trend,
    )

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_datasets(YEARS, data_dir=data_dir)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        # Remove empty datasets after filtering
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets to generate plots")

    plot_top_taxa_trend(datasets, output_dir / "top_taxa_trend")

    # Newcomer plots (comparing 2017 to 2025)
    if "2017" in datasets and "2025" in datasets:
        plot_taxa_newcomers(
            datasets["2017"],
            datasets["2025"],
            "Order",
            output_dir / "taxa_newcomers_order.png",
            "tab:green",
        )
        plot_taxa_newcomers(
            datasets["2017"],
            datasets["2025"],
            "Family",
            output_dir / "taxa_newcomers_family.png",
            "tab:blue",
        )

    click.echo("Taxa analysis complete.")


@analysis.command()
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/toxprot"),
    show_default=True,
    help="Directory containing processed CSV files.",
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
    "--top-n",
    type=int,
    default=15,
    show_default=True,
    help="Number of top protein families to display.",
)
@click.pass_context
def families(ctx, data_dir, output_dir, top_n):
    """Run protein family distribution analysis.

    Generates sequence length histograms, stacked bar charts of protein
    families, and summary statistics tables.

    \b
    Examples:
        toxprot analysis families
        toxprot analysis -d all families
        toxprot analysis families --top-n 20
    """
    from .analyze_protein_families import (
        generate_summary_table,
        plot_sequence_length_histogram,
        plot_stacked_bar_protein_families,
    )

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)
    protein_families_dir = output_dir / "protein_families"
    protein_families_dir.mkdir(parents=True, exist_ok=True)

    # Load datasets
    years = ["2005", "2015", "2025"]
    click.echo("Loading ToxProt datasets...")
    datasets = {}
    for year in years:
        filepath = data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)
            click.echo(f"  Loaded {year}: {len(datasets[year]):,} entries")
        else:
            click.echo(f"  Warning: {filepath} not found, skipping")

    if definition != "all":
        click.echo(f"Filtering by definition: {definition}")
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets")

    # Define output paths
    length_hist_path = output_dir / "sequence_length_distribution.png"
    families_bar_path = protein_families_dir / "top_families_stacked_bar.png"
    summary_table_path = output_dir / "dataset_summary_statistics.png"

    click.echo("\nGenerating visualizations...")

    if "2017" in datasets and "2025" in datasets:
        click.echo("  Sequence length distribution histogram...")
        plot_sequence_length_histogram(datasets["2017"], datasets["2025"], length_hist_path)

        click.echo(f"  Top {top_n} protein families stacked bar chart...")
        plot_stacked_bar_protein_families(
            datasets["2017"], datasets["2025"], families_bar_path, top_n=top_n
        )

    click.echo("  Summary statistics table...")
    generate_summary_table(datasets, summary_table_path)

    click.echo("\nDone!")


@analysis.command()
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/toxprot"),
    show_default=True,
    help="Directory containing processed CSV files.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("figures/ptm"),
    show_default=True,
    help="Directory to save output figures.",
)
@click.option(
    "--years",
    type=str,
    default="2005,2015,2025",
    show_default=True,
    help="Comma-separated years for comparison (e.g., '2005,2015,2025').",
)
@click.pass_context
def ptm(ctx, data_dir, output_dir, years):
    """Run PTM (post-translational modification) analysis.

    Creates a two-panel overview figure comparing PTM frequencies and
    distributions across multiple time points (default: 2005, 2015, 2025).

    \b
    Examples:
        toxprot analysis ptm
        toxprot analysis -d all ptm
        toxprot analysis ptm --years 2010,2020
    """
    from .analyze_ptm import create_combined_figure, load_data

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    year_list = [int(y.strip()) for y in years.split(",")]

    datasets = {}
    for year in year_list:
        filepath = data_dir / f"toxprot_{year}.csv"
        if not filepath.exists():
            raise click.ClickException(f"File not found: {filepath}")
        datasets[year] = load_data(filepath)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets after filtering")

    click.echo(f"Loaded {sum(len(df) for df in datasets.values()):,} entries")
    create_combined_figure(datasets, output_dir)
    click.echo(f"Saved to {output_dir}")


@analysis.command()
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/toxprot"),
    show_default=True,
    help="Directory containing processed CSV files.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=Path("figures/definitions"),
    show_default=True,
    help="Directory to save output figures.",
)
@click.option(
    "--year",
    type=int,
    default=2025,
    show_default=True,
    help="Year of dataset to analyze.",
)
def definitions(data_dir, output_dir, year):
    """Compare ToxProt definition categories.

    Creates a two-panel figure comparing entries by their ToxProt definition:
    - Panel A: Alluvial diagram showing taxonomic coverage
    - Panel B: Venn diagram showing entry overlap

    \b
    Examples:
        toxprot analysis definitions
        toxprot analysis definitions --year 2024
        toxprot analysis definitions -o figures/custom_dir
    """
    from .analyze_definitions import create_definition_comparison_figure

    filepath = data_dir / f"toxprot_{year}.csv"
    if not filepath.exists():
        raise click.ClickException(f"File not found: {filepath}")

    df = pd.read_csv(filepath)
    click.echo(f"Loaded {len(df):,} entries from {year}")

    create_definition_comparison_figure(df, output_dir)
    click.echo(f"Saved to {output_dir}")
