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
@click.option(
    "--level",
    "-l",
    type=click.Choice(["Phylum", "Class", "Order", "Family", "all"]),
    default="Order",
    show_default=True,
    help="Taxonomic level for alluvial plot.",
)
@click.option(
    "--skip-trend",
    is_flag=True,
    help="Skip generating the top taxa trend plot.",
)
@click.pass_context
def taxa(ctx, data_dir, output_dir, level, skip_trend):
    """Run taxonomic distribution analysis.

    Generates alluvial diagrams showing taxa evolution across decade steps
    (2005, 2015, 2025) at different taxonomic levels.

    \b
    Examples:
        toxprot analysis taxa
        toxprot analysis taxa --level Family
        toxprot analysis taxa -l Order --skip-trend
        toxprot analysis -d all taxa
    """
    from .analyze_taxa import (
        YEARS,
        load_datasets,
        plot_newcomers_alluvial,
        plot_top_taxa_trend,
    )

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_datasets(YEARS, data_dir=data_dir)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets to generate plots")

    click.echo(f"Loaded {len(datasets)} datasets (definition: {definition})")

    if not skip_trend:
        click.echo("Generating top taxa trend plot...")
        plot_top_taxa_trend(datasets, output_dir / "top_taxa_trend")

    # Determine which levels to generate
    all_levels = ["Phylum", "Class", "Order", "Family"]
    levels_to_generate = all_levels if level == "all" else [level]

    for lvl in levels_to_generate:
        click.echo(f"Generating alluvial plot ({lvl})...")
        plot_newcomers_alluvial(
            datasets,
            output_dir / f"taxa_newcomers_alluvial_{lvl.lower()}.png",
            taxa_level=lvl,
        )

    click.echo("Taxa analysis complete.")


@analysis.command("length")
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
@click.pass_context
def length(ctx, data_dir, output_dir):
    """Run sequence length distribution analysis.

    Generates a histogram comparing sequence lengths across ToxProt
    time points (2005, 2015, 2025).

    \b
    Examples:
        toxprot analysis length
        toxprot analysis -d all length
        toxprot analysis length -o figures/custom
    """
    from .analyze_sequence_length import plot_sequence_length_histogram

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load datasets
    years = [2005, 2015, 2025]
    datasets = {}
    for year in years:
        filepath = data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets")

    output_path = output_dir / "sequence_length_distribution.png"
    plot_sequence_length_histogram(datasets, output_path)
    click.echo(f"Saved {output_path} (definition: {definition})")


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
    default=10,
    show_default=True,
    help="Number of top protein families to display.",
)
@click.pass_context
def families(ctx, data_dir, output_dir, top_n):
    """Run protein family distribution analysis.

    Generates stacked bar and alluvial charts comparing protein family
    distributions across 2005, 2015, and 2025.

    \b
    Examples:
        toxprot analysis families
        toxprot analysis -d all families
        toxprot analysis families --top-n 10
    """
    from .analyze_protein_families import (
        plot_alluvial_protein_families,
        plot_stacked_bar_protein_families,
    )

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)
    protein_families_dir = output_dir / "protein_families"
    protein_families_dir.mkdir(parents=True, exist_ok=True)

    # Load datasets for comparison (2005, 2015, 2025)
    years = [2005, 2015, 2025]
    datasets = {}
    for year in years:
        filepath = data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets for comparison")

    # Generate stacked bar chart
    families_bar_path = protein_families_dir / "top_families_stacked_bar.png"
    plot_stacked_bar_protein_families(datasets, families_bar_path, top_n=top_n)
    click.echo(f"Saved {families_bar_path} (definition: {definition})")

    # Generate alluvial plot
    alluvial_path = protein_families_dir / "top_families_alluvial.png"
    plot_alluvial_protein_families(datasets, alluvial_path, top_n=top_n)
    click.echo(f"Saved {alluvial_path} (definition: {definition})")


@analysis.command("summary")
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
@click.pass_context
def summary(ctx, data_dir, output_dir):
    """Generate summary statistics table.

    Creates a table comparing key statistics across ToxProt datasets
    (2005, 2015, 2025).

    \b
    Examples:
        toxprot analysis summary
        toxprot analysis -d all summary
        toxprot analysis summary -o figures/custom
    """
    from .analyze_summary import generate_summary_table

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load datasets
    years = [2005, 2015, 2025]
    datasets = {}
    for year in years:
        filepath = data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets")

    output_path = output_dir / "dataset_summary_statistics.png"
    generate_summary_table(datasets, output_path)
    click.echo(f"Saved {output_path} (definition: {definition})")


@analysis.command("go")
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
    default=Path("figures/go_terms"),
    show_default=True,
    help="Directory to save output figures.",
)
@click.option(
    "--top-n",
    type=int,
    default=5,
    show_default=True,
    help="Number of top GO terms to display.",
)
@click.pass_context
def go(ctx, data_dir, output_dir, top_n):
    """Run GO term distribution analysis.

    Generates multiple figures showing GO term growth patterns:
    - Raw counts (explicitly annotated terms)
    - Propagated counts (with GO hierarchy)
    - Annotation depth over time (specificity)
    - New vs growing terms
    - Growth trajectories

    \b
    Examples:
        toxprot analysis go
        toxprot analysis -d all go
        toxprot analysis go --top-n 20
    """
    from .analyze_go_terms import (
        ALL_YEARS,
        generate_all_figures,
        load_datasets,
        load_go_hierarchy,
    )

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    click.echo("Loading GO hierarchy...")
    go_graph = load_go_hierarchy()

    click.echo("Loading ToxProt datasets...")
    datasets = load_datasets(ALL_YEARS, data_dir)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets")

    click.echo(f"Loaded {len(datasets)} years")
    click.echo("Generating figures...")
    generate_all_figures(datasets, output_dir, go_graph, top_n=top_n)

    # Clean up old individual figures
    old_files = [
        "go_terms_comparison.png",
        "go_terms_depth.png",
        "go_terms_depth_distribution.png",
        "go_terms_depth_shift.png",
        "go_terms_new_vs_growth.png",
        "go_terms_propagated.png",
        "go_terms_trends.png",
        "go_terms_coverage_depth.png",
        "go_terms_category_trends.png",
    ]
    for f in old_files:
        old_path = output_dir / f
        if old_path.exists():
            old_path.unlink()

    click.echo(f"Saved figures to {output_dir} (definition: {definition})")


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


@analysis.command("habitat")
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
    default=Path("figures/habitat"),
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
def habitat(ctx, data_dir, output_dir, top_n):
    """Analyze habitat distributions for ToxProt entries.

    Creates a two-panel figure:
    - Panel A: Taxa distribution by habitat (terrestrial vs marine)
    - Panel B: Dual-habitat protein family evolution across 2005, 2015, 2025

    \b
    Examples:
        toxprot analysis habitat
        toxprot analysis -d all habitat
        toxprot analysis habitat --top-n 20
    """
    from .analyze_habitat import YEARS, load_datasets, plot_habitat_combined

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_datasets(YEARS, data_dir)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 3:
        raise click.ClickException("Need at least 3 datasets (2005, 2015, 2025)")

    click.echo(f"Loaded {len(datasets)} datasets (definition: {definition})")

    # Use 2025 data for Panel A
    df_2025 = datasets[max(datasets.keys())]

    output_path = output_dir / "habitat.png"
    plot_habitat_combined(df_2025, datasets, output_path, top_n=top_n)
    click.echo(f"Saved {output_path}")


@analysis.command("source-tissue")
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
    default=Path("figures/source_tissue"),
    show_default=True,
    help="Directory to save output figures.",
)
@click.option(
    "--top-n",
    type=int,
    default=10,
    show_default=True,
    help="Number of top tissues to display individually.",
)
@click.pass_context
def source_tissue(ctx, data_dir, output_dir, top_n):
    """Analyze source tissue evolution over time.

    Creates an alluvial diagram showing tissue flows across 2005, 2015, 2025.
    Tissues are split from semicolon-separated values and counted individually.

    \b
    Examples:
        toxprot analysis source-tissue
        toxprot analysis -d all source-tissue
        toxprot analysis source-tissue --top-n 5
    """
    from .analyze_source_tissue import YEARS, load_datasets, plot_source_tissue_alluvial

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_datasets(YEARS, data_dir=data_dir)

    if definition != "all":
        datasets = {year: filter_by_definition(df, definition) for year, df in datasets.items()}
        datasets = {year: df for year, df in datasets.items() if len(df) > 0}

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets")

    click.echo(f"Loaded {len(datasets)} datasets (definition: {definition})")

    output_path = output_dir / "source_tissue_alluvial.png"
    plot_source_tissue_alluvial(datasets, output_path, top_n=top_n)
    click.echo(f"Saved {output_path}")


@analysis.command("protein-evidence")
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
    default=Path("figures/protein_evidence"),
    show_default=True,
    help="Directory to save output figures.",
)
@click.pass_context
def protein_evidence(ctx, data_dir, output_dir):
    """Analyze protein evidence category evolution.

    Creates an alluvial diagram showing protein evidence level transitions
    across 2008, 2015, and 2025 (PE levels introduced mid-2007).

    Shows flows between categories and removed proteins at each time point.

    \b
    Examples:
        toxprot analysis protein-evidence
        toxprot analysis -d all protein-evidence
        toxprot analysis protein-evidence -o figures/custom
    """
    from .plot_protein_evidence_sankey import (
        YEARS,
        filter_by_definition,
        load_datasets,
        plot_protein_evidence_alluvial,
    )

    definition = ctx.obj["definition"]
    output_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_datasets(YEARS, data_dir=data_dir)

    if len(datasets) < 2:
        raise click.ClickException("Need at least 2 datasets")

    # Apply definition filter for display
    filtered = {y: filter_by_definition(df, definition) for y, df in datasets.items()}

    click.echo(f"Loaded {len(datasets)} datasets (definition: {definition})")
    for year, df in filtered.items():
        click.echo(f"  {year}: {len(df):,} proteins")

    output_path = output_dir / "protein_evidence_sankey.png"
    plot_protein_evidence_alluvial(datasets, output_path, years=YEARS, definition=definition)
    click.echo(f"Saved {output_path}")


@analysis.command()
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data/processed/toxprot"),
    show_default=True,
    help="Directory containing processed CSV files.",
)
@click.pass_context
def pipeline(ctx, data_dir):
    """Run all analysis commands with default parameters.

    Generates all figures: summary, taxa, families, length, habitat,
    source-tissue, ptm, go, protein-evidence, and definitions.

    Each command uses its default output directory under figures/.

    \b
    Examples:
        toxprot analysis pipeline
        toxprot analysis -d all pipeline
        toxprot analysis pipeline --data-dir data/custom
    """
    commands = [
        "summary",
        "taxa",
        "families",
        "length",
        "habitat",
        "source-tissue",
        "ptm",
        "go",
        "protein-evidence",
        "definitions",
    ]

    click.echo(f"Running all analyses (definition: {ctx.obj['definition']})")

    for cmd_name in commands:
        click.echo(f"\n{'='*60}")
        click.echo(f"Running: {cmd_name}")
        click.echo("=" * 60)
        ctx.invoke(analysis.commands[cmd_name], data_dir=data_dir)

    click.echo(f"\nAll analyses complete. Figures saved to figures/")


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
