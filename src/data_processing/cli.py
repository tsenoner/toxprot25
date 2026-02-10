"""CLI commands for data processing."""

import logging
import sys
from pathlib import Path

import click

from ..config import ALL_YEARS, DATA_DIR, INTERIM_DIR, MAX_YEAR, MIN_YEAR, RAW_DIR
from .config import PipelineConfig
from .pipeline import run_pipeline, setup_logging


def validate_years(ctx, param, value):
    """Validate and parse years argument."""
    if not value:
        return list(ALL_YEARS)

    years = []
    for item in value:
        if "-" in item:
            # Range like "2020-2025"
            try:
                start, end = map(int, item.split("-"))
                years.extend(range(start, end + 1))
            except ValueError:
                raise click.BadParameter(f"Invalid year range: {item}") from None
        else:
            try:
                years.append(int(item))
            except ValueError:
                raise click.BadParameter(f"Invalid year: {item}") from None

    # Validate all years
    for year in years:
        if year < MIN_YEAR or year > MAX_YEAR:
            raise click.BadParameter(f"Year {year} out of range ({MIN_YEAR}-{MAX_YEAR})")

    return sorted(set(years))


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
def data():
    """Data processing commands."""
    pass


@data.command()
@click.option(
    "--years",
    "-y",
    multiple=True,
    callback=validate_years,
    help="Years to process. Can be single years (2020) or ranges (2020-2025). "
    "Can be specified multiple times. Default: all years (2005-2025).",
)
@click.option(
    "--raw-dir",
    type=click.Path(path_type=Path),
    default=RAW_DIR,
    show_default=True,
    help="Directory for raw .dat files.",
)
@click.option(
    "--interim-dir",
    type=click.Path(path_type=Path),
    default=INTERIM_DIR,
    show_default=True,
    help="Directory for intermediate .tsv files.",
)
@click.option(
    "--processed-dir",
    type=click.Path(path_type=Path),
    default=DATA_DIR,
    show_default=True,
    help="Directory for final output files.",
)
@click.option(
    "--data-dir",
    type=click.Path(path_type=Path),
    default=Path("data"),
    show_default=True,
    help="Base data directory for auxiliary files.",
)
@click.option(
    "--keep-raw/--delete-raw",
    default=False,
    show_default=True,
    help="Keep raw XML files after parsing (default: delete to save space).",
)
@click.option(
    "--keep-tsv/--delete-tsv",
    default=True,
    show_default=True,
    help="Keep intermediate .tsv files after cleaning.",
)
@click.option(
    "--force/--no-force",
    "-f",
    default=False,
    show_default=True,
    help="Force reprocessing even if output exists.",
)
@click.option(
    "--verbose/--quiet",
    "-v/-q",
    default=False,
    help="Enable verbose logging.",
)
def pipeline(
    years, raw_dir, interim_dir, processed_dir, data_dir, keep_raw, keep_tsv, force, verbose
):
    """Run the unified data processing pipeline.

    Downloads UniProt releases, parses them to extract toxin proteins,
    and cleans the data for analysis. Processes year-by-year to minimize
    disk usage (~2.5GB instead of ~50GB).

    \b
    Examples:
        # Process all years (2005-2025)
        toxprot data pipeline

        # Process specific years
        toxprot data pipeline -y 2020 -y 2021 -y 2022

        # Process a range of years
        toxprot data pipeline -y 2020-2025

        # Force reprocess a single year
        toxprot data pipeline -y 2025 --force

        # Keep intermediate files for debugging
        toxprot data pipeline -y 2025 --keep-raw --keep-tsv
    """
    setup_logging(verbose=verbose)

    config = PipelineConfig(
        years=years,
        raw_dir=raw_dir,
        interim_dir=interim_dir,
        processed_dir=processed_dir,
        data_dir=data_dir,
        delete_raw_files=not keep_raw,
        delete_tsv_files=not keep_tsv,
        skip_existing=not force,
    )

    results = run_pipeline(config)

    # Exit with error code if any failures
    if any(not r.success for r in results):
        sys.exit(1)


@data.command()
@click.option(
    "--years",
    "-y",
    multiple=True,
    callback=validate_years,
    help="Years to download. Default: all years (2005-2025).",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=RAW_DIR,
    show_default=True,
    help="Directory to save downloaded files.",
)
@click.option(
    "--keep-archives/--no-keep-archives",
    default=False,
    show_default=True,
    help="Keep the .tar.gz archives after extraction.",
)
@click.option(
    "--list-only",
    is_flag=True,
    help="List available releases without downloading.",
)
def download(years, output_dir, keep_archives, list_only):
    """Download UniProt Swiss-Prot releases.

    Downloads the first release of each year from the UniProt FTP server.
    Extracts .dat files from the archives.

    \b
    Examples:
        # Download all releases
        toxprot data download

        # Download specific years
        toxprot data download -y 2020-2025

        # List available releases
        toxprot data download --list-only
    """
    from .download_uniprot_releases import RELEASES, download_release

    if list_only:
        click.echo("Available releases:")
        for year, (release_dir, archive) in sorted(RELEASES.items()):
            click.echo(f"  {year}: {release_dir}/knowledgebase/{archive}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)

    click.echo(f"Output directory: {output_dir}")
    click.echo(f"Years to download: {min(years)}-{max(years)}")
    click.echo("=" * 60)

    success_count = 0
    for year in sorted(years):
        click.echo(f"\n[{year}]")
        if download_release(year, output_dir, keep_archives):
            success_count += 1

    click.echo("\n" + "=" * 60)
    click.echo(f"Downloaded {success_count}/{len(years)} releases")


@data.command()
@click.argument(
    "input_files",
    nargs=-1,
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True, path_type=Path),
    help="Directory containing SwissProt .xml files.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=INTERIM_DIR,
    show_default=True,
    help="Directory for output files.",
)
@click.option(
    "--years",
    "-y",
    multiple=True,
    callback=validate_years,
    help="Years to process (looks for {year}_sprot.xml in input-dir).",
)
@click.option(
    "--delete-input/--keep-input",
    default=False,
    show_default=True,
    help="Delete input .xml files after successful processing.",
)
def parse(input_files, input_dir, output_dir, years, delete_input):
    """Parse Swiss-Prot .xml files to extract toxin proteins.

    Filters for Metazoa organisms with venom tissue specificity.
    Outputs TSV files with protein metadata and sequences.

    \b
    Examples:
        # Parse a single file
        toxprot data parse data/raw/uniprot_releases/2025_sprot.xml

        # Parse all files in a directory
        toxprot data parse -i data/raw/uniprot_releases

        # Parse specific years
        toxprot data parse -i data/raw/uniprot_releases -y 2020-2025
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
    )
    logger = logging.getLogger(__name__)

    from .parse_sprot_dat import get_output_filename
    from .parse_sprot_xml import process_xml_file

    # Determine input files
    files_to_process = []

    if years and years != ALL_YEARS:
        # Process specific years from input-dir
        search_dir = input_dir or RAW_DIR
        for year in years:
            input_path = search_dir / f"{year}_sprot.xml"
            if input_path.exists():
                files_to_process.append(input_path)
            else:
                logger.warning(f"File not found: {input_path}")
    elif input_dir:
        # Process all .xml files in directory
        files_to_process = sorted(input_dir.glob("*_sprot.xml"))
    elif input_files:
        # Process specified files
        files_to_process = list(input_files)
    else:
        raise click.UsageError("Specify input files, --input-dir, or --years")

    if not files_to_process:
        raise click.ClickException("No input files found")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("SwissProt XML Parser for Tox-Prot")
    logger.info("=" * 60)
    logger.info(f"Input files: {len(files_to_process)}")
    logger.info(f"Output directory: {output_dir}")

    # Process each file
    successful = 0
    for input_path in files_to_process:
        output_filename = get_output_filename(input_path)
        output_path = output_dir / output_filename

        logger.info(f"Processing: {input_path.name}")

        if process_xml_file(input_path, output_path):
            successful += 1
            if delete_input:
                logger.info(f"Deleting input file: {input_path}")
                input_path.unlink()

    logger.info("=" * 60)
    logger.info(f"Processed {successful}/{len(files_to_process)} files successfully")


@data.command()
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True, path_type=Path),
    default=INTERIM_DIR,
    show_default=True,
    help="Directory containing toxprot_*.tsv files.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    default=DATA_DIR,
    show_default=True,
    help="Directory for output files.",
)
@click.option(
    "--years",
    "-y",
    multiple=True,
    callback=validate_years,
    help="Specific years to process. Default: all found in input-dir.",
)
@click.option(
    "--data-dir",
    type=click.Path(exists=True, path_type=Path),
    default=Path("data"),
    show_default=True,
    help="Base data directory for habitat mappings.",
)
def clean(input_dir, output_dir, years, data_dir):
    """Clean and process parsed ToxProt data.

    Adds taxonomy information, standardizes protein family names,
    classifies habitats, and creates FASTA files with signal peptide removal.

    \b
    Examples:
        # Clean all parsed files
        toxprot data clean

        # Clean specific years
        toxprot data clean -y 2020-2025
    """
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
    )
    logger = logging.getLogger(__name__)

    from .clean_data import (
        add_habitat_classification,
        create_fasta_file,
        get_year_from_filename,
        initialize_taxdb,
        process_dataframe_with_taxonomy,
        process_toxprot_tsv,
        update_protfams,
    )

    # Define paths for habitat mappings
    habitat_mapping_path = data_dir / "raw" / "marine_terrestrial.json"
    habitat_detailed_path = data_dir / "raw" / "habitat_detailed.json"

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find input files
    if years and years != ALL_YEARS:
        input_files = [input_dir / f"toxprot_{year}.tsv" for year in years]
        input_files = [f for f in input_files if f.exists()]
    else:
        input_files = sorted(input_dir.glob("toxprot_*.tsv"))

    if not input_files:
        raise click.ClickException(f"No toxprot_*.tsv files found in {input_dir}")

    logger.info("=" * 60)
    logger.info("ToxProt Data Cleaning Pipeline")
    logger.info("=" * 60)
    logger.info(f"Input directory: {input_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Files to process: {len(input_files)}")

    # Initialize taxonomy database
    logger.info("Initializing taxonomy database...")
    taxdb = initialize_taxdb()

    successful = 0
    for tsv_path in input_files:
        year = get_year_from_filename(tsv_path)
        logger.info(f"[{year}] Processing {tsv_path.name}...")

        # Step 1: Process TSV to CSV and FASTA
        process_toxprot_tsv(tsv_path, update_protfams, create_fasta_file)

        # Intermediate files
        interim_csv_path = tsv_path.with_suffix(".csv")
        interim_fasta_path = tsv_path.with_suffix(".fasta")

        # Final output paths
        processed_csv_path = output_dir / f"toxprot_{year}.csv"
        processed_fasta_path = output_dir / f"toxprot_{year}.fasta"

        # Step 2: Add taxonomic information
        df = process_dataframe_with_taxonomy(interim_csv_path, processed_csv_path, taxdb)

        # Step 3: Add habitat classification
        df = add_habitat_classification(
            processed_csv_path, habitat_mapping_path, habitat_detailed_path
        )

        # Move FASTA to output directory
        interim_fasta_path.rename(processed_fasta_path)

        # Clean up intermediate CSV
        interim_csv_path.unlink()

        logger.info(f"    Entries: {len(df)}")
        successful += 1

    logger.info("=" * 60)
    logger.info(f"Processed {successful}/{len(input_files)} files successfully")
