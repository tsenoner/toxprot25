#!/usr/bin/env python3
"""
Unified ToxProt Data Processing Pipeline.

Orchestrates the download, parsing, and cleaning of UniProt releases year-by-year
to minimize disk usage (~2.5GB instead of ~50GB). Each year processes through:
download .dat → parse to .tsv → clean to .csv/.fasta → delete intermediates.

Usage:
    # Process all years (2005-2025)
    toxprot data pipeline

    # Process specific years
    toxprot data pipeline -y 2020-2025

    # Force reprocess a single year
    toxprot data pipeline -y 2025 --force

    # Keep intermediate files for debugging
    toxprot data pipeline -y 2025 --keep-dat --keep-tsv
"""

import logging
import sys

from .clean_data import (
    add_habitat_classification,
    create_fasta_file,
    initialize_taxdb,
    process_dataframe_with_taxonomy,
    process_toxprot_tsv,
    update_protfams,
)
from .config import PipelineConfig, YearResult
from .download_uniprot_releases import download_release
from .parse_sprot_dat import (
    PTMVocabulary,
    process_swissprot_file,
)

# Module logger
logger = logging.getLogger(__name__)


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for the pipeline.

    Args:
        verbose: If True, set DEBUG level; otherwise INFO level.
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s [%(name)s] - %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
        ],
    )


def is_year_complete(year: int, config: PipelineConfig) -> bool:
    """Check if a year has been fully processed.

    A year is considered complete if both the final CSV and FASTA files exist.

    Args:
        year: The year to check.
        config: Pipeline configuration.

    Returns:
        True if both output files exist, False otherwise.
    """
    csv_path = config.processed_dir / f"toxprot_{year}.csv"
    fasta_path = config.processed_dir / f"toxprot_{year}.fasta"
    return csv_path.exists() and fasta_path.exists()


def process_single_year(
    year: int,
    config: PipelineConfig,
    ptm_vocab: dict | None = None,
    taxdb=None,
) -> YearResult:
    """Process a single year through all pipeline stages.

    Stages:
        1. Download .dat file (skip if exists)
        2. Parse to .tsv (filter Metazoa + venom tissue)
        3. Delete .dat file (if configured)
        4. Clean → final .csv + .fasta
        5. Delete intermediate .tsv (if configured)

    Args:
        year: The year to process.
        config: Pipeline configuration.
        ptm_vocab: Pre-loaded PTM vocabulary (optional).
        taxdb: Pre-initialized taxonomy database (optional).

    Returns:
        YearResult with processing outcome.
    """
    logger.info(f"[{year}] Starting processing...")

    # File paths
    dat_path = config.raw_dir / f"{year}_sprot.dat"
    tsv_path = config.interim_dir / f"toxprot_{year}.tsv"
    csv_path = config.processed_dir / f"toxprot_{year}.csv"
    fasta_path = config.processed_dir / f"toxprot_{year}.fasta"

    try:
        # Stage 1: Download
        logger.info(f"[{year}] Stage 1: Download")
        if not dat_path.exists():
            if not download_release(year, config.raw_dir, keep_archive=False):
                return YearResult(
                    year=year,
                    success=False,
                    stage_completed="download",
                    error="Download failed",
                )
        else:
            logger.info(f"[{year}]   .dat file exists, skipping download")

        # Stage 2: Parse
        logger.info(f"[{year}] Stage 2: Parse")
        if not tsv_path.exists():
            config.interim_dir.mkdir(parents=True, exist_ok=True)
            success = process_swissprot_file(dat_path, tsv_path, ptm_vocab)
            if not success:
                return YearResult(
                    year=year,
                    success=False,
                    stage_completed="parse",
                    error="Parsing failed",
                )
        else:
            logger.info(f"[{year}]   .tsv file exists, skipping parse")

        # Stage 3: Delete .dat file (if configured)
        if config.delete_dat_files and dat_path.exists():
            logger.info(f"[{year}] Deleting .dat file to save disk space")
            dat_path.unlink()

        # Stage 4: Clean and process
        logger.info(f"[{year}] Stage 3: Clean")
        config.processed_dir.mkdir(parents=True, exist_ok=True)

        # Process TSV to CSV and FASTA (in interim directory as intermediate)
        process_toxprot_tsv(tsv_path, update_protfams, create_fasta_file)

        # Intermediate files created next to input
        interim_csv_path = tsv_path.with_suffix(".csv")
        interim_fasta_path = tsv_path.with_suffix(".fasta")

        # Add taxonomic information
        logger.info(f"[{year}]   Adding taxonomy...")
        df = process_dataframe_with_taxonomy(interim_csv_path, csv_path, taxdb)

        # Add habitat classification
        logger.info(f"[{year}]   Adding habitat classification...")
        habitat_mapping_path = config.data_dir / "raw" / "marine_terrestrial.json"
        habitat_detailed_path = config.data_dir / "raw" / "habitat_detailed.json"
        df = add_habitat_classification(csv_path, habitat_mapping_path, habitat_detailed_path)

        # Move FASTA to output directory
        if interim_fasta_path.exists():
            if fasta_path.exists():
                fasta_path.unlink()
            interim_fasta_path.rename(fasta_path)

        # Clean up intermediate CSV
        if interim_csv_path.exists():
            interim_csv_path.unlink()

        entries_count = len(df) if df is not None else 0

        # Stage 5: Delete intermediate TSV (if configured)
        if config.delete_tsv_files and tsv_path.exists():
            logger.info(f"[{year}] Deleting intermediate .tsv file")
            tsv_path.unlink()

        logger.info(f"[{year}] Complete! {entries_count} entries")
        return YearResult(
            year=year,
            success=True,
            stage_completed="complete",
            entries_count=entries_count,
        )

    except Exception as e:
        logger.error(f"[{year}] Error: {e}", exc_info=True)
        return YearResult(
            year=year,
            success=False,
            stage_completed="unknown",
            error=str(e),
        )


def run_pipeline(config: PipelineConfig) -> list[YearResult]:
    """Run the full pipeline for all configured years.

    Args:
        config: Pipeline configuration.

    Returns:
        List of YearResult objects, one per year.
    """
    # Validate configuration
    errors = config.validate()
    if errors:
        logger.error("Configuration errors:")
        for error in errors:
            logger.error(f"  - {error}")
        return []

    # Create directories
    config.raw_dir.mkdir(parents=True, exist_ok=True)
    config.interim_dir.mkdir(parents=True, exist_ok=True)
    config.processed_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("ToxProt Unified Data Processing Pipeline")
    logger.info("=" * 60)
    logger.info(f"Years to process: {min(config.years)}-{max(config.years)} ({len(config.years)} years)")
    logger.info(f"Raw directory: {config.raw_dir}")
    logger.info(f"Interim directory: {config.interim_dir}")
    logger.info(f"Processed directory: {config.processed_dir}")
    logger.info(f"Delete .dat files: {config.delete_dat_files}")
    logger.info(f"Delete .tsv files: {config.delete_tsv_files}")
    logger.info(f"Skip existing: {config.skip_existing}")
    logger.info("=" * 60)

    # Initialize shared resources
    logger.info("Initializing PTM vocabulary...")
    vocab_manager = PTMVocabulary(str(config.data_dir / "raw"))
    vocab_manager.ensure_ptmlist_available()
    vocab_manager.load_ptmlist()
    ptm_vocab = vocab_manager.ptm_vocab

    logger.info("Initializing taxonomy database...")
    taxdb = initialize_taxdb()

    # Process each year
    results = []
    for year in sorted(config.years):
        # Check if already complete
        if config.skip_existing and is_year_complete(year, config):
            logger.info(f"[{year}] Already complete, skipping")
            # Read existing file to get entry count
            csv_path = config.processed_dir / f"toxprot_{year}.csv"
            try:
                import pandas as pd
                df = pd.read_csv(csv_path)
                entries_count = len(df)
            except Exception:
                entries_count = None
            results.append(YearResult(
                year=year,
                success=True,
                stage_completed="complete",
                entries_count=entries_count,
            ))
            continue

        result = process_single_year(year, config, ptm_vocab, taxdb)
        results.append(result)

    # Summary
    logger.info("=" * 60)
    logger.info("Pipeline Summary")
    logger.info("=" * 60)
    successful = sum(1 for r in results if r.success)
    total_entries = sum(r.entries_count or 0 for r in results if r.success)
    logger.info(f"Successful: {successful}/{len(results)} years")
    logger.info(f"Total entries: {total_entries}")

    # Report failures
    failures = [r for r in results if not r.success]
    if failures:
        logger.warning("Failed years:")
        for r in failures:
            logger.warning(f"  {r}")

    logger.info("=" * 60)
    return results


if __name__ == "__main__":
    # When run as module, use the click CLI
    from .cli import pipeline as pipeline_cmd
    pipeline_cmd()
