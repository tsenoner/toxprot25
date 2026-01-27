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


class YearLogger:
    """Logger wrapper that automatically prefixes messages with [year]."""

    def __init__(self, year: int, base_logger: logging.Logger):
        self.year = year
        self.logger = base_logger

    def _fmt(self, msg: str) -> str:
        return f"[{self.year}] {msg}"

    def debug(self, msg: str) -> None:
        self.logger.debug(self._fmt(msg))

    def info(self, msg: str) -> None:
        self.logger.info(self._fmt(msg))

    def warning(self, msg: str) -> None:
        self.logger.warning(self._fmt(msg))

    def error(self, msg: str) -> None:
        self.logger.error(self._fmt(msg))


def setup_logging(verbose: bool = False) -> None:
    """Configure logging for the pipeline.

    Args:
        verbose: If True, set DEBUG level; otherwise WARNING level.
    """
    # Set root logger to WARNING by default, DEBUG if verbose
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG if verbose else logging.WARNING)

    # Configure handler
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s %(levelname)s [%(name)s] - %(message)s" if verbose
        else "%(message)s",
        datefmt="%H:%M:%S"
    ))
    root_logger.addHandler(handler)

    # Set pipeline logger to INFO by default, DEBUG if verbose
    pipeline_logger = logging.getLogger("src.data_processing.pipeline")
    pipeline_logger.setLevel(logging.DEBUG if verbose else logging.INFO)


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
    log = YearLogger(year, logger)

    # File paths
    dat_path = config.raw_dir / f"{year}_sprot.dat"
    tsv_path = config.interim_dir / f"toxprot_{year}.tsv"
    csv_path = config.processed_dir / f"toxprot_{year}.csv"
    fasta_path = config.processed_dir / f"toxprot_{year}.fasta"

    try:
        # Stage 1: Download
        if not dat_path.exists():
            if not download_release(year, config.raw_dir, keep_archive=False):
                return YearResult(year=year, success=False, stage_completed="download", error="Download failed")
        else:
            log.debug("Using existing .dat file")

        # Stage 2: Parse
        if not tsv_path.exists() or not config.skip_existing:
            config.interim_dir.mkdir(parents=True, exist_ok=True)
            if not process_swissprot_file(dat_path, tsv_path, ptm_vocab):
                return YearResult(year=year, success=False, stage_completed="parse", error="Parsing failed")
        else:
            log.debug("Using existing .tsv file")

        # Delete .dat file to save space
        if config.delete_dat_files and dat_path.exists():
            dat_path.unlink()

        # Stage 3: Clean and process
        config.processed_dir.mkdir(parents=True, exist_ok=True)
        process_toxprot_tsv(tsv_path, update_protfams, create_fasta_file)

        # Add taxonomy and habitat classification
        interim_csv_path = tsv_path.with_suffix(".csv")
        interim_fasta_path = tsv_path.with_suffix(".fasta")

        log.debug("Adding taxonomy and habitat data")
        df = process_dataframe_with_taxonomy(interim_csv_path, csv_path, taxdb)

        habitat_mapping_path = config.data_dir / "raw" / "marine_terrestrial.json"
        habitat_detailed_path = config.data_dir / "raw" / "habitat_detailed.json"
        df = add_habitat_classification(csv_path, habitat_mapping_path, habitat_detailed_path)

        # Move FASTA and cleanup
        if interim_fasta_path.exists():
            if fasta_path.exists():
                fasta_path.unlink()
            interim_fasta_path.rename(fasta_path)
        if interim_csv_path.exists():
            interim_csv_path.unlink()
        if config.delete_tsv_files and tsv_path.exists():
            tsv_path.unlink()

        entries_count = len(df) if df is not None else 0
        log.info(f"✓ {entries_count:,} entries")

        return YearResult(year=year, success=True, stage_completed="complete", entries_count=entries_count)

    except Exception as e:
        log.error(f"Error: {e}")
        logger.debug(f"[{year}] Full traceback:", exc_info=True)
        return YearResult(year=year, success=False, stage_completed="unknown", error=str(e))


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

    year_range = f"{min(config.years)}-{max(config.years)}"
    logger.info(f"ToxProt Pipeline: Processing {len(config.years)} years ({year_range})")
    logger.debug("=" * 60)
    logger.debug(f"Directories: raw={config.raw_dir}, interim={config.interim_dir}, processed={config.processed_dir}")
    logger.debug(f"Options: delete_dat={config.delete_dat_files}, delete_tsv={config.delete_tsv_files}, skip_existing={config.skip_existing}")
    logger.debug("=" * 60)

    # Initialize shared resources
    logger.debug("Initializing resources (PTM vocabulary, taxonomy database)...")
    vocab_manager = PTMVocabulary(str(config.data_dir / "raw"))
    vocab_manager.ensure_ptmlist_available()
    vocab_manager.load_ptmlist()
    ptm_vocab = vocab_manager.ptm_vocab
    taxdb = initialize_taxdb()

    # Process each year
    results = []
    for year in sorted(config.years):
        if config.skip_existing and is_year_complete(year, config):
            # Read cached entry count
            csv_path = config.processed_dir / f"toxprot_{year}.csv"
            try:
                import pandas as pd
                entries_count = len(pd.read_csv(csv_path))
                logger.info(f"[{year}] ✓ {entries_count:,} entries (cached)")
            except Exception:
                entries_count = None
                logger.info(f"[{year}] ✓ cached")

            results.append(YearResult(year=year, success=True, stage_completed="complete", entries_count=entries_count))
            continue

        results.append(process_single_year(year, config, ptm_vocab, taxdb))

    # Summary
    successful = sum(1 for r in results if r.success)
    total_entries = sum(r.entries_count or 0 for r in results if r.success)
    logger.info("")
    logger.info(f"Complete: {successful}/{len(results)} years, {total_entries:,} total entries")

    # Report failures
    failures = [r for r in results if not r.success]
    if failures:
        logger.warning(f"Failed: {len(failures)} years")
        for r in failures:
            logger.warning(f"  [{r.year}] {r.error}")

    logger.debug("=" * 60)
    return results


if __name__ == "__main__":
    # When run as module, use the click CLI
    from .cli import pipeline as pipeline_cmd
    pipeline_cmd()
