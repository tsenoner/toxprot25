"""Data processing module for ToxProt25.

This module provides tools for processing UniProt Swiss-Prot releases to extract
toxin-related protein data. The main entry point is the unified pipeline.

Main components:
    - config: Configuration dataclasses (PipelineConfig, YearResult)
    - pipeline: Unified processing pipeline (run_pipeline, process_single_year)
    - download_uniprot_releases: Download UniProt releases
    - parse_sprot_dat: Parse Swiss-Prot .dat files
    - clean_data: Clean and process parsed data
"""

from .config import PipelineConfig, YearResult
from .pipeline import is_year_complete, process_single_year, run_pipeline

__all__ = [
    "PipelineConfig",
    "YearResult",
    "run_pipeline",
    "process_single_year",
    "is_year_complete",
]
