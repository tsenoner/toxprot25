"""Central project constants for ToxProt analysis pipeline.

Single source of truth for years, directory paths, and other shared
constants used across analysis and data processing modules.
"""

from pathlib import Path

# Year constants
MIN_YEAR = 2005
MAX_YEAR = 2025
ALL_YEARS = list(range(MIN_YEAR, MAX_YEAR + 1))
COMPARISON_YEARS = [2005, 2015, 2025]
PE_COMPARISON_YEARS = [2008, 2015, 2025]

# Directory paths
DATA_DIR = Path("data/processed/toxprot")
FIGURES_DIR = Path("figures")
RAW_DIR = Path("data/raw/uniprot_releases")
INTERIM_DIR = Path("data/interim/toxprot_parsed")
