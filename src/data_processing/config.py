"""Configuration dataclasses for the ToxProt data processing pipeline."""

from dataclasses import dataclass, field
from pathlib import Path

from ..config import ALL_YEARS, DATA_DIR, INTERIM_DIR, MAX_YEAR, MIN_YEAR, RAW_DIR


@dataclass
class PipelineConfig:
    """Configuration for the ToxProt data processing pipeline.

    Attributes:
        years: List of years to process (default: 2005-2025)
        raw_dir: Directory for raw UniProt .dat files
        interim_dir: Directory for intermediate parsed TSV files
        processed_dir: Directory for final processed CSV/FASTA files
        delete_raw_files: Whether to delete raw files (.xml/.dat) after parsing (saves disk space)
        delete_tsv_files: Whether to delete intermediate TSV files after cleaning
        skip_existing: Whether to skip years that are already processed (resume support)
        data_dir: Base data directory for auxiliary files (PTM vocabulary, habitat mappings)
    """

    years: list[int] = field(default_factory=lambda: list(ALL_YEARS))
    raw_dir: Path = field(default_factory=lambda: Path(RAW_DIR))
    interim_dir: Path = field(default_factory=lambda: Path(INTERIM_DIR))
    processed_dir: Path = field(default_factory=lambda: Path(DATA_DIR))
    delete_raw_files: bool = True
    delete_tsv_files: bool = False
    skip_existing: bool = True
    data_dir: Path = field(default_factory=lambda: Path("data"))

    def __post_init__(self):
        """Convert string paths to Path objects if needed."""
        if isinstance(self.raw_dir, str):
            self.raw_dir = Path(self.raw_dir)
        if isinstance(self.interim_dir, str):
            self.interim_dir = Path(self.interim_dir)
        if isinstance(self.processed_dir, str):
            self.processed_dir = Path(self.processed_dir)
        if isinstance(self.data_dir, str):
            self.data_dir = Path(self.data_dir)

    def validate(self) -> list[str]:
        """Validate configuration and return list of errors."""
        errors = []

        # Validate years
        if not self.years:
            errors.append("No years specified")
        else:
            for year in self.years:
                if not isinstance(year, int) or year < MIN_YEAR or year > MAX_YEAR:
                    errors.append(f"Invalid year: {year} (must be {MIN_YEAR}-{MAX_YEAR})")

        # Check data directory exists
        if not self.data_dir.exists():
            errors.append(f"Data directory does not exist: {self.data_dir}")

        # Check habitat mapping files exist
        habitat_mapping = self.data_dir / "raw" / "marine_terrestrial.json"
        habitat_detailed = self.data_dir / "raw" / "habitat_detailed.json"
        if not habitat_mapping.exists():
            errors.append(f"Habitat mapping not found: {habitat_mapping}")
        if not habitat_detailed.exists():
            errors.append(f"Detailed habitat mapping not found: {habitat_detailed}")

        return errors


@dataclass
class YearResult:
    """Result of processing a single year through the pipeline.

    Attributes:
        year: The year that was processed
        success: Whether all stages completed successfully
        stage_completed: The last stage that completed ('download', 'parse', 'clean', 'complete')
        entries_count: Number of entries in the final output (if successful)
        error: Error message if processing failed
    """

    year: int
    success: bool
    stage_completed: str  # 'download', 'parse', 'clean', 'complete'
    entries_count: int | None = None
    error: str | None = None

    def __str__(self) -> str:
        if self.success:
            return f"Year {self.year}: SUCCESS ({self.entries_count} entries)"
        else:
            return f"Year {self.year}: FAILED at {self.stage_completed} - {self.error}"
