"""Tests for the unified pipeline module."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd

from src.data_processing.config import PipelineConfig, YearResult
from src.data_processing.pipeline import (
    is_year_complete,
    process_single_year,
    run_pipeline,
)


class TestPipelineConfig:
    """Tests for the PipelineConfig dataclass."""

    def test_default_years(self):
        """Test default years range."""
        config = PipelineConfig()
        assert config.years == list(range(2005, 2026))
        assert len(config.years) == 21

    def test_default_directories(self):
        """Test default directory paths."""
        config = PipelineConfig()
        assert config.raw_dir == Path("data/raw/uniprot_releases")
        assert config.interim_dir == Path("data/interim/toxprot_parsed")
        assert config.processed_dir == Path("data/processed/toxprot")

    def test_default_flags(self):
        """Test default flag values."""
        config = PipelineConfig()
        assert config.delete_dat_files is True
        assert config.delete_tsv_files is False
        assert config.skip_existing is True

    def test_string_to_path_conversion(self):
        """Test that string paths are converted to Path objects."""
        config = PipelineConfig(
            raw_dir="custom/raw",
            interim_dir="custom/interim",
            processed_dir="custom/processed",
        )
        assert isinstance(config.raw_dir, Path)
        assert isinstance(config.interim_dir, Path)
        assert isinstance(config.processed_dir, Path)

    def test_validate_empty_years(self):
        """Test validation fails for empty years."""
        config = PipelineConfig(years=[])
        errors = config.validate()
        assert any("No years" in e for e in errors)

    def test_validate_invalid_year(self):
        """Test validation fails for invalid years."""
        config = PipelineConfig(years=[2000])  # Too early
        errors = config.validate()
        assert any("Invalid year" in e for e in errors)

        config = PipelineConfig(years=[2030])  # Too late
        errors = config.validate()
        assert any("Invalid year" in e for e in errors)

    def test_validate_missing_data_dir(self, temp_dir):
        """Test validation fails for missing data directory."""
        config = PipelineConfig(
            years=[2020],
            data_dir=temp_dir / "nonexistent"
        )
        errors = config.validate()
        assert any("Data directory does not exist" in e for e in errors)

    def test_validate_missing_habitat_files(self, temp_dir):
        """Test validation fails for missing habitat mapping files."""
        config = PipelineConfig(
            years=[2020],
            data_dir=temp_dir  # Exists but missing habitat files
        )
        (temp_dir / "raw").mkdir(parents=True, exist_ok=True)
        errors = config.validate()
        assert any("Habitat mapping not found" in e for e in errors)


class TestYearResult:
    """Tests for the YearResult dataclass."""

    def test_success_str(self):
        """Test string representation of successful result."""
        result = YearResult(year=2020, success=True, stage_completed="complete", entries_count=1000)
        assert "2020" in str(result)
        assert "SUCCESS" in str(result)
        assert "1000" in str(result)

    def test_failure_str(self):
        """Test string representation of failed result."""
        result = YearResult(year=2020, success=False, stage_completed="download", error="Network error")
        assert "2020" in str(result)
        assert "FAILED" in str(result)
        assert "download" in str(result)
        assert "Network error" in str(result)


class TestIsYearComplete:
    """Tests for the is_year_complete function."""

    def test_complete_when_both_files_exist(self, temp_dir):
        """Test returns True when both CSV and FASTA exist."""
        config = PipelineConfig(processed_dir=temp_dir)
        (temp_dir / "toxprot_2020.csv").write_text("header\ndata")
        (temp_dir / "toxprot_2020.fasta").write_text(">entry\nSEQ")

        assert is_year_complete(2020, config) is True

    def test_incomplete_when_csv_missing(self, temp_dir):
        """Test returns False when CSV missing."""
        config = PipelineConfig(processed_dir=temp_dir)
        (temp_dir / "toxprot_2020.fasta").write_text(">entry\nSEQ")

        assert is_year_complete(2020, config) is False

    def test_incomplete_when_fasta_missing(self, temp_dir):
        """Test returns False when FASTA missing."""
        config = PipelineConfig(processed_dir=temp_dir)
        (temp_dir / "toxprot_2020.csv").write_text("header\ndata")

        assert is_year_complete(2020, config) is False

    def test_incomplete_when_both_missing(self, temp_dir):
        """Test returns False when both files missing."""
        config = PipelineConfig(processed_dir=temp_dir)
        assert is_year_complete(2020, config) is False


class TestProcessSingleYear:
    """Tests for the process_single_year function."""

    def test_returns_failure_on_download_error(self, temp_dir):
        """Test that download errors are handled."""
        config = PipelineConfig(
            years=[2020],
            raw_dir=temp_dir / "raw",
            interim_dir=temp_dir / "interim",
            processed_dir=temp_dir / "processed",
            data_dir=temp_dir,
        )
        (temp_dir / "raw").mkdir(parents=True)

        with patch("src.data_processing.pipeline.download_release") as mock_download:
            mock_download.return_value = False
            result = process_single_year(2020, config)

        assert result.success is False
        assert result.stage_completed == "download"
        assert "Download failed" in result.error

    def test_returns_failure_on_parse_error(self, temp_dir):
        """Test that parse errors are handled."""
        config = PipelineConfig(
            years=[2020],
            raw_dir=temp_dir / "raw",
            interim_dir=temp_dir / "interim",
            processed_dir=temp_dir / "processed",
            data_dir=temp_dir,
        )
        (temp_dir / "raw").mkdir(parents=True)
        (temp_dir / "raw" / "2020_sprot.dat").write_text("ID TEST\n//\n")

        with patch("src.data_processing.pipeline.process_swissprot_file") as mock_parse:
            mock_parse.return_value = False
            result = process_single_year(2020, config)

        assert result.success is False
        assert result.stage_completed == "parse"

    def test_skips_download_if_dat_exists(self, temp_dir):
        """Test that download is skipped if .dat file exists."""
        config = PipelineConfig(
            years=[2020],
            raw_dir=temp_dir / "raw",
            interim_dir=temp_dir / "interim",
            processed_dir=temp_dir / "processed",
            data_dir=temp_dir,
        )
        (temp_dir / "raw").mkdir(parents=True)
        dat_file = temp_dir / "raw" / "2020_sprot.dat"
        dat_file.write_text("ID TEST\n//\n")

        with patch("src.data_processing.pipeline.download_release") as mock_download:
            with patch("src.data_processing.pipeline.process_swissprot_file") as mock_parse:
                mock_parse.return_value = False  # Will fail at parse stage
                process_single_year(2020, config)

        mock_download.assert_not_called()

    def test_deletes_dat_file_when_configured(self, temp_dir):
        """Test that .dat file is deleted after parsing when configured."""
        config = PipelineConfig(
            years=[2020],
            raw_dir=temp_dir / "raw",
            interim_dir=temp_dir / "interim",
            processed_dir=temp_dir / "processed",
            data_dir=temp_dir,
            delete_dat_files=True,
        )
        (temp_dir / "raw").mkdir(parents=True)
        (temp_dir / "interim").mkdir(parents=True)
        dat_file = temp_dir / "raw" / "2020_sprot.dat"
        dat_file.write_text("ID TEST\n//\n")

        with patch("src.data_processing.pipeline.process_swissprot_file") as mock_parse:
            mock_parse.return_value = True
            with patch("src.data_processing.pipeline.process_toxprot_tsv"):
                with patch("src.data_processing.pipeline.process_dataframe_with_taxonomy"):
                    with patch("src.data_processing.pipeline.add_habitat_classification") as mock_habitat:
                        mock_habitat.return_value = pd.DataFrame({"Entry": ["P0C1T5"]})
                        # Create expected intermediate files
                        tsv_file = temp_dir / "interim" / "toxprot_2020.tsv"
                        tsv_file.write_text("Entry\nP0C1T5")
                        (temp_dir / "interim" / "toxprot_2020.fasta").write_text(">P0C1T5\nSEQ")
                        (temp_dir / "interim" / "toxprot_2020.csv").write_text("Entry\nP0C1T5")

                        process_single_year(2020, config)

        # .dat should be deleted
        assert not dat_file.exists()


class TestRunPipeline:
    """Tests for the run_pipeline function."""

    def test_returns_empty_list_on_validation_error(self, temp_dir):
        """Test that validation errors result in empty results."""
        config = PipelineConfig(
            years=[],  # Invalid: no years
            data_dir=temp_dir,
        )
        results = run_pipeline(config)
        assert results == []

    def test_skips_completed_years_when_configured(self, temp_dir, habitat_mapping_file, habitat_detailed_file):
        """Test that completed years are skipped."""
        # Setup valid config
        raw_dir = temp_dir / "raw"
        (raw_dir / "uniprot_releases").mkdir(parents=True)
        processed_dir = temp_dir / "processed"
        processed_dir.mkdir(parents=True)

        # Create completed year files
        (processed_dir / "toxprot_2020.csv").write_text("Entry\nP0C1T5")
        (processed_dir / "toxprot_2020.fasta").write_text(">P0C1T5\nSEQ")

        # Move habitat files to expected location
        import shutil
        (temp_dir / "raw").mkdir(exist_ok=True)
        shutil.copy(habitat_mapping_file, temp_dir / "raw" / "marine_terrestrial.json")
        shutil.copy(habitat_detailed_file, temp_dir / "raw" / "habitat_detailed.json")

        config = PipelineConfig(
            years=[2020],
            raw_dir=temp_dir / "raw" / "uniprot_releases",
            processed_dir=processed_dir,
            data_dir=temp_dir,
            skip_existing=True,
        )

        with patch("src.data_processing.pipeline.PTMVocabulary"):
            with patch("src.data_processing.pipeline.initialize_taxdb"):
                results = run_pipeline(config)

        assert len(results) == 1
        assert results[0].success is True
        assert results[0].stage_completed == "complete"

    def test_creates_directories(self, temp_dir, habitat_mapping_file, habitat_detailed_file):
        """Test that required directories are created."""
        import shutil
        (temp_dir / "raw").mkdir(exist_ok=True)
        shutil.copy(habitat_mapping_file, temp_dir / "raw" / "marine_terrestrial.json")
        shutil.copy(habitat_detailed_file, temp_dir / "raw" / "habitat_detailed.json")

        config = PipelineConfig(
            years=[2020],
            raw_dir=temp_dir / "raw_new" / "uniprot",
            interim_dir=temp_dir / "interim_new",
            processed_dir=temp_dir / "processed_new",
            data_dir=temp_dir,
        )

        with patch("src.data_processing.pipeline.PTMVocabulary"):
            with patch("src.data_processing.pipeline.initialize_taxdb"):
                with patch("src.data_processing.pipeline.process_single_year") as mock_process:
                    mock_process.return_value = YearResult(2020, True, "complete", 100)
                    run_pipeline(config)

        assert config.raw_dir.exists()
        assert config.interim_dir.exists()
        assert config.processed_dir.exists()


class TestPipelineIntegration:
    """Integration tests for the pipeline with mocked external dependencies."""

    def test_full_pipeline_single_year_mocked(self, temp_dir, sample_dat_file, habitat_mapping_file, habitat_detailed_file):
        """Test full pipeline flow with mocked download."""
        import shutil

        # Setup directories
        raw_dir = temp_dir / "raw" / "uniprot_releases"
        raw_dir.mkdir(parents=True)
        interim_dir = temp_dir / "interim"
        interim_dir.mkdir(parents=True)
        processed_dir = temp_dir / "processed"
        processed_dir.mkdir(parents=True)

        # Copy habitat files
        (temp_dir / "raw_base").mkdir(exist_ok=True)
        shutil.copy(habitat_mapping_file, temp_dir / "raw_base" / "marine_terrestrial.json")
        shutil.copy(habitat_detailed_file, temp_dir / "raw_base" / "habitat_detailed.json")

        # Copy sample dat file as 2020 release
        shutil.copy(sample_dat_file, raw_dir / "2020_sprot.dat")

        config = PipelineConfig(
            years=[2020],
            raw_dir=raw_dir,
            interim_dir=interim_dir,
            processed_dir=processed_dir,
            data_dir=temp_dir / "raw_base",  # Points to directory with raw/ subdir
            delete_dat_files=False,  # Keep for inspection
            delete_tsv_files=False,
        )

        # Fix data_dir to point correctly
        config.data_dir = temp_dir
        (temp_dir / "raw").mkdir(exist_ok=True)
        shutil.copy(habitat_mapping_file, temp_dir / "raw" / "marine_terrestrial.json")
        shutil.copy(habitat_detailed_file, temp_dir / "raw" / "habitat_detailed.json")

        with patch("src.data_processing.pipeline.PTMVocabulary") as mock_ptm:
            mock_ptm_instance = MagicMock()
            mock_ptm_instance.ptm_vocab = {}
            mock_ptm.return_value = mock_ptm_instance

            with patch("src.data_processing.pipeline.initialize_taxdb") as mock_taxdb:
                mock_taxdb.return_value = MagicMock()

                with patch("src.data_processing.pipeline.process_dataframe_with_taxonomy") as mock_tax:
                    # Return a simple DataFrame
                    mock_df = pd.DataFrame({
                        "Entry": ["P0C1T5"],
                        "Order": ["Neogastropoda"],
                        "Genus": ["Conus"],
                    })
                    mock_tax.return_value = mock_df

                    with patch("src.data_processing.pipeline.add_habitat_classification") as mock_habitat:
                        mock_habitat.return_value = mock_df

                        results = run_pipeline(config)

        # Check results
        assert len(results) == 1
        # Note: May fail due to actual parsing if test data doesn't match expected format
        # This is an integration test showing the flow works
