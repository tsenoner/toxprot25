"""Tests for the data cleaning module."""

from pathlib import Path

import pandas as pd

from src.data_processing.clean_data import (
    create_fasta_file,
    determine_habitat,
    determine_habitat_detailed,
    get_year_from_filename,
    update_protfams,
)


class TestUpdateProtfams:
    """Tests for the update_protfams function."""

    def test_splits_on_period(self):
        """Test that family names are split on period."""
        df = pd.DataFrame({"Protein families": ["Conotoxin M superfamily. More info."]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Conotoxin M superfamily"

    def test_splits_on_comma(self):
        """Test that family names are split on comma."""
        df = pd.DataFrame({"Protein families": ["Conotoxin M superfamily, subfamily A"]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Conotoxin M superfamily"

    def test_splits_on_semicolon(self):
        """Test that family names are split on semicolon."""
        df = pd.DataFrame({"Protein families": ["Conotoxin M superfamily; more info"]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Conotoxin M superfamily"

    def test_corrects_i1_superfamily(self):
        """Test I1 superfamily correction."""
        df = pd.DataFrame({"Protein families": ["I1 superfamily"]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Conotoxin I1 superfamily"

    def test_corrects_o1_superfamily(self):
        """Test O1 superfamily correction."""
        df = pd.DataFrame({"Protein families": ["O1 superfamily"]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Conotoxin O1 superfamily"

    def test_corrects_o2_superfamily(self):
        """Test O2 superfamily correction."""
        df = pd.DataFrame({"Protein families": ["O2 superfamily"]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Conotoxin O2 superfamily"

    def test_corrects_conotoxin_m_family(self):
        """Test Conotoxin M family to superfamily correction."""
        df = pd.DataFrame({"Protein families": ["Conotoxin M family"]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Conotoxin M superfamily"

    def test_corrects_bradykinin_hyphen(self):
        """Test Bradykinin-potentiating peptide family correction."""
        df = pd.DataFrame({"Protein families": ["Bradykinin- potentiating peptide family"]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Bradykinin-potentiating peptide family"

    def test_preserves_correct_names(self):
        """Test that already correct names are preserved."""
        df = pd.DataFrame({"Protein families": ["Conotoxin M superfamily"]})
        result = update_protfams(df)
        assert result["Protein families"].iloc[0] == "Conotoxin M superfamily"

    def test_handles_empty_values(self):
        """Test handling of empty/NaN values."""
        df = pd.DataFrame({"Protein families": ["", pd.NA, "Conotoxin M superfamily"]})
        result = update_protfams(df)
        assert (
            pd.isna(result["Protein families"].iloc[0]) or result["Protein families"].iloc[0] == ""
        )


class TestCreateFastaFile:
    """Tests for the create_fasta_file function."""

    def test_creates_fasta_file(self, temp_dir):
        """Test that FASTA file is created."""
        df = pd.DataFrame({"Entry": ["P0C1T5"], "Sequence": ["MKLTCVLVVALLLLVPATTI"]})
        fasta_path = temp_dir / "test.fasta"
        create_fasta_file(df, "Entry", "Sequence", fasta_path)
        assert fasta_path.exists()

    def test_fasta_format_correct(self, temp_dir):
        """Test that FASTA format is correct."""
        df = pd.DataFrame({"Entry": ["P0C1T5"], "Sequence": ["MKLTCVLVVALLLLVPATTI"]})
        fasta_path = temp_dir / "test.fasta"
        create_fasta_file(df, "Entry", "Sequence", fasta_path)

        content = fasta_path.read_text()
        assert content.startswith(">P0C1T5")
        assert "MKLTCVLVVALLLLVPATTI" in content

    def test_removes_signal_peptide(self, temp_dir):
        """Test that signal peptide is removed when range provided."""
        df = pd.DataFrame(
            {
                "Entry": ["P0C1T5"],
                "Sequence": ["MKLTCVLVVALLLLVPATTIREST"],  # 24 chars total
                "Signal peptide (range)": ["1-10"],  # Remove first 10 amino acids
            }
        )
        fasta_path = temp_dir / "test.fasta"
        create_fasta_file(df, "Entry", "Sequence", fasta_path, "Signal peptide (range)")

        content = fasta_path.read_text()
        # Sequence should start after position 10
        assert "LLLLVPATTIREST" in content

    def test_preserves_sequence_without_signal(self, temp_dir):
        """Test that sequence is preserved when no signal peptide range."""
        df = pd.DataFrame(
            {
                "Entry": ["P0C1T5"],
                "Sequence": ["MKLTCVLVVALLLLVPATTI"],
                "Signal peptide (range)": [""],
            }
        )
        fasta_path = temp_dir / "test.fasta"
        create_fasta_file(df, "Entry", "Sequence", fasta_path, "Signal peptide (range)")

        content = fasta_path.read_text()
        assert "MKLTCVLVVALLLLVPATTI" in content

    def test_skips_empty_sequences(self, temp_dir):
        """Test that empty sequences are skipped."""
        df = pd.DataFrame({"Entry": ["P0C1T5", "P0C1T6"], "Sequence": ["MKLTCVLVVALLLLVPATTI", ""]})
        fasta_path = temp_dir / "test.fasta"
        create_fasta_file(df, "Entry", "Sequence", fasta_path)

        content = fasta_path.read_text()
        assert ">P0C1T5" in content
        assert ">P0C1T6" not in content

    def test_creates_parent_directory(self, temp_dir):
        """Test that parent directory is created if needed."""
        df = pd.DataFrame({"Entry": ["P0C1T5"], "Sequence": ["MKLTCVLVVALLLLVPATTI"]})
        fasta_path = temp_dir / "subdir" / "test.fasta"
        create_fasta_file(df, "Entry", "Sequence", fasta_path)
        assert fasta_path.exists()


class TestDetermineHabitat:
    """Tests for the determine_habitat function."""

    def test_clear_terrestrial_order(self, habitat_mapping_file):
        """Test clear terrestrial order classification."""
        import json

        with open(habitat_mapping_file) as f:
            mapping = json.load(f)

        row = {"Order": "Squamata", "Genus": "Naja"}
        assert determine_habitat(row, mapping) == "terrestrial"

    def test_clear_marine_order(self, habitat_mapping_file):
        """Test clear marine order classification."""
        import json

        with open(habitat_mapping_file) as f:
            mapping = json.load(f)

        row = {"Order": "Neogastropoda", "Genus": "Conus"}
        assert determine_habitat(row, mapping) == "marine"

    def test_ambiguous_order_marine_genus(self, habitat_mapping_file):
        """Test ambiguous order with marine genus."""
        import json

        with open(habitat_mapping_file) as f:
            mapping = json.load(f)

        row = {"Order": "Anguilliformes", "Genus": "Gymnothorax"}
        assert determine_habitat(row, mapping) == "marine"

    def test_unknown_order(self, habitat_mapping_file):
        """Test unknown order returns 'unknown'."""
        import json

        with open(habitat_mapping_file) as f:
            mapping = json.load(f)

        row = {"Order": "UnknownOrder", "Genus": "UnknownGenus"}
        assert determine_habitat(row, mapping) == "unknown"

    def test_missing_genus(self, habitat_mapping_file):
        """Test handling of missing genus."""
        import json

        with open(habitat_mapping_file) as f:
            mapping = json.load(f)

        row = {"Order": "Squamata"}  # No Genus key
        assert determine_habitat(row, mapping) == "terrestrial"


class TestDetermineHabitatDetailed:
    """Tests for the determine_habitat_detailed function."""

    def test_clear_terrestrial_order(self, habitat_detailed_file):
        """Test clear terrestrial order classification."""
        import json

        with open(habitat_detailed_file) as f:
            mapping = json.load(f)

        row = {"Order": "Squamata", "Genus": "Naja"}
        assert determine_habitat_detailed(row, mapping) == "terrestrial"

    def test_clear_marine_order(self, habitat_detailed_file):
        """Test clear marine order classification."""
        import json

        with open(habitat_detailed_file) as f:
            mapping = json.load(f)

        row = {"Order": "Neogastropoda", "Genus": "Conus"}
        assert determine_habitat_detailed(row, mapping) == "marine"

    def test_ambiguous_order_freshwater_genus(self, habitat_detailed_file):
        """Test ambiguous order with freshwater genus."""
        import json

        with open(habitat_detailed_file) as f:
            mapping = json.load(f)

        row = {"Order": "Anguilliformes", "Genus": "Anguilla"}
        assert determine_habitat_detailed(row, mapping) == "freshwater"


class TestGetYearFromFilename:
    """Tests for the get_year_from_filename function."""

    def test_extracts_year(self):
        """Test extracting year from filename."""
        assert get_year_from_filename(Path("toxprot_2020.tsv")) == "2020"
        assert get_year_from_filename(Path("toxprot_2005.tsv")) == "2005"
        assert get_year_from_filename(Path("toxprot_2025.tsv")) == "2025"

    def test_handles_full_path(self):
        """Test handling full path."""
        assert get_year_from_filename(Path("/data/interim/toxprot_2020.tsv")) == "2020"

    def test_fallback_for_non_standard_name(self):
        """Test fallback for non-standard filename."""
        assert get_year_from_filename(Path("custom_name.tsv")) == "custom_name"
