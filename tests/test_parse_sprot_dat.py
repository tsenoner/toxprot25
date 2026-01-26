"""Tests for the Swiss-Prot parser module."""

from pathlib import Path

from src.data_processing.parse_sprot_dat import (
    SwissProtParser,
    get_output_filename,
    process_swissprot_file,
)


class TestSwissProtParser:
    """Tests for the SwissProtParser class."""

    def test_parse_position_modern_format(self):
        """Test parsing modern '..' format positions."""
        parser = SwissProtParser()
        assert parser.parse_position("41..77") == "41-77"
        assert parser.parse_position("1..22") == "1-22"
        assert parser.parse_position("100..200") == "100-200"

    def test_parse_position_same_start_end(self):
        """Test that same start/end returns single position."""
        parser = SwissProtParser()
        assert parser.parse_position("41..41") == "41"
        assert parser.parse_position("100..100") == "100"

    def test_parse_position_with_question_marks(self):
        """Test parsing positions with unknown values."""
        parser = SwissProtParser()
        assert parser.parse_position("?..77") == "?-77"
        assert parser.parse_position("41..?") == "41-?"
        # When both are ?, they're treated as same position
        assert parser.parse_position("?..?") == "?"

    def test_parse_position_space_separated(self):
        """Test parsing old space-separated format."""
        parser = SwissProtParser()
        assert parser.parse_position("41     77") == "41-77"
        assert parser.parse_position("1 22") == "1-22"

    def test_parse_position_empty(self):
        """Test parsing empty position string."""
        parser = SwissProtParser()
        assert parser.parse_position("") == ""

    def test_parse_entry_extracts_accession(self, sample_swissprot_entry):
        """Test that entry accession is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert entry_data["Entry"] == "P0C1T5"

    def test_parse_entry_extracts_entry_name(self, sample_swissprot_entry):
        """Test that entry name is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert entry_data["Entry Name"] == "TOXB1_CONMA"

    def test_parse_entry_extracts_organism(self, sample_swissprot_entry):
        """Test that organism is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert "Conus magus" in entry_data["Organism"]

    def test_parse_entry_extracts_taxonomy_id(self, sample_swissprot_entry):
        """Test that taxonomy ID is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert entry_data["Organism (ID)"] == "6491"

    def test_parse_entry_extracts_protein_name(self, sample_swissprot_entry):
        """Test that protein name is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert "Mu-conotoxin SmIIIA" in entry_data["Protein names"]

    def test_parse_entry_extracts_length(self, sample_swissprot_entry):
        """Test that sequence length is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert entry_data["Length"] == "85"

    def test_parse_entry_extracts_function(self, sample_swissprot_entry):
        """Test that function annotation is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert "sodium channels" in entry_data["Function [CC]"]

    def test_parse_entry_extracts_tissue_specificity(self, sample_swissprot_entry):
        """Test that tissue specificity is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert "venom duct" in entry_data["Tissue specificity"]

    def test_parse_entry_extracts_signal_peptide(self, sample_swissprot_entry):
        """Test that signal peptide is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert entry_data["Signal peptide"] == "Yes"
        assert entry_data["Signal peptide (range)"] == "1-22"

    def test_parse_entry_extracts_disulfide_bonds(self, sample_swissprot_entry):
        """Test that disulfide bonds are extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert "41-77" in entry_data["Disulfide bond"]
        assert "48-63" in entry_data["Disulfide bond"]

    def test_parse_entry_extracts_interpro(self, sample_swissprot_entry):
        """Test that InterPro IDs are extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert "IPR004214" in entry_data["InterPro"]

    def test_parse_entry_extracts_pfam(self, sample_swissprot_entry):
        """Test that Pfam IDs are extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert "PF02950" in entry_data["Pfam"]

    def test_parse_entry_extracts_go_terms(self, sample_swissprot_entry):
        """Test that GO terms are extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        assert "GO:0005576" in entry_data["Gene Ontology (GO)"]
        assert "GO:0090729" in entry_data["Gene Ontology (GO)"]

    def test_parse_entry_extracts_sequence(self, sample_swissprot_entry):
        """Test that sequence is extracted correctly."""
        parser = SwissProtParser()
        entry_data, _ = parser.parse_entry(sample_swissprot_entry)
        # Sequence should have no spaces
        assert " " not in entry_data["Sequence"]
        assert entry_data["Sequence"].startswith("MKLTCVLVVA")

    def test_parse_entry_meets_criteria_metazoa_venom(self, sample_swissprot_entry):
        """Test that Metazoa + venom tissue meets criteria."""
        parser = SwissProtParser()
        _, meets_criteria = parser.parse_entry(sample_swissprot_entry)
        assert meets_criteria is True

    def test_parse_entry_rejects_non_metazoa(self, sample_swissprot_entry_non_metazoa):
        """Test that non-Metazoa entry is rejected."""
        parser = SwissProtParser()
        _, meets_criteria = parser.parse_entry(sample_swissprot_entry_non_metazoa)
        assert meets_criteria is False


class TestProcessSwissProtFile:
    """Tests for the process_swissprot_file function."""

    def test_process_creates_output_file(self, sample_dat_file, temp_dir):
        """Test that processing creates output TSV file."""
        output_file = temp_dir / "output.tsv"
        result = process_swissprot_file(sample_dat_file, output_file)
        assert result is True
        assert output_file.exists()

    def test_process_output_has_header(self, sample_dat_file, temp_dir):
        """Test that output file has correct header."""
        output_file = temp_dir / "output.tsv"
        process_swissprot_file(sample_dat_file, output_file)

        with open(output_file) as f:
            header = f.readline().strip()

        assert "Entry" in header
        assert "Organism" in header
        assert "Sequence" in header

    def test_process_filters_entries(self, temp_dir, sample_swissprot_entry, sample_swissprot_entry_non_metazoa):
        """Test that only matching entries are included in output."""
        # Create a file with both matching and non-matching entries
        dat_path = temp_dir / "mixed.dat"
        with open(dat_path, "w") as f:
            for line in sample_swissprot_entry:
                f.write(line)
            f.write("//\n")
            for line in sample_swissprot_entry_non_metazoa:
                f.write(line)
            f.write("//\n")

        output_file = temp_dir / "output.tsv"
        process_swissprot_file(dat_path, output_file)

        # Read output and count entries (excluding header)
        with open(output_file) as f:
            lines = f.readlines()
        assert len(lines) == 2  # header + 1 matching entry

    def test_process_returns_false_for_missing_file(self, temp_dir):
        """Test that processing returns False for non-existent input."""
        result = process_swissprot_file(
            temp_dir / "nonexistent.dat",
            temp_dir / "output.tsv"
        )
        assert result is False

    def test_process_creates_output_directory(self, sample_dat_file, temp_dir):
        """Test that output directory is created if needed."""
        output_file = temp_dir / "subdir" / "output.tsv"
        result = process_swissprot_file(sample_dat_file, output_file)
        assert result is True
        assert output_file.exists()


class TestGetOutputFilename:
    """Tests for the get_output_filename helper function."""

    def test_extracts_year_from_filename(self):
        """Test extracting year from standard filename."""
        result = get_output_filename(Path("2005_sprot.dat"))
        assert result == "toxprot_2005.tsv"

    def test_handles_different_years(self):
        """Test handling various years."""
        assert get_output_filename(Path("2020_sprot.dat")) == "toxprot_2020.tsv"
        assert get_output_filename(Path("2025_sprot.dat")) == "toxprot_2025.tsv"

    def test_handles_path_with_directory(self):
        """Test handling full path with directory."""
        result = get_output_filename(Path("/data/raw/2015_sprot.dat"))
        assert result == "toxprot_2015.tsv"

    def test_fallback_for_non_standard_name(self):
        """Test fallback for non-standard filename."""
        result = get_output_filename(Path("custom_name.dat"))
        assert result == "toxprot_custom_name.tsv"
