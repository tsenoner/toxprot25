"""Tests for the Swiss-Prot XML parser module."""

import json
import xml.etree.ElementTree as ET
from pathlib import Path

from src.data_processing.clean_data import extract_ptm_summary
from src.data_processing.parse_sprot_dat import PTMVocabulary
from src.data_processing.parse_sprot_xml import (
    SwissProtXMLParser,
    detect_namespace,
    process_xml_file,
)


def _parse_entry(xml_str):
    """Helper: parse an XML entry string and return (entry_data, meets_criteria)."""
    elem = ET.fromstring(xml_str)
    # Detect namespace from the element tag
    tag = elem.tag
    ns = ""
    if tag.startswith("{"):
        ns_uri = tag[1 : tag.index("}")]
        ns = f"{{{ns_uri}}}"
    parser = SwissProtXMLParser(ns=ns)
    return parser.parse_entry(elem)


class TestSwissProtXMLParser:
    """Tests for the SwissProtXMLParser class."""

    def test_parse_entry_extracts_accession(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert entry_data["Entry"] == "P0C1T5"

    def test_parse_entry_extracts_entry_name(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert entry_data["Entry Name"] == "TOXB1_CONMA"

    def test_parse_entry_extracts_organism(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert "Conus magus" in entry_data["Organism"]

    def test_parse_entry_extracts_taxonomy_id(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert entry_data["Organism (ID)"] == "6491"

    def test_parse_entry_extracts_protein_names(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert "Mu-conotoxin SmIIIA" in entry_data["Protein names"]
        assert "SmIIIA" in entry_data["Protein names"]

    def test_parse_entry_extracts_length(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert entry_data["Length"] == "85"

    def test_parse_entry_extracts_mass(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert entry_data["Mass"] == "8547"

    def test_parse_entry_extracts_function(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert "sodium channels" in entry_data["Function [CC]"]

    def test_parse_entry_extracts_tissue_specificity(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert "venom duct" in entry_data["Tissue specificity"]

    def test_parse_entry_extracts_signal_peptide(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert entry_data["Signal peptide"] == "Yes"
        assert entry_data["Signal peptide (range)"] == "1-22"

    def test_parse_entry_extracts_ptm_features_json(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        ptm_features = json.loads(entry_data["PTM_Features"])

        # Check disulfide bonds (now uses readable names)
        assert "Disulfide bond" in ptm_features
        disulfids = ptm_features["Disulfide bond"]
        positions = [d["pos"] for d in disulfids]
        assert "41-77" in positions
        assert "48-63" in positions

        # Check modified residues (now uses readable names)
        assert "Modified residue" in ptm_features
        mod_res = ptm_features["Modified residue"]
        assert any(m["pos"] == "49" and m["note"] == "4-hydroxyproline" for m in mod_res)
        assert any(m["pos"] == "52" and m["note"] == "4-hydroxyproline" for m in mod_res)

    def test_parse_entry_ptm_features_no_summary(self, sample_xml_entry):
        """PTM_Features JSON should contain raw data only, no summary."""
        entry_data, _ = _parse_entry(sample_xml_entry)
        ptm_features = json.loads(entry_data["PTM_Features"])

        # Summary is now computed in clean_data.py, not stored in PTM_Features
        assert "summary" not in ptm_features
        assert "mod_res_keywords" not in ptm_features

    def test_parse_entry_extracts_interpro(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert "IPR004214" in entry_data["InterPro"]

    def test_parse_entry_extracts_pfam(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert "PF02950" in entry_data["Pfam"]

    def test_parse_entry_extracts_go_terms(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert "GO:0005576" in entry_data["Gene Ontology (GO)"]
        assert "GO:0090729" in entry_data["Gene Ontology (GO)"]
        assert "extracellular region" in entry_data["Gene Ontology (cellular component)"]
        assert "toxin activity" in entry_data["Gene Ontology (molecular function)"]

    def test_parse_entry_extracts_protein_existence(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert "1: Evidence at protein level" in entry_data["Protein existence"]

    def test_parse_entry_extracts_sequence(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert " " not in entry_data["Sequence"]
        assert entry_data["Sequence"].startswith("MKLTCVLVVA")

    def test_parse_entry_meets_criteria_metazoa_venom(self, sample_xml_entry):
        _, meets_criteria = _parse_entry(sample_xml_entry)
        assert meets_criteria is True

    def test_parse_entry_toxprot_definition_both(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        assert entry_data["ToxProt definition"] == "both"

    def test_parse_entry_rejects_non_metazoa(self, sample_xml_entry_non_metazoa):
        _, meets_criteria = _parse_entry(sample_xml_entry_non_metazoa)
        assert meets_criteria is False

    def test_parse_entry_cleans_protein_families(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        # "Belongs to the conotoxin M superfamily." → "Conotoxin M superfamily"
        assert entry_data["Protein families"] == "Conotoxin M superfamily"

    def test_parse_entry_resolves_evidence(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        ptm_features = json.loads(entry_data["PTM_Features"])
        # Disulfide bond at 41-77 has evidence="1" which maps to ECO:0000269
        disulfid_41_77 = next(d for d in ptm_features["Disulfide bond"] if d["pos"] == "41-77")
        assert "ECO:0000269" in disulfid_41_77["evidence"]

    def test_parse_entry_extracts_ptm_keywords(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        ptm_keywords = entry_data["PTM Keywords"]
        # Keywords should be sorted alphabetically and joined with "; "
        assert "Disulfide bond" in ptm_keywords
        assert "Hydroxylation" in ptm_keywords
        # Sorted: Disulfide bond comes before Hydroxylation
        assert ptm_keywords == "Disulfide bond; Hydroxylation"

    def test_parse_entry_ptm_keywords_excludes_non_ptm(self, sample_xml_entry):
        entry_data, _ = _parse_entry(sample_xml_entry)
        ptm_keywords = entry_data["PTM Keywords"]
        # "Toxin" (KW-0800) and "Secreted" (KW-0964) are not PTM keywords
        assert "Toxin" not in ptm_keywords
        assert "Secreted" not in ptm_keywords


class TestPTMVocabularyEnhancements:
    """Tests for PTMVocabulary description_to_keyword mapping."""

    def test_description_to_keyword_mapping(self, sample_ptmlist_file):
        vocab = PTMVocabulary(str(sample_ptmlist_file.parent))
        vocab.ptmlist_path = sample_ptmlist_file
        vocab.load_ptmlist()

        assert vocab.description_to_keyword["4-hydroxyproline"] == "Hydroxylation"
        assert vocab.description_to_keyword["Alanine amide"] == "Amidation"
        assert vocab.description_to_keyword["4-carboxyglutamate"] == "Gamma-carboxyglutamic acid"
        assert vocab.description_to_keyword["Phosphoserine"] == "Phosphoprotein"

    def test_resolve_mod_res_keyword(self, sample_ptmlist_file):
        vocab = PTMVocabulary(str(sample_ptmlist_file.parent))
        vocab.ptmlist_path = sample_ptmlist_file
        vocab.load_ptmlist()

        assert vocab.resolve_mod_res_keyword("4-hydroxyproline") == "Hydroxylation"
        assert vocab.resolve_mod_res_keyword("Alanine amide") == "Amidation"

    def test_resolve_mod_res_keyword_strips_qualifier(self, sample_ptmlist_file):
        vocab = PTMVocabulary(str(sample_ptmlist_file.parent))
        vocab.ptmlist_path = sample_ptmlist_file
        vocab.load_ptmlist()

        assert vocab.resolve_mod_res_keyword("Phosphoserine; by host") == "Phosphoprotein"

    def test_unknown_description_returns_other(self, sample_ptmlist_file):
        vocab = PTMVocabulary(str(sample_ptmlist_file.parent))
        vocab.ptmlist_path = sample_ptmlist_file
        vocab.load_ptmlist()

        assert vocab.resolve_mod_res_keyword("UnknownModification") == "Other"


class TestCleanDataPTMExtraction:
    """Tests for PTM summary extraction functions in clean_data.py."""

    def test_extract_ptm_summary_basic(self):
        ptm_json = json.dumps({
            "Disulfide bond": [{"pos": "41-77", "note": "", "evidence": ""}],
            "Modified residue": [
                {"pos": "49", "note": "4-hydroxyproline", "evidence": ""},
                {"pos": "52", "note": "4-hydroxyproline", "evidence": ""},
            ],
        })
        summary = extract_ptm_summary(ptm_json)
        assert "Disulfide bond:1" in summary
        # Modified residue notes are resolved to keywords
        assert "Hydroxylation:2" in summary
        assert "Modified residue" not in summary

    def test_extract_ptm_summary_empty(self):
        assert extract_ptm_summary("") == ""
        assert extract_ptm_summary(None) == ""

    def test_extract_ptm_summary_invalid_json(self):
        assert extract_ptm_summary("not valid json") == ""


class TestNamespaceDetection:
    """Tests for XML namespace detection."""

    def test_namespace_detection_http(self, sample_xml_file):
        ns = detect_namespace(sample_xml_file)
        assert ns == "{http://uniprot.org/uniprot}"

    def test_namespace_detection_https(self, temp_dir):
        xml_path = temp_dir / "https_test.xml"
        content = '''<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="https://uniprot.org/uniprot">
  <entry dataset="Swiss-Prot">
    <accession>P00001</accession>
    <name>TEST</name>
    <sequence length="10" mass="1000">MAAAAA</sequence>
  </entry>
</uniprot>'''
        with open(xml_path, "w") as f:
            f.write(content)
        ns = detect_namespace(xml_path)
        assert ns == "{https://uniprot.org/uniprot}"


class TestProcessXMLFile:
    """Tests for the process_xml_file function."""

    def test_process_xml_creates_output_file(self, sample_xml_file, temp_dir):
        output_file = temp_dir / "output.tsv"
        result = process_xml_file(sample_xml_file, output_file)
        assert result is True
        assert output_file.exists()

    def test_process_xml_output_has_header(self, sample_xml_file, temp_dir):
        output_file = temp_dir / "output.tsv"
        process_xml_file(sample_xml_file, output_file)

        with open(output_file) as f:
            header = f.readline().strip()

        assert "Entry" in header
        assert "Organism" in header
        assert "Sequence" in header
        assert "PTM_Features" in header

    def test_process_xml_filters_entries(self, sample_xml_file_mixed, temp_dir):
        output_file = temp_dir / "output.tsv"
        process_xml_file(sample_xml_file_mixed, output_file)

        with open(output_file) as f:
            lines = f.readlines()
        # header + 1 matching entry (non-Metazoa rejected)
        assert len(lines) == 2

    def test_process_xml_returns_false_for_missing_file(self, temp_dir):
        result = process_xml_file(
            temp_dir / "nonexistent.xml",
            temp_dir / "output.tsv",
        )
        assert result is False

    def test_process_xml_creates_output_directory(self, sample_xml_file, temp_dir):
        output_file = temp_dir / "subdir" / "output.tsv"
        result = process_xml_file(sample_xml_file, output_file)
        assert result is True
        assert output_file.exists()

    def test_process_xml_with_ptm_features(self, sample_xml_file, temp_dir):
        output_file = temp_dir / "output.tsv"
        result = process_xml_file(sample_xml_file, output_file)
        assert result is True

        with open(output_file) as f:
            header = f.readline().strip().split("\t")
            data = f.readline().strip().split("\t")

        ptm_features_idx = header.index("PTM_Features")
        ptm_features = json.loads(data[ptm_features_idx])
        # Uses readable names now
        assert "Disulfide bond" in ptm_features
        assert "Modified residue" in ptm_features
        assert len(ptm_features["Disulfide bond"]) == 2
        assert len(ptm_features["Modified residue"]) == 2
