"""Pytest configuration and shared fixtures for ToxProt tests."""

import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_swissprot_entry():
    """Sample Swiss-Prot entry lines for testing parser."""
    return [
        "ID   TOXB1_CONMA             Reviewed;          85 AA.\n",
        "AC   P0C1T5;\n",
        "DT   07-JUN-2005, integrated into UniProtKB/Swiss-Prot.\n",
        "DE   RecName: Full=Mu-conotoxin SmIIIA;\n",
        "DE   AltName: Full=SmIIIA;\n",
        "DE   Flags: Precursor;\n",
        "OS   Conus magus (Magician's cone).\n",
        "OC   Eukaryota; Metazoa; Spiralia; Lophotrochozoa; Mollusca; Gastropoda;\n",
        "OC   Caenogastropoda; Hypsogastropoda; Neogastropoda; Conoidea; Conidae;\n",
        "OC   Conus.\n",
        "OX   NCBI_TaxID=6491;\n",
        "CC   -!- FUNCTION: Mu-conotoxins block voltage-gated sodium channels.\n",
        "CC   -!- TISSUE SPECIFICITY: Expressed by the venom duct.\n",
        "CC   -!- SIMILARITY: Belongs to the conotoxin M superfamily.\n",
        "KW   Toxin; Ion channel impairing toxin; Voltage-gated sodium channel impairing\n",
        "KW   toxin.\n",
        "FT   SIGNAL          1..22\n",
        "FT   DISULFID        41..77\n",
        "FT                   /note=\"Interchain\"\n",
        "FT   DISULFID        48..63\n",
        "PE   1: Evidence at protein level;\n",
        "DR   InterPro; IPR004214; Conotoxin.\n",
        "DR   Pfam; PF02950; Toxin_2; 1.\n",
        "DR   GO; GO:0005576; C:extracellular region; IEA:UniProtKB-SubCell.\n",
        "DR   GO; GO:0090729; F:toxin activity; IEA:UniProtKB-KW.\n",
        "SQ   SEQUENCE   85 AA;  8547 MW;  ABC123 CRC64;\n",
        "     MKLTCVLVVA LLLLVPATTI SASGDGRCCK GKRECNNPPC KGKGCSSPKC\n",
        "     WPGCC\n",
    ]


@pytest.fixture
def sample_swissprot_entry_non_metazoa():
    """Sample Swiss-Prot entry for non-Metazoa organism."""
    return [
        "ID   PROT_ECOLI             Reviewed;          100 AA.\n",
        "AC   P12345;\n",
        "DE   RecName: Full=Test protein;\n",
        "OS   Escherichia coli.\n",
        "OC   Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales;\n",
        "OC   Enterobacteriaceae; Escherichia.\n",
        "OX   NCBI_TaxID=562;\n",
        "CC   -!- FUNCTION: Test function.\n",
        "SQ   SEQUENCE   100 AA;  10000 MW;  ABC123 CRC64;\n",
        "     MKKKKKKKKK AAAAAAAAAA BBBBBBBBBB CCCCCCCCCC\n",
    ]


@pytest.fixture
def sample_dat_file(temp_dir, sample_swissprot_entry):
    """Create a sample .dat file with test entries."""
    dat_path = temp_dir / "test_sprot.dat"
    with open(dat_path, "w") as f:
        for line in sample_swissprot_entry:
            f.write(line)
        f.write("//\n")
    return dat_path


@pytest.fixture
def sample_tsv_file(temp_dir):
    """Create a sample toxprot TSV file for testing clean_data."""
    tsv_path = temp_dir / "toxprot_2020.tsv"
    data = {
        "Entry": ["P0C1T5", "P0C1T6"],
        "Organism": ["Conus magus", "Conus striatus"],
        "Organism (ID)": [6491, 6492],
        "Protein families": ["Conotoxin M superfamily. More info here.", "I1 superfamily"],
        "Length": [85, 90],
        "Fragment": ["", "fragment"],
        "Toxic dose": ["", "LD50=1mg/kg"],
        "Post-translational modification": ["Disulfide bond", ""],
        "PTM Summary": ["DISULFID:2", ""],
        "Sequence": ["MKLTCVLVVALLLLVPATTISASGDGRCCKGKRECNNPPCKGKGCSSPKCWPGCC", "MKKKKKKKKK"],
        "Signal peptide (range)": ["1-22", ""],
        "Protein existence": ["1: Evidence at protein level", "2: Evidence at transcript level"],
        "Gene Ontology (GO)": ["GO:0005576", "GO:0005576"],
        "Gene Ontology (biological process)": ["", ""],
        "Gene Ontology (cellular component)": ["GO:0005576 (extracellular region)", ""],
        "Gene Ontology (molecular function)": ["GO:0090729 (toxin activity)", ""],
    }
    df = pd.DataFrame(data)
    df.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


@pytest.fixture
def sample_csv_file(temp_dir):
    """Create a sample processed CSV file for habitat testing."""
    csv_path = temp_dir / "toxprot_2020.csv"
    data = {
        "Entry": ["P0C1T5", "P0C1T6", "P0C1T7"],
        "Organism (ID)": [6491, 6492, 6493],
        "Protein families": ["Conotoxin M superfamily", "Conotoxin I1 superfamily", "Snake toxin"],
        "Order": ["Neogastropoda", "Neogastropoda", "Squamata"],
        "Genus": ["Conus", "Conus", "Naja"],
    }
    df = pd.DataFrame(data)
    df.to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def habitat_mapping_file(temp_dir):
    """Create a sample habitat mapping JSON file."""
    import json
    mapping = {
        "clear_orders": {
            "terrestrial": ["Squamata", "Scorpiones", "Araneae"],
            "marine": ["Neogastropoda", "Octopoda"]
        },
        "ambiguous_orders": {
            "Anguilliformes": {
                "marine": ["Gymnothorax"],
                "freshwater": ["Anguilla"]
            }
        }
    }
    mapping_path = temp_dir / "marine_terrestrial.json"
    with open(mapping_path, "w") as f:
        json.dump(mapping, f)
    return mapping_path


@pytest.fixture
def habitat_detailed_file(temp_dir):
    """Create a sample detailed habitat mapping JSON file."""
    import json
    mapping = {
        "clear_orders": {
            "terrestrial": ["Squamata", "Scorpiones", "Araneae"],
            "marine": ["Neogastropoda", "Octopoda"],
            "freshwater": [],
            "estuarine": []
        },
        "ambiguous_orders": {
            "Anguilliformes": {
                "marine": ["Gymnothorax"],
                "freshwater": ["Anguilla"],
                "terrestrial": [],
                "estuarine": []
            }
        }
    }
    detailed_path = temp_dir / "habitat_detailed.json"
    with open(detailed_path, "w") as f:
        json.dump(mapping, f)
    return detailed_path
