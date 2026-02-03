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
        "RC   TISSUE=Venom duct;\n",
        "CC   -!- FUNCTION: Mu-conotoxins block voltage-gated sodium channels.\n",
        "CC   -!- TISSUE SPECIFICITY: Expressed by the venom duct.\n",
        "CC   -!- SIMILARITY: Belongs to the conotoxin M superfamily.\n",
        "KW   Toxin; Ion channel impairing toxin; Voltage-gated sodium channel impairing\n",
        "KW   toxin.\n",
        "FT   SIGNAL          1..22\n",
        "FT   DISULFID        41..77\n",
        'FT                   /note="Interchain"\n',
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
    import json

    tsv_path = temp_dir / "toxprot_2020.tsv"
    # PTM_Features JSON for first entry (uses readable names, no summary)
    ptm_features_1 = json.dumps(
        {
            "Disulfide bond": [{"pos": "41-77", "note": "", "evidence": ""}],
        }
    )
    data = {
        "Entry": ["P0C1T5", "P0C1T6"],
        "Organism": ["Conus magus", "Conus striatus"],
        "Organism (ID)": [6491, 6492],
        "Protein families": ["Conotoxin M superfamily. More info here.", "I1 superfamily"],
        "Length": [85, 90],
        "Fragment": ["", "fragment"],
        "Toxic dose": ["", "LD50=1mg/kg"],
        "PTM_Features": [ptm_features_1, ""],
        "PTM Keywords": ["Disulfide bond", ""],
        "Sequence": ["MKLTCVLVVALLLLVPATTISASGDGRCCKGKRECNNPPCKGKGCSSPKCWPGCC", "MKKKKKKKKK"],
        "Signal peptide (range)": ["1-22", ""],
        "Protein existence": ["1: Evidence at protein level", "2: Evidence at transcript level"],
        "ToxProt definition": ["both", "venom_tissue"],
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
            "marine": ["Neogastropoda", "Octopoda"],
        },
        "ambiguous_orders": {
            "Anguilliformes": {"marine": ["Gymnothorax"], "freshwater": ["Anguilla"]}
        },
    }
    mapping_path = temp_dir / "marine_terrestrial.json"
    with open(mapping_path, "w") as f:
        json.dump(mapping, f)
    return mapping_path


@pytest.fixture
def sample_xml_entry():
    """Minimal UniProt XML entry string for a toxin protein."""
    return """<entry dataset="Swiss-Prot" created="2005-06-07" modified="2024-01-24" version="42" xmlns="http://uniprot.org/uniprot">
  <accession>P0C1T5</accession>
  <accession>P0C1T6</accession>
  <name>TOXB1_CONMA</name>
  <protein>
    <recommendedName>
      <fullName>Mu-conotoxin SmIIIA</fullName>
    </recommendedName>
    <alternativeName>
      <fullName>SmIIIA</fullName>
    </alternativeName>
  </protein>
  <organism>
    <name type="scientific">Conus magus</name>
    <name type="common">Magician's cone</name>
    <dbReference type="NCBI Taxonomy" id="6491"/>
    <lineage>
      <taxon>Eukaryota</taxon>
      <taxon>Metazoa</taxon>
      <taxon>Spiralia</taxon>
      <taxon>Lophotrochozoa</taxon>
      <taxon>Mollusca</taxon>
      <taxon>Gastropoda</taxon>
    </lineage>
  </organism>
  <reference key="1">
    <source>
      <tissue>Venom duct</tissue>
    </source>
  </reference>
  <comment type="function">
    <text>Mu-conotoxins block voltage-gated sodium channels.</text>
  </comment>
  <comment type="tissue specificity">
    <text>Expressed by the venom duct.</text>
  </comment>
  <comment type="similarity">
    <text>Belongs to the conotoxin M superfamily.</text>
  </comment>
  <keyword id="KW-0800">Toxin</keyword>
  <keyword id="KW-0964">Secreted</keyword>
  <keyword id="KW-1015">Disulfide bond</keyword>
  <keyword id="KW-0379">Hydroxylation</keyword>
  <feature type="signal peptide">
    <location>
      <begin position="1"/>
      <end position="22"/>
    </location>
  </feature>
  <feature type="disulfide bond" evidence="1">
    <location>
      <begin position="41"/>
      <end position="77"/>
    </location>
  </feature>
  <feature type="disulfide bond">
    <location>
      <begin position="48"/>
      <end position="63"/>
    </location>
  </feature>
  <feature type="modified residue" description="4-hydroxyproline" evidence="2">
    <location>
      <position position="49"/>
    </location>
  </feature>
  <feature type="modified residue" description="4-hydroxyproline">
    <location>
      <position position="52"/>
    </location>
  </feature>
  <evidence key="1" type="ECO:0000269">
    <source><dbReference type="PubMed" id="12345"/></source>
  </evidence>
  <evidence key="2" type="ECO:0000305"/>
  <dbReference type="InterPro" id="IPR004214"/>
  <dbReference type="Pfam" id="PF02950"/>
  <dbReference type="GO" id="GO:0005576">
    <property type="term" value="C:extracellular region"/>
    <property type="evidence" value="ECO:0000501"/>
  </dbReference>
  <dbReference type="GO" id="GO:0090729">
    <property type="term" value="F:toxin activity"/>
    <property type="evidence" value="ECO:0000501"/>
  </dbReference>
  <proteinExistence type="evidence at protein level"/>
  <sequence length="85" mass="8547" checksum="ABC123" modified="2005-06-07" version="1">
MKLTCVLVVALLLLVPATTISASGDGRCCKGKRECNNPPCKGKGCSSPKCWPGCC
  </sequence>
</entry>"""


@pytest.fixture
def sample_xml_entry_non_metazoa():
    """Minimal UniProt XML entry string for a non-Metazoa organism."""
    return """<entry dataset="Swiss-Prot" xmlns="http://uniprot.org/uniprot">
  <accession>P12345</accession>
  <name>PROT_ECOLI</name>
  <protein>
    <recommendedName>
      <fullName>Test protein</fullName>
    </recommendedName>
  </protein>
  <organism>
    <name type="scientific">Escherichia coli</name>
    <dbReference type="NCBI Taxonomy" id="562"/>
    <lineage>
      <taxon>Bacteria</taxon>
      <taxon>Proteobacteria</taxon>
    </lineage>
  </organism>
  <comment type="function">
    <text>Test function.</text>
  </comment>
  <proteinExistence type="evidence at protein level"/>
  <sequence length="100" mass="10000" checksum="DEF456" modified="2005-01-01" version="1">
MKKKKKKKKKAAAAAAAAAABBBBBBBBBBCCCCCCCCCC
  </sequence>
</entry>"""


@pytest.fixture
def sample_xml_file(temp_dir, sample_xml_entry):
    """Create a sample .xml file wrapping the test entry in a uniprot document."""
    xml_path = temp_dir / "test_sprot.xml"
    content = f"""<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot">
{sample_xml_entry}
</uniprot>"""
    with open(xml_path, "w") as f:
        f.write(content)
    return xml_path


@pytest.fixture
def sample_xml_file_mixed(temp_dir, sample_xml_entry, sample_xml_entry_non_metazoa):
    """Create a .xml file with both matching and non-matching entries."""
    xml_path = temp_dir / "mixed_sprot.xml"
    content = f"""<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot">
{sample_xml_entry}
{sample_xml_entry_non_metazoa}
</uniprot>"""
    with open(xml_path, "w") as f:
        f.write(content)
    return xml_path


@pytest.fixture
def sample_ptmlist_file(temp_dir):
    """Create a minimal ptmlist.txt for testing PTMVocabulary."""
    ptmlist_path = temp_dir / "ptmlist.txt"
    content = """ID   4-hydroxyproline
AC   PTM-0001
FT   MOD_RES
KW   Hydroxylation.
//
ID   Alanine amide
AC   PTM-0057
FT   MOD_RES
KW   Amidation.
//
ID   4-carboxyglutamate
AC   PTM-0100
FT   MOD_RES
KW   Gamma-carboxyglutamic acid.
//
ID   Phosphoserine
AC   PTM-0200
FT   MOD_RES
KW   Phosphoprotein.
//
"""
    with open(ptmlist_path, "w") as f:
        f.write(content)
    return ptmlist_path


@pytest.fixture
def habitat_detailed_file(temp_dir):
    """Create a sample detailed habitat mapping JSON file."""
    import json

    mapping = {
        "clear_orders": {
            "terrestrial": ["Squamata", "Scorpiones", "Araneae"],
            "marine": ["Neogastropoda", "Octopoda"],
            "freshwater": [],
            "estuarine": [],
        },
        "ambiguous_orders": {
            "Anguilliformes": {
                "marine": ["Gymnothorax"],
                "freshwater": ["Anguilla"],
                "terrestrial": [],
                "estuarine": [],
            }
        },
    }
    detailed_path = temp_dir / "habitat_detailed.json"
    with open(detailed_path, "w") as f:
        json.dump(mapping, f)
    return detailed_path
