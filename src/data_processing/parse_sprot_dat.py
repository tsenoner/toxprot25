#!/usr/bin/env python3
"""
UniProtKB/SwissProt parser for Tox-Prot (animal toxin annotation project)


Features:
- Extracts protein metadata, PTMs, signal peptides, and functional annotations
- Handles PTM vocabulary (ptmlist.txt) automatically
- Outputs results in TSV format
"""

import argparse
import logging
import re
import sys
import urllib.request
from pathlib import Path

from tqdm import tqdm

# Module logger - configuration handled by caller or main()
logger = logging.getLogger(__name__)


class PTMVocabulary:
    """Manages PTM controlled vocabulary from ptmlist.txt"""

    PTMLIST_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt"

    PTM_FEATURES = {
        "MOD_RES": "Modified residue",
        "CARBOHYD": "Glycosylation",
        "DISULFID": "Disulfide bond",
        "CROSSLNK": "Cross-link",
        "LIPID": "Lipidation",
    }

    def __init__(self, data_dir="data/raw"):
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.ptmlist_path = self.data_dir / "ptmlist.txt"
        self.ptm_vocab = {}

    def ensure_ptmlist_available(self):
        """Download ptmlist.txt if it doesn't exist"""
        if self.ptmlist_path.exists():
            return True

        try:
            logger.info(f"Downloading ptmlist.txt from {self.PTMLIST_URL}...")
            urllib.request.urlretrieve(self.PTMLIST_URL, self.ptmlist_path)
            logger.info(f"Successfully downloaded to {self.ptmlist_path}")
            return True
        except Exception as e:
            logger.warning(f"Could not download ptmlist.txt: {e}")
            logger.warning("The parser will continue without PTM vocabulary reference")
            logger.warning(f"You can manually download from: {self.PTMLIST_URL}")
            return False

    def load_ptmlist(self):
        """Load and parse ptmlist.txt"""
        if not self.ptmlist_path.exists():
            logger.warning("ptmlist.txt not available, skipping vocabulary loading")
            return

        try:
            with open(self.ptmlist_path, encoding="utf-8") as f:
                current_entry = {}
                for line in f:
                    line = line.rstrip()
                    if line.startswith("ID   "):
                        if current_entry:
                            self.ptm_vocab[current_entry.get("ID")] = current_entry
                        current_entry = {"ID": line[5:].strip()}
                    elif line.startswith("AC   "):
                        current_entry["AC"] = line[5:].strip()
                    elif line.startswith("FT   "):
                        current_entry["FT"] = line[5:].strip()
                    elif line.startswith("//"):
                        if current_entry:
                            self.ptm_vocab[current_entry.get("ID")] = current_entry
                            current_entry = {}
            logger.info(
                f"Loaded {len(self.ptm_vocab)} PTM definitions from {self.ptmlist_path}"
            )
        except Exception as e:
            logger.error(f"Error loading ptmlist.txt: {e}")


class SwissProtParser:
    """SwissProt parser for extracting protein information and features"""

    # Define header fields
    HEADERS = [
        "Entry",
        "Entry Name",
        "Organism",
        "Organism (ID)",
        "Protein names",
        "Protein families",
        "Length",
        "Mass",
        "Mass spectrometry",
        "Fragment",
        "Function [CC]",
        "Tissue specificity",
        "Toxic dose",
        "Disulfide bond",
        "Glycosylation",
        "Modified residue",
        "Cross-link",
        "Lipidation",
        "PTM Summary",
        "Post-translational modification",
        "Signal peptide",
        "Signal peptide (range)",
        "Propeptide",
        "InterPro",
        "Pfam",
        "KEGG",
        "Gene Ontology (GO)",
        "Gene Ontology (biological process)",
        "Gene Ontology (cellular component)",
        "Gene Ontology (molecular function)",
        "Protein existence",
        "ToxProt definition",
        "Sequence",
    ]

    def __init__(self, ptm_vocab=None):
        self.ptm_vocab = ptm_vocab

    def parse_position(self, position_str):
        """Parse position information from feature line"""
        # Handle modern ".." format (e.g., "41..77")
        if ".." in position_str:
            match = re.match(r"(\d+|\?)\.\.(\d+|\?)", position_str)
            if match:
                start = match.group(1) if match.group(1) != "?" else "?"
                end = match.group(2) if match.group(2) != "?" else "?"
                # Return single position for interchain bonds (same start/end)
                if start == end:
                    return start
                return f"{start}-{end}"

        # Handle space-separated format (e.g., "41     77" or "41 77")
        parts = position_str.split()
        if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            # Return single position for interchain bonds (same start/end)
            if parts[0] == parts[1]:
                return parts[0]
            return f"{parts[0]}-{parts[1]}"
        elif len(parts) >= 1:
            match = re.match(r"(\d+)", parts[0])
            if match:
                return match.group(1)

        return ""

    def extract_ptm_details(self, lines, start_idx):
        """Extract PTM note and evidence from feature continuation lines"""
        note = ""
        evidence_list = []
        i = start_idx + 1

        # Process all continuation lines for this feature
        while i < len(lines) and lines[i].startswith("FT                   "):
            line_content = lines[i][21:]  # Get content after "FT                   "

            # Extract note - format: /note="content..."
            if line_content.strip().startswith('/note="'):
                # Start extracting note content (without the '/note="' prefix)
                note_parts = [line_content.strip()[7:]]  # Remove '/note="'

                # Check if note ends on this line
                if not note_parts[0].rstrip().endswith('"'):
                    # Multi-line note - keep reading
                    i += 1
                    while i < len(lines) and lines[i].startswith(
                        "FT                   "
                    ):
                        next_content = lines[i][21:].strip()
                        # Check if this is a new qualifier
                        if next_content.startswith("/"):
                            # Don't consume this line, it's a new qualifier
                            i -= 1
                            break
                        note_parts.append(next_content)
                        # Check if note ends on this line
                        if next_content.rstrip().endswith('"'):
                            break
                        i += 1

                # Join all parts and remove trailing quote
                note = " ".join(note_parts).rstrip().rstrip('"')

            # Extract evidence - format: /evidence="ECO:..."
            elif line_content.strip().startswith('/evidence="'):
                # Start extracting evidence content (without the '/evidence="' prefix)
                evidence_parts = [line_content.strip()[11:]]  # Remove '/evidence="'

                # Check if evidence ends on this line
                if not evidence_parts[0].rstrip().endswith('"'):
                    # Multi-line evidence - keep reading
                    i += 1
                    while i < len(lines) and lines[i].startswith(
                        "FT                   "
                    ):
                        next_content = lines[i][21:].strip()
                        # Check if this is a new qualifier
                        if next_content.startswith("/"):
                            # Don't consume this line, it's a new qualifier
                            i -= 1
                            break
                        evidence_parts.append(next_content)
                        # Check if evidence ends on this line
                        if next_content.rstrip().endswith('"'):
                            break
                        i += 1

                # Join all parts, remove trailing quote, and add to list
                evidence_str = "".join(evidence_parts).rstrip().rstrip('"')
                if evidence_str:
                    evidence_list.append(evidence_str)

            i += 1

        # Combine multiple evidence codes
        evidence = ", ".join(evidence_list) if evidence_list else ""

        return note, evidence

    def parse_entry(self, entry_lines):
        """Parse a SwissProt entry and extract all protein information and features"""
        entry_data = {field: "" for field in self.HEADERS}

        # Initialize variables
        current_section = None
        is_metazoa = False
        has_kw0800_toxin = False
        go_terms = {"P": [], "C": [], "F": []}
        processing_de_flags = False
        in_sequence_block = False

        # PTM storage
        ptm_features = {
            "MOD_RES": [],
            "CARBOHYD": [],
            "DISULFID": [],
            "CROSSLNK": [],
            "LIPID": [],
        }

        # Process each line
        for idx, line in enumerate(entry_lines):
            content = line[5:].strip()

            # Entry ID and length
            if line.startswith("ID"):
                in_sequence_block = False
                match = re.search(r"^(\S+)", content)
                if match:
                    entry_data["Entry Name"] = match.group(1)
                length_match = re.search(r"(\d+) AA", content)
                if length_match:
                    entry_data["Length"] = length_match.group(1)

            # Accession
            elif line.startswith("AC"):
                in_sequence_block = False
                if not entry_data["Entry"]:
                    primary_acc = content.split(";")[0].strip()
                    entry_data["Entry"] = primary_acc

            # Organism
            elif line.startswith("OS"):
                in_sequence_block = False
                if not entry_data["Organism"]:
                    entry_data["Organism"] = content.rstrip(".")
                else:
                    entry_data["Organism"] += " " + content.rstrip(".")

            # Check Metazoa
            elif line.startswith("OC"):
                in_sequence_block = False
                if "Metazoa" in content:
                    is_metazoa = True

            # Taxonomy ID
            elif line.startswith("OX"):
                in_sequence_block = False
                match = re.search(r"NCBI_TaxID=(\d+)", content)
                if match:
                    entry_data["Organism (ID)"] = match.group(1)

            # Keywords
            elif line.startswith("KW"):
                in_sequence_block = False
                keywords_text = content.rstrip(".")
                keywords = [kw.strip() for kw in keywords_text.split(";")]
                if "Toxin" in keywords:
                    has_kw0800_toxin = True

            # Protein names and fragments
            elif line.startswith("DE"):
                in_sequence_block = False
                if "Flags:" in content:
                    processing_de_flags = True
                    if "Fragment" in content:
                        entry_data["Fragment"] = "fragment"
                elif processing_de_flags and line.startswith("DE   "):
                    if "Fragment" in content:
                        entry_data["Fragment"] = "fragment"
                    if not content.endswith(";"):
                        processing_de_flags = False
                elif any(
                    x in content
                    for x in ["RecName: Full=", "AltName: Full=", "SubName: Full="]
                ):
                    processing_de_flags = False
                    name = re.search(
                        r"(?:RecName|AltName|SubName): Full=([^;]+)", content
                    )
                    if name:
                        if entry_data["Protein names"]:
                            entry_data["Protein names"] += "; " + name.group(1)
                        else:
                            entry_data["Protein names"] = name.group(1)
                # Handle old (pre-2014) simple DE format: "DE   Bucain."
                # Also handles "[Contains: ...]" and "[Includes: ...]" annotations
                elif not entry_data["Protein names"] and not any(
                    x in content
                    for x in [
                        "RecName:",
                        "AltName:",
                        "SubName:",
                        "Flags:",
                    ]
                ):
                    name = content.rstrip(".")
                    # Extract name before [Contains: or [Includes: if present
                    if "[Contains:" in name:
                        name = name.split("[Contains:")[0].strip()
                    if "[Includes:" in name:
                        name = name.split("[Includes:")[0].strip()
                    # Handle fragment notation in old format
                    if "(Fragment)" in name:
                        entry_data["Fragment"] = "fragment"
                        name = name.replace("(Fragment)", "").strip()
                    name = name.rstrip(".")
                    if name:
                        entry_data["Protein names"] = name

            # Comments
            elif line.startswith("CC"):
                in_sequence_block = False
                if "-!- CAUTION:" in content:
                    current_section = None
                elif "-!- FUNCTION:" in content:
                    current_section = "Function [CC]"
                    entry_data[current_section] = content.replace(
                        "-!- FUNCTION:", ""
                    ).strip()
                elif "-!- TISSUE SPECIFICITY:" in content:
                    current_section = "Tissue specificity"
                    entry_data[current_section] = content.replace(
                        "-!- TISSUE SPECIFICITY:", ""
                    ).strip()
                elif "-!- TOXIC DOSE:" in content:
                    current_section = "Toxic dose"
                    entry_data[current_section] = content.replace(
                        "-!- TOXIC DOSE:", ""
                    ).strip()
                elif "-!- PTM:" in content:
                    current_section = "Post-translational modification"
                    entry_data[current_section] = content.replace(
                        "-!- PTM:", ""
                    ).strip()
                elif "-!- MASS SPECTROMETRY:" in content:
                    current_section = "Mass spectrometry"
                    entry_data[current_section] = content.replace(
                        "-!- MASS SPECTROMETRY:", ""
                    ).strip()
                elif "-!- SIMILARITY:" in content:
                    current_section = "Protein families"
                    entry_data[current_section] = content.replace(
                        "-!- SIMILARITY:", ""
                    ).strip()
                elif content.startswith("-!-"):
                    current_section = None
                elif line.startswith("CC       ") and current_section:
                    entry_data[current_section] += " " + content

            # Feature table
            elif line.startswith("FT"):
                in_sequence_block = False
                feature_key = content.split()[0] if content else ""

                if feature_key == "SIGNAL":
                    entry_data["Signal peptide"] = "Yes"
                    range_match = re.search(r"(\d+)\.\.(\d+)", content)
                    if range_match:
                        entry_data["Signal peptide (range)"] = (
                            f"{range_match.group(1)}-{range_match.group(2)}"
                        )
                    else:
                        parts = content.split()
                        if (
                            len(parts) >= 3
                            and parts[1].isdigit()
                            and parts[2].isdigit()
                        ):
                            entry_data["Signal peptide (range)"] = (
                                f"{parts[1]}-{parts[2]}"
                            )

                elif feature_key == "PROPEP":
                    entry_data["Propeptide"] = "Yes"

                # PTM features
                elif feature_key in ptm_features:
                    position = self.parse_position(
                        content.split(None, 1)[1]
                        if len(content.split(None, 1)) > 1
                        else ""
                    )
                    note, evidence = self.extract_ptm_details(entry_lines, idx)

                    ptm_entry = {
                        "position": position,
                        "note": note,
                        "evidence": evidence,
                        "raw": content,
                    }
                    ptm_features[feature_key].append(ptm_entry)

            # Sequence
            elif line.startswith("SQ"):
                in_sequence_block = True
                entry_data["Sequence"] = ""
                match = re.search(r"(\d+) MW", content)
                if match:
                    entry_data["Mass"] = match.group(1)
                if "Fragment" in content and not entry_data["Fragment"]:
                    entry_data["Fragment"] = "fragment"

            # Database cross-references
            elif line.startswith("DR"):
                in_sequence_block = False
                if "InterPro" in content:
                    ipr_id = content.split(";")[1].strip()
                    entry_data["InterPro"] = (
                        (entry_data["InterPro"] + "; " + ipr_id)
                        if entry_data["InterPro"]
                        else ipr_id
                    )
                elif "Pfam" in content:
                    pfam_id = content.split(";")[1].strip()
                    entry_data["Pfam"] = (
                        (entry_data["Pfam"] + "; " + pfam_id)
                        if entry_data["Pfam"]
                        else pfam_id
                    )
                elif "KEGG" in content:
                    kegg_id = content.split(";")[1].strip()
                    entry_data["KEGG"] = (
                        (entry_data["KEGG"] + "; " + kegg_id)
                        if entry_data["KEGG"]
                        else kegg_id
                    )
                elif "GO;" in content:
                    go_match = re.search(r"GO:(\d+); ([CPF]):(.+?);", content)
                    if go_match:
                        go_id = "GO:" + go_match.group(1)
                        go_type = go_match.group(2)
                        go_desc = go_match.group(3).strip()

                        entry_data["Gene Ontology (GO)"] = (
                            (entry_data["Gene Ontology (GO)"] + "; " + go_id)
                            if entry_data["Gene Ontology (GO)"]
                            else go_id
                        )

                        if go_type == "P":
                            go_terms["P"].append(f"{go_id} ({go_desc})")
                        elif go_type == "C":
                            go_terms["C"].append(f"{go_id} ({go_desc})")
                        elif go_type == "F":
                            go_terms["F"].append(f"{go_id} ({go_desc})")

            # Protein existence
            elif line.startswith("PE"):
                in_sequence_block = False
                entry_data["Protein existence"] = content

            # Sequence data
            elif in_sequence_block and line.startswith("  "):
                sequence_part = line.strip().replace(" ", "")
                entry_data["Sequence"] += sequence_part

        # Process GO terms
        if go_terms["P"]:
            entry_data["Gene Ontology (biological process)"] = "; ".join(go_terms["P"])
        if go_terms["C"]:
            entry_data["Gene Ontology (cellular component)"] = "; ".join(go_terms["C"])
        if go_terms["F"]:
            entry_data["Gene Ontology (molecular function)"] = "; ".join(go_terms["F"])

        # Process PTM data
        self._process_ptm_data(entry_data, ptm_features)

        # Clean protein families
        self._clean_protein_families(entry_data)

        # Check criteria
        has_venom_tissue = "venom" in entry_data["Tissue specificity"].lower()
        meets_criteria = is_metazoa and (has_venom_tissue or has_kw0800_toxin)

        # Set ToxProt definition
        if has_venom_tissue and has_kw0800_toxin:
            entry_data["ToxProt definition"] = "both"
        elif has_venom_tissue:
            entry_data["ToxProt definition"] = "venom_tissue"
        elif has_kw0800_toxin:
            entry_data["ToxProt definition"] = "kw_toxin"

        return entry_data, meets_criteria

    def _process_ptm_data(self, entry_data, ptm_features):
        """Process and format PTM data"""
        ptm_summary = []

        # Disulfide bonds
        if ptm_features["DISULFID"]:
            detailed_list = []
            for ptm in ptm_features["DISULFID"]:
                detail = ptm["position"]
                if ptm["note"] and ptm["evidence"]:
                    detail += f" ({ptm['note']}) {{{ptm['evidence']}}}"
                elif ptm["note"]:
                    detail += f" ({ptm['note']})"
                elif ptm["evidence"]:
                    detail += f" {{{ptm['evidence']}}}"
                detailed_list.append(detail)
            entry_data["Disulfide bond"] = "; ".join(detailed_list)
            ptm_summary.append(f"DISULFID:{len(ptm_features['DISULFID'])}")

        # Glycosylation
        if ptm_features["CARBOHYD"]:
            detailed_list = []
            for ptm in ptm_features["CARBOHYD"]:
                detail = ptm["position"]
                if ptm["note"] and ptm["evidence"]:
                    detail += f" ({ptm['note']}) {{{ptm['evidence']}}}"
                elif ptm["note"]:
                    detail += f" ({ptm['note']})"
                elif ptm["evidence"]:
                    detail += f" {{{ptm['evidence']}}}"
                detailed_list.append(detail)
            entry_data["Glycosylation"] = "; ".join(detailed_list)
            ptm_summary.append(f"CARBOHYD:{len(ptm_features['CARBOHYD'])}")

        # Modified residues
        if ptm_features["MOD_RES"]:
            detailed_list = []
            for ptm in ptm_features["MOD_RES"]:
                detail = ptm["position"]
                if ptm["note"] and ptm["evidence"]:
                    detail += f" ({ptm['note']}) {{{ptm['evidence']}}}"
                elif ptm["note"]:
                    detail += f" ({ptm['note']})"
                elif ptm["evidence"]:
                    detail += f" {{{ptm['evidence']}}}"
                detailed_list.append(detail)
            entry_data["Modified residue"] = "; ".join(detailed_list)
            ptm_summary.append(f"MOD_RES:{len(ptm_features['MOD_RES'])}")

        # Cross-links
        if ptm_features["CROSSLNK"]:
            detailed_list = []
            for ptm in ptm_features["CROSSLNK"]:
                detail = ptm["position"]
                if ptm["note"] and ptm["evidence"]:
                    detail += f" ({ptm['note']}) {{{ptm['evidence']}}}"
                elif ptm["note"]:
                    detail += f" ({ptm['note']})"
                elif ptm["evidence"]:
                    detail += f" {{{ptm['evidence']}}}"
                detailed_list.append(detail)
            entry_data["Cross-link"] = "; ".join(detailed_list)
            ptm_summary.append(f"CROSSLNK:{len(ptm_features['CROSSLNK'])}")

        # Lipidation
        if ptm_features["LIPID"]:
            detailed_list = []
            for ptm in ptm_features["LIPID"]:
                detail = ptm["position"]
                if ptm["note"] and ptm["evidence"]:
                    detail += f" ({ptm['note']}) {{{ptm['evidence']}}}"
                elif ptm["note"]:
                    detail += f" ({ptm['note']})"
                elif ptm["evidence"]:
                    detail += f" {{{ptm['evidence']}}}"
                detailed_list.append(detail)
            entry_data["Lipidation"] = "; ".join(detailed_list)
            ptm_summary.append(f"LIPID:{len(ptm_features['LIPID'])}")

        # PTM summary
        if ptm_summary:
            entry_data["PTM Summary"] = "; ".join(ptm_summary)

    def _clean_protein_families(self, entry_data):
        """Clean up protein families field"""
        if entry_data["Protein families"]:
            family_text = re.sub(
                r"\s*\{ECO:[^\}]+\}\.?", "", entry_data["Protein families"]
            ).strip()
            match_belongs = re.search(
                r"belongs to the\s+(.+)", family_text, re.IGNORECASE
            )
            if match_belongs:
                processed_family_text = match_belongs.group(1)
            else:
                processed_family_text = re.sub(
                    r"^Belongs to the\s+", "", family_text, flags=re.IGNORECASE
                )

            if processed_family_text.endswith("."):
                processed_family_text = processed_family_text[:-1]

            stripped_text = processed_family_text.strip()
            if stripped_text:
                entry_data["Protein families"] = (
                    stripped_text[0].upper() + stripped_text[1:]
                )
            else:
                entry_data["Protein families"] = stripped_text


def process_swissprot_file(input_file, output_file, ptm_vocab=None):
    """Process SwissProt file and extract matching entries.

    Args:
        input_file: Path to the input .dat file.
        output_file: Path to the output .tsv file.
        ptm_vocab: Optional PTM vocabulary dictionary.

    Returns:
        True if processing was successful, False otherwise.
        The number of entries can be retrieved by checking the output file.
    """
    parser = SwissProtParser(ptm_vocab)
    entries_found = []

    try:
        if not Path(input_file).exists():
            logger.error(f"Input file {input_file} not found")
            return False

        # Ensure output directory exists
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)

        # Count total entries
        logger.debug("Counting total entries...")
        total_entries = 0
        with open(input_file, encoding="utf-8") as f:
            for line in f:
                if line.startswith("//"):
                    total_entries += 1
        logger.debug(f"Found {total_entries:,} entries")

        # Process entries (show progress bar only if logging at DEBUG level)
        show_progress = logger.isEnabledFor(logging.DEBUG)
        with open(input_file, encoding="utf-8") as f:
            current_entry = []
            pbar = tqdm(total=total_entries, desc="Processing entries", disable=not show_progress)

            for line in f:
                if line.startswith("//"):
                    pbar.update(1)
                    if current_entry:
                        entry_data, meets_criteria = parser.parse_entry(current_entry)
                        if meets_criteria:
                            entries_found.append(entry_data)
                        current_entry = []
                else:
                    current_entry.append(line)

            pbar.close()

        # Write TSV output
        logger.debug(f"Writing {len(entries_found):,} entries to {output_file}")
        with open(output_file, "w", encoding="utf-8") as out:
            out.write("\t".join(parser.HEADERS) + "\n")
            for entry in entries_found:
                row = [entry.get(field, "") for field in parser.HEADERS]
                out.write("\t".join(row) + "\n")

        logger.debug(f"Processing complete: {len(entries_found):,} matching entries")
        return True

    except Exception as e:
        logger.error(f"Error processing file: {e}", exc_info=True)
        return False


def get_output_filename(input_path: Path) -> str:
    """Generate output filename from input filename (e.g., 2005_sprot.dat -> toxprot_2005.tsv)"""
    stem = input_path.stem  # e.g., "2005_sprot"
    # Extract year from filename
    match = re.match(r"(\d{4})", stem)
    if match:
        year = match.group(1)
        return f"toxprot_{year}.tsv"
    # Fallback: use input stem
    return f"toxprot_{stem}.tsv"


def main():
    # Configure logging for standalone execution
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s - %(message)s",
        datefmt="%H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
        ],
    )

    parser = argparse.ArgumentParser(
        description="SwissProt parser for Tox-Prot (animal toxin annotation project)",
        epilog="""Examples:
  # Process single file
  python parse_sprot_dat.py data/raw/uniprot_releases/2005_sprot.dat

  # Process multiple files
  python parse_sprot_dat.py data/raw/uniprot_releases/2005_sprot.dat data/raw/uniprot_releases/2010_sprot.dat

  # Process all files in directory
  python parse_sprot_dat.py --input-dir data/raw/uniprot_releases

  # Process specific years
  python parse_sprot_dat.py --years 2005 2010 2015 2020 2025""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input_files",
        type=str,
        nargs="*",
        help="Path(s) to SwissProt data file(s)",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Directory containing SwissProt .dat files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/interim/toxprot_parsed"),
        help="Directory for output files (default: data/interim/toxprot_parsed)",
    )
    parser.add_argument(
        "--years",
        type=int,
        nargs="+",
        help="Years to process (looks for {year}_sprot.dat in input-dir)",
    )
    parser.add_argument(
        "--data-dir",
        type=str,
        default="data/raw",
        help="Directory for PTM vocabulary files (default: data/raw)",
    )
    parser.add_argument(
        "--skip-ptmlist",
        action="store_true",
        help="Skip downloading/loading ptmlist.txt",
    )
    parser.add_argument(
        "--delete-input",
        action="store_true",
        help="Delete input .dat files after successful processing (to save disk space)",
    )

    args = parser.parse_args()

    # Determine input files
    input_files = []

    if args.years:
        # Process specific years from input-dir
        input_dir = args.input_dir or Path("data/raw/uniprot_releases")
        for year in args.years:
            input_path = input_dir / f"{year}_sprot.dat"
            if input_path.exists():
                input_files.append(input_path)
            else:
                logger.warning(f"File not found: {input_path}")
    elif args.input_dir:
        # Process all .dat files in directory
        input_files = sorted(args.input_dir.glob("*_sprot.dat"))
    elif args.input_files:
        # Process specified files
        input_files = [Path(f) for f in args.input_files]
    else:
        parser.error("Specify input files, --input-dir, or --years")

    if not input_files:
        logger.error("No input files found")
        return

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("=" * 60)
    logger.info("SwissProt Parser for Tox-Prot")
    logger.info("=" * 60)
    logger.info(f"Input files: {len(input_files)}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(
        "Query: (taxonomy_id:33208) AND (cc_tissue_specificity:venom OR keyword:KW-0800)"
    )
    logger.info("=" * 60)

    # Initialize PTM vocabulary (once for all files)
    ptm_vocab = None
    if not args.skip_ptmlist:
        vocab_manager = PTMVocabulary(args.data_dir)
        vocab_manager.ensure_ptmlist_available()
        vocab_manager.load_ptmlist()
        ptm_vocab = vocab_manager.ptm_vocab

    # Process each file
    successful = 0
    for input_path in input_files:
        if not input_path.exists():
            logger.warning(f"Skipping non-existent file: {input_path}")
            continue

        output_filename = get_output_filename(input_path)
        output_path = args.output_dir / output_filename

        logger.info("-" * 60)
        logger.info(f"Processing: {input_path.name}")
        logger.info(f"Output: {output_path}")

        if process_swissprot_file(input_path, output_path, ptm_vocab):
            successful += 1
            if args.delete_input:
                logger.info(f"Deleting input file: {input_path}")
                input_path.unlink()

    logger.info("=" * 60)
    logger.info(f"Processing complete! {successful}/{len(input_files)} file(s) processed successfully")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
