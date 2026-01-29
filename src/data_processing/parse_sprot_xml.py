#!/usr/bin/env python3
"""
UniProtKB/SwissProt XML parser for Tox-Prot (animal toxin annotation project).

Streaming XML parser that extracts structured PTM data. Raw PTM data is stored
in a consolidated PTM_Features JSON column; summaries and keyword resolution
are computed during the clean_data.py processing step.

Features:
- Streaming iterparse for memory-efficient processing of large XML files
- Namespace auto-detection (http:// for ≤2024, https:// for 2025+)
- Consolidated PTM_Features JSON column containing raw PTM data
- PTM Keywords extraction from <keyword> elements
"""

import json
import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from tqdm import tqdm

logger = logging.getLogger(__name__)


def load_ptm_keyword_ids(path: str = "data/raw/ptm_keywords.tsv") -> set:
    """Load set of PTM keyword IDs from ptm_keywords.tsv.

    Args:
        path: Path to the ptm_keywords.tsv file.

    Returns:
        Set of PTM keyword IDs (e.g., {"KW-0379", "KW-0027", ...}).
    """
    try:
        df = pd.read_csv(path, sep="\t")
        return set(df["Keyword ID"].tolist())
    except FileNotFoundError:
        logger.warning(f"PTM keywords file not found: {path}")
        return set()


# Load PTM keyword IDs at module level
PTM_KEYWORD_IDS = load_ptm_keyword_ids()

# XML feature type → readable PTM name for JSON output
FEATURE_TYPE_MAP = {
    "modified residue": "Modified residue",
    "glycosylation site": "Glycosylation",
    "disulfide bond": "Disulfide bond",
    "cross-link": "Cross-link",
    "lipid moiety-binding region": "Lipidation",
}


class SwissProtXMLParser:
    """SwissProt XML parser for extracting protein information and features."""

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
        "PTM_Features",
        "PTM Keywords",
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

    # PTM readable names for JSON structure (order matters for consistent output)
    PTM_NAMES = ["Disulfide bond", "Glycosylation", "Modified residue", "Cross-link", "Lipidation"]

    def __init__(self, ns=""):
        """Initialize parser with XML namespace.

        Args:
            ns: XML namespace string (e.g., "{http://uniprot.org/uniprot}").
        """
        self.ns = ns

    def _tag(self, local_name):
        """Return a fully qualified tag name."""
        return f"{self.ns}{local_name}"

    def _find_text(self, elem, path):
        """Find element and return its text, or empty string."""
        parts = path.split("/")
        current = elem
        for part in parts:
            current = current.find(self._tag(part))
            if current is None:
                return ""
        return (current.text or "").strip()

    def _findall(self, elem, path):
        """Find all elements matching a tag path."""
        parts = path.split("/")
        if len(parts) == 1:
            return elem.findall(self._tag(parts[0]))
        # Navigate to parent, then findall on last part
        current = elem
        for part in parts[:-1]:
            current = current.find(self._tag(part))
            if current is None:
                return []
        return current.findall(self._tag(parts[-1]))

    def parse_entry(self, entry_elem):
        """Parse an XML <entry> element and extract all protein information.

        Args:
            entry_elem: An xml.etree.ElementTree Element for an <entry>.

        Returns:
            Tuple of (entry_data dict, meets_criteria bool).
        """
        entry_data = {field: "" for field in self.HEADERS}

        # --- Accession ---
        accessions = entry_elem.findall(self._tag("accession"))
        if accessions:
            entry_data["Entry"] = (accessions[0].text or "").strip()

        # --- Entry Name ---
        name_elem = entry_elem.find(self._tag("name"))
        if name_elem is not None:
            entry_data["Entry Name"] = (name_elem.text or "").strip()

        # --- Organism ---
        organism_elem = entry_elem.find(self._tag("organism"))
        is_metazoa = False
        if organism_elem is not None:
            sci_name = organism_elem.find(self._tag("name") + "[@type='scientific']")
            # Fallback for old format (≤2005): <name type="full">
            if sci_name is None:
                sci_name = organism_elem.find(self._tag("name") + "[@type='full']")
            if sci_name is not None:
                org_str = (sci_name.text or "").strip()
                # Include common name to match .dat format: "Species (Common name)"
                common_name = organism_elem.find(self._tag("name") + "[@type='common']")
                if common_name is not None and common_name.text:
                    org_str += f" ({common_name.text.strip()})"
                entry_data["Organism"] = org_str

            # Taxonomy ID
            tax_ref = organism_elem.find(self._tag("dbReference") + "[@type='NCBI Taxonomy']")
            if tax_ref is not None:
                entry_data["Organism (ID)"] = tax_ref.get("id", "")

            # Check lineage for Metazoa
            lineage = organism_elem.find(self._tag("lineage"))
            if lineage is not None:
                for taxon in lineage.findall(self._tag("taxon")):
                    if (taxon.text or "").strip() == "Metazoa":
                        is_metazoa = True
                        break

        # --- Protein names ---
        protein_elem = entry_elem.find(self._tag("protein"))
        if protein_elem is not None:
            names = []
            # Modern format: <recommendedName><fullName>...</fullName></recommendedName>
            for name_group_tag in ["recommendedName", "alternativeName", "submittedName"]:
                for name_group in protein_elem.findall(self._tag(name_group_tag)):
                    full_name = name_group.find(self._tag("fullName"))
                    if full_name is not None and full_name.text:
                        names.append(full_name.text.strip())
            # Old format (≤2005): <protein><name>Text</name></protein>
            if not names:
                for name_elem in protein_elem.findall(self._tag("name")):
                    if name_elem.text:
                        names.append(name_elem.text.strip())
            entry_data["Protein names"] = "; ".join(names)

            # Fragment: <protein type="fragment">
            if protein_elem.get("type") == "fragment":
                entry_data["Fragment"] = "fragment"

        # --- Sequence / Fragment ---
        seq_elem = entry_elem.find(self._tag("sequence"))
        if seq_elem is not None:
            frag = seq_elem.get("fragment")
            if frag:
                entry_data["Fragment"] = "fragment"
            entry_data["Length"] = seq_elem.get("length", "")
            entry_data["Mass"] = seq_elem.get("mass", "")
            entry_data["Sequence"] = (seq_elem.text or "").replace("\n", "").replace(" ", "")

        # --- Keywords ---
        has_kw0800_toxin = False
        ptm_keywords = []
        for kw_elem in entry_elem.findall(self._tag("keyword")):
            kw_id = kw_elem.get("id", "")
            if kw_id == "KW-0800":
                has_kw0800_toxin = True
            if kw_id in PTM_KEYWORD_IDS:
                kw_text = (kw_elem.text or "").strip()
                if kw_text:
                    ptm_keywords.append(kw_text)
        entry_data["PTM Keywords"] = "; ".join(sorted(ptm_keywords))

        # --- Comments ---
        comment_type_map = {
            "function": "Function [CC]",
            "tissue specificity": "Tissue specificity",
            "toxic dose": "Toxic dose",
            "similarity": "Protein families",
        }
        for comment in entry_elem.findall(self._tag("comment")):
            ctype = comment.get("type", "")

            # Mass spectrometry uses attributes, not <text>
            if ctype == "mass spectrometry":
                ms_parts = []
                mass = comment.get("mass", "")
                method = comment.get("method", "")
                if mass:
                    ms_parts.append(f"MW={mass}")
                if method:
                    ms_parts.append(f"METHOD={method}")
                # Extract range from <location>
                loc = comment.find(self._tag("location"))
                if loc is not None:
                    begin = loc.find(self._tag("begin"))
                    end = loc.find(self._tag("end"))
                    if begin is not None and end is not None:
                        b = begin.get("position", "")
                        e = end.get("position", "")
                        if b and e:
                            ms_parts.append(f"RANGE={b}-{e}")
                # Extract note
                note_elem = comment.find(self._tag("note"))
                if note_elem is not None and note_elem.text:
                    ms_parts.append(f"NOTE={note_elem.text.strip()}")
                # Also check for <text> (newer format)
                text_elem = comment.find(self._tag("text"))
                if text_elem is not None and text_elem.text:
                    ms_val = text_elem.text.strip()
                elif ms_parts:
                    ms_val = "; ".join(ms_parts) + "."
                else:
                    continue
                if entry_data["Mass spectrometry"]:
                    entry_data["Mass spectrometry"] += " " + ms_val
                else:
                    entry_data["Mass spectrometry"] = ms_val
                continue

            header_key = comment_type_map.get(ctype)
            if header_key:
                text_elem = comment.find(self._tag("text"))
                if text_elem is not None and text_elem.text:
                    text_val = text_elem.text.strip()
                    # Append status attribute if present (e.g., "By similarity", "Probable")
                    status = comment.get("status") or text_elem.get("status")
                    if status:
                        text_val += f" ({status.capitalize()})"
                    # Ensure trailing period to match .dat format
                    if text_val and not text_val.endswith("."):
                        text_val += "."
                    if entry_data[header_key]:
                        entry_data[header_key] += " " + text_val
                    else:
                        entry_data[header_key] = text_val

        # Clean protein families
        self._clean_protein_families(entry_data)

        # --- Build evidence map ---
        evidence_map = {}
        for ev in entry_elem.findall(self._tag("evidence")):
            key = ev.get("key", "")
            etype = ev.get("type", "")
            if key and etype:
                evidence_map[key] = etype

        # --- Features ---
        # Use readable names as keys for PTM features
        ptm_features = {name: [] for name in self.PTM_NAMES}

        for feature in entry_elem.findall(self._tag("feature")):
            ftype = feature.get("type", "")

            # Signal peptide
            if ftype == "signal peptide":
                entry_data["Signal peptide"] = "Yes"
                loc = feature.find(self._tag("location"))
                if loc is not None:
                    begin = loc.find(self._tag("begin"))
                    end = loc.find(self._tag("end"))
                    if begin is not None and end is not None:
                        b = begin.get("position", "?")
                        e = end.get("position", "?")
                        entry_data["Signal peptide (range)"] = f"{b}-{e}"
                continue

            # Propeptide
            if ftype == "propeptide":
                entry_data["Propeptide"] = "Yes"
                continue

            # PTM features
            ptm_name = FEATURE_TYPE_MAP.get(ftype)
            if ptm_name:
                description = feature.get("description", "")
                evidence_str = self._resolve_evidence(feature.get("evidence", ""), evidence_map)
                position = self._extract_position(feature)

                ptm_features[ptm_name].append(
                    {
                        "position": position,
                        "note": description,
                        "evidence": evidence_str,
                    }
                )

        # Process PTM data
        self._process_ptm_data(entry_data, ptm_features)

        # --- Database cross-references ---
        go_terms = {"P": [], "C": [], "F": []}
        for dbref in entry_elem.findall(self._tag("dbReference")):
            dbtype = dbref.get("type", "")
            dbid = dbref.get("id", "")

            if dbtype == "InterPro":
                entry_data["InterPro"] = (
                    (entry_data["InterPro"] + "; " + dbid) if entry_data["InterPro"] else dbid
                )
            elif dbtype == "Pfam":
                entry_data["Pfam"] = (
                    (entry_data["Pfam"] + "; " + dbid) if entry_data["Pfam"] else dbid
                )
            elif dbtype == "KEGG":
                entry_data["KEGG"] = (
                    (entry_data["KEGG"] + "; " + dbid) if entry_data["KEGG"] else dbid
                )
            elif dbtype == "GO":
                entry_data["Gene Ontology (GO)"] = (
                    (entry_data["Gene Ontology (GO)"] + "; " + dbid)
                    if entry_data["Gene Ontology (GO)"]
                    else dbid
                )
                # Parse GO term category
                for prop in dbref.findall(self._tag("property")):
                    if prop.get("type") == "term":
                        term_val = prop.get("value", "")
                        if term_val and ":" in term_val:
                            go_type = term_val[0]  # P, C, or F
                            go_desc = term_val[2:]  # after "X:"
                            if go_type in go_terms:
                                go_terms[go_type].append(f"{dbid} ({go_desc})")

        # Process GO terms
        if go_terms["P"]:
            entry_data["Gene Ontology (biological process)"] = "; ".join(go_terms["P"])
        if go_terms["C"]:
            entry_data["Gene Ontology (cellular component)"] = "; ".join(go_terms["C"])
        if go_terms["F"]:
            entry_data["Gene Ontology (molecular function)"] = "; ".join(go_terms["F"])

        # --- Protein existence ---
        pe_elem = entry_elem.find(self._tag("proteinExistence"))
        if pe_elem is not None:
            pe_type = pe_elem.get("type", "")
            # Map XML type to dat-style format
            pe_map = {
                "evidence at protein level": "1: Evidence at protein level",
                "evidence at transcript level": "2: Evidence at transcript level",
                "inferred from homology": "3: Inferred from homology",
                "predicted": "4: Predicted",
                "uncertain": "5: Uncertain",
            }
            entry_data["Protein existence"] = pe_map.get(pe_type, pe_type)

        # --- Criteria check ---
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

    def _resolve_evidence(self, evidence_attr, evidence_map):
        """Resolve evidence attribute keys to ECO codes.

        Args:
            evidence_attr: Space-separated evidence keys from XML (e.g., "1 2 3").
            evidence_map: Dict mapping key → ECO type string.

        Returns:
            Comma-separated ECO codes string, or empty string.
        """
        if not evidence_attr:
            return ""
        codes = []
        for key in evidence_attr.split():
            eco = evidence_map.get(key, "")
            if eco:
                codes.append(eco)
        return ", ".join(codes)

    def _extract_position(self, feature_elem):
        """Extract position string from a feature element's <location>."""
        loc = feature_elem.find(self._tag("location"))
        if loc is None:
            return ""

        # Single position
        pos = loc.find(self._tag("position"))
        if pos is not None:
            return pos.get("position", "")

        # Range
        begin = loc.find(self._tag("begin"))
        end = loc.find(self._tag("end"))
        if begin is not None and end is not None:
            b = begin.get("position", "?")
            e = end.get("position", "?")
            if b == e:
                return b
            return f"{b}-{e}"

        return ""

    def _process_ptm_data(self, entry_data, ptm_features):
        """Process and format PTM data into a single PTM_Features JSON field.

        Stores raw PTM data only. Summaries and keyword resolution are computed
        during the clean_data.py processing step.

        The JSON structure:
        {
            "Disulfide bond": [{"pos": "42-146", "note": "Interchain", "evidence": "ECO:0000269"}, ...],
            "Modified residue": [{"pos": "49", "note": "4-hydroxyproline", "evidence": "ECO:0000305"}, ...],
            "Glycosylation": [...],
            "Cross-link": [...],
            "Lipidation": [...]
        }
        """
        result = {}

        for ptm_name in self.PTM_NAMES:
            features = ptm_features.get(ptm_name, [])
            if not features:
                continue

            ptm_list = [
                {
                    "pos": ptm["position"],
                    "note": ptm["note"],
                    "evidence": ptm["evidence"],
                }
                for ptm in features
            ]
            result[ptm_name] = ptm_list

        entry_data["PTM_Features"] = json.dumps(result, separators=(",", ":")) if result else ""

    def _clean_protein_families(self, entry_data):
        """Clean up protein families field (same logic as .dat parser)."""
        if entry_data["Protein families"]:
            family_text = re.sub(
                r"\s*\{ECO:[^\}]+\}\.?", "", entry_data["Protein families"]
            ).strip()
            match_belongs = re.search(r"belongs to the\s+(.+)", family_text, re.IGNORECASE)
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
                entry_data["Protein families"] = stripped_text[0].upper() + stripped_text[1:]
            else:
                entry_data["Protein families"] = stripped_text


def detect_namespace(input_file):
    """Detect the XML namespace from the root element.

    UniProt uses 'http://uniprot.org/uniprot' for ≤2024 releases and
    'https://uniprot.org/uniprot' for 2025+ releases.

    Args:
        input_file: Path to the XML file.

    Returns:
        Namespace string in "{uri}" format, or "" if none found.
    """
    for _event, elem in ET.iterparse(str(input_file), events=("start",)):
        tag = elem.tag
        if tag.startswith("{"):
            ns_uri = tag[1 : tag.index("}")]
            return f"{{{ns_uri}}}"
        return ""
    return ""


def process_xml_file(input_file, output_file):
    """Process a UniProt XML file and extract matching entries.

    Streaming parser that processes entries one at a time to minimize memory.
    Extracts raw PTM data; summaries and keyword resolution are done in clean_data.py.

    Args:
        input_file: Path to the input .xml file.
        output_file: Path to the output .tsv file.

    Returns:
        True if processing was successful, False otherwise.
    """
    try:
        input_file = Path(input_file)
        output_file = Path(output_file)

        if not input_file.exists():
            logger.error(f"Input file {input_file} not found")
            return False

        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Detect namespace
        logger.debug("Detecting XML namespace...")
        ns = detect_namespace(input_file)
        logger.debug(f"Namespace: {ns or '(none)'}")

        parser = SwissProtXMLParser(ns=ns)
        entries_found = []

        # Streaming parse (single pass, no pre-scan for count)
        show_progress = logger.isEnabledFor(logging.DEBUG)
        entry_tag = f"{ns}entry"
        root = None
        processed = 0

        context = ET.iterparse(str(input_file), events=("start", "end"))
        pbar = tqdm(desc="Processing entries", unit=" entries", disable=not show_progress)

        for event, elem in context:
            if event == "start" and root is None:
                root = elem
            elif event == "end" and elem.tag == entry_tag:
                processed += 1
                pbar.update(1)
                entry_data, meets_criteria = parser.parse_entry(elem)
                if meets_criteria:
                    entries_found.append(entry_data)
                root.clear()

        pbar.close()
        logger.debug(f"Processed {processed:,} total entries")

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
        logger.error(f"Error processing XML file: {e}", exc_info=True)
        return False
