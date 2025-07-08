#!/usr/bin/env python3
"""
SwissProt Parser for Venom/Toxin Proteins

This script parses a UniProtKB/SwissProt data file to extract entries that match
the query: (taxonomy_id:33208) AND ((cc_tissue_specificity:venom) OR (keyword:KW-0800))
https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/release-2017_11/

It extracts specific fields and outputs them in TSV format.
"""

import re
import argparse
from pathlib import Path
from tqdm import tqdm

# Define header fields
headers = [
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
    "Sequence",
]


def parse_entry(entry_lines):
    """Parse a SwissProt entry and extract the required fields."""
    entry_data = {field: "" for field in headers}

    # Initialize variables
    current_section = None
    is_metazoa = False
    has_toxin_keyword = False
    go_terms = {"P": [], "C": [], "F": []}
    processing_de_flags = False
    in_sequence_block = False  # Flag to indicate if we are processing sequence lines

    # Process each line in the entry
    for line in entry_lines:
        content = line[5:].strip()

        # Extract the entry ID
        if line.startswith("ID"):
            in_sequence_block = False  # End sequence block
            match = re.search(r"^(\S+)", content)
            if match:
                entry_data["Entry Name"] = match.group(1)
                # Extract length
                length_match = re.search(r"(\d+) AA", content)
                if length_match:
                    entry_data["Length"] = length_match.group(1)

        # Extract accession
        elif line.startswith("AC"):
            in_sequence_block = False  # End sequence block
            if not entry_data["Entry"]:  # Only take the first accession
                primary_acc = content.split(";")[0].strip()
                entry_data["Entry"] = primary_acc

        # Extract organism information
        elif line.startswith("OS"):
            in_sequence_block = False  # End sequence block
            if not entry_data["Organism"]:
                entry_data["Organism"] = content.rstrip(".")
            else:
                entry_data["Organism"] += " " + content.rstrip(".")

        # Check if organism is Metazoa
        elif line.startswith("OC"):
            in_sequence_block = False  # End sequence block
            if "Metazoa" in content:
                is_metazoa = True

        # Extract taxonomy ID
        elif line.startswith("OX"):
            in_sequence_block = False  # End sequence block
            match = re.search(r"NCBI_TaxID=(\d+)", content)
            if match:
                entry_data["Organism (ID)"] = match.group(1)

        # Extract protein names and check for fragment flags
        elif line.startswith("DE"):
            in_sequence_block = False  # End sequence block
            # Check for Fragment in DE Flags section
            if "Flags:" in content:
                processing_de_flags = True
                if "Fragment" in content:
                    entry_data["Fragment"] = "fragment"
            elif processing_de_flags and line.startswith("DE   "):
                # Continue processing DE Flags on subsequent lines
                if "Fragment" in content:
                    entry_data["Fragment"] = "fragment"
                # Stop processing if we hit a line that doesn't continue the flags
                if not content.endswith(";"):
                    processing_de_flags = False
            elif (
                "RecName: Full=" in content
                or "AltName: Full=" in content
                or "SubName: Full=" in content
            ):
                # No longer processing flags
                processing_de_flags = False
                name = re.search(r"(?:RecName|AltName|SubName): Full=([^;]+)", content)
                if name:
                    if entry_data["Protein names"]:
                        entry_data["Protein names"] += "; " + name.group(1)
                    else:
                        entry_data["Protein names"] = name.group(1)

        # Extract comments for function, tissue specificity, similarity, etc.
        elif line.startswith("CC"):
            in_sequence_block = False  # End sequence block
            # If we hit a CAUTION section, reset current_section to stop adding content
            if "-!- CAUTION:" in content:
                current_section = None
            # Process other sections
            elif "-!- FUNCTION:" in content:
                current_section = "Function [CC]"
                function_text = content.replace("-!- FUNCTION:", "").strip()
                entry_data[current_section] = function_text
            elif "-!- TISSUE SPECIFICITY:" in content:
                current_section = "Tissue specificity"
                tissue_text = content.replace("-!- TISSUE SPECIFICITY:", "").strip()
                entry_data[current_section] = tissue_text
            elif "-!- TOXIC DOSE:" in content:
                current_section = "Toxic dose"
                toxic_text = content.replace("-!- TOXIC DOSE:", "").strip()
                entry_data[current_section] = toxic_text
            elif "-!- PTM:" in content:
                current_section = "Post-translational modification"
                ptm_text = content.replace("-!- PTM:", "").strip()
                entry_data[current_section] = ptm_text
            elif "-!- MASS SPECTROMETRY:" in content:
                current_section = "Mass spectrometry"
                mass_spec_text = content.replace("-!- MASS SPECTROMETRY:", "").strip()
                entry_data[current_section] = mass_spec_text
            elif "-!- SIMILARITY:" in content:
                current_section = "Protein families"
                similarity_text = content.replace("-!- SIMILARITY:", "").strip()
                # Store the relatively raw text; detailed parsing is deferred.
                entry_data[current_section] = similarity_text
            # If it's a new CC -!- section header that we don't explicitly handle,
            # reset current_section to prevent appending its content to the previous section.
            elif content.startswith("-!-"):
                current_section = None
            elif line.startswith("CC       ") and current_section:
                # Continuation of a previous comment section
                entry_data[current_section] += " " + content

        # Extract feature information for signal peptide, propeptide, disulfide bonds
        elif line.startswith("FT"):
            in_sequence_block = False  # End sequence block
            feature_key = content.split()[0] if content else ""

            if feature_key == "SIGNAL":
                entry_data["Signal peptide"] = "Yes"
                # Try to parse new format first (e.g., "1..21")
                range_match_new = re.search(r"(\d+)\.\.(\d+)", content)
                if range_match_new:
                    entry_data["Signal peptide (range)"] = (
                        f"{range_match_new.group(1)}-{range_match_new.group(2)}"
                    )
                else:
                    # Try to parse old format (e.g., "1     21")
                    parts = content.split()
                    if len(parts) >= 3 and parts[1].isdigit() and parts[2].isdigit():
                        entry_data["Signal peptide (range)"] = f"{parts[1]}-{parts[2]}"
            elif feature_key == "PROPEP":
                entry_data["Propeptide"] = "Yes"
            elif feature_key == "DISULFID":
                if entry_data["Disulfide bond"]:
                    entry_data["Disulfide bond"] += "; " + content
                else:
                    entry_data["Disulfide bond"] = content
            elif feature_key == "CARBOHYD":
                if entry_data["Glycosylation"]:
                    entry_data["Glycosylation"] += "; " + content
                else:
                    entry_data["Glycosylation"] = content

        # Extract sequence information for mass
        elif line.startswith("SQ"):
            in_sequence_block = True  # Start of sequence block
            entry_data["Sequence"] = ""  # Initialize sequence string for this entry

            # Existing logic to parse Mass and Fragment from SQ line content
            match = re.search(r"(\d+) MW", content)
            if match:
                entry_data["Mass"] = match.group(1)

            # Also check for fragment in SQ line as a backup
            if "Fragment" in content and not entry_data["Fragment"]:
                entry_data["Fragment"] = "fragment"

        # Extract database cross-references
        elif line.startswith("DR"):
            in_sequence_block = False  # End sequence block
            if "InterPro" in content:
                ipr_id = content.split(";")[1].strip()
                if entry_data["InterPro"]:
                    entry_data["InterPro"] += "; " + ipr_id
                else:
                    entry_data["InterPro"] = ipr_id
            elif "Pfam" in content:
                pfam_id = content.split(";")[1].strip()
                if entry_data["Pfam"]:
                    entry_data["Pfam"] += "; " + pfam_id
                else:
                    entry_data["Pfam"] = pfam_id
            elif "KEGG" in content:
                kegg_id = content.split(";")[1].strip()
                if entry_data["KEGG"]:
                    entry_data["KEGG"] += "; " + kegg_id
                else:
                    entry_data["KEGG"] = kegg_id
            elif "GO;" in content:
                go_match = re.search(r"GO:(\d+); ([CPF]):(.+?);", content)
                if go_match:
                    go_id = "GO:" + go_match.group(1)
                    go_type = go_match.group(2)
                    go_desc = go_match.group(3).strip()

                    # Add to overall GO terms
                    if entry_data["Gene Ontology (GO)"]:
                        entry_data["Gene Ontology (GO)"] += "; " + go_id
                    else:
                        entry_data["Gene Ontology (GO)"] = go_id

                    # Add to specific GO category
                    if go_type == "P":
                        go_terms["P"].append(f"{go_id} ({go_desc})")
                    elif go_type == "C":
                        go_terms["C"].append(f"{go_id} ({go_desc})")
                    elif go_type == "F":
                        go_terms["F"].append(f"{go_id} ({go_desc})")

        # Extract keywords
        elif line.startswith("KW"):
            in_sequence_block = False  # End sequence block
            if "Toxin" in content:
                has_toxin_keyword = True

        # Extract protein existence
        elif line.startswith("PE"):
            in_sequence_block = False  # End sequence block
            entry_data["Protein existence"] = content

        # Handle sequence data lines if in sequence block
        # These lines start with spaces and follow an SQ line.
        elif in_sequence_block and line.startswith("  "):
            # For sequence data lines, process the raw line content directly
            sequence_part = line.strip().replace(" ", "")
            entry_data["Sequence"] += sequence_part

    # Process GO terms
    if go_terms["P"]:
        entry_data["Gene Ontology (biological process)"] = "; ".join(go_terms["P"])
    if go_terms["C"]:
        entry_data["Gene Ontology (cellular component)"] = "; ".join(go_terms["C"])
    if go_terms["F"]:
        entry_data["Gene Ontology (molecular function)"] = "; ".join(go_terms["F"])

    # Clean up protein families field - remove ECO codes and similar annotations
    if entry_data["Protein families"]:
        # Remove ECO codes and similar annotations
        family_text = re.sub(
            r"\s*\{ECO:[^\}]+\}\.?", "", entry_data["Protein families"]
        ).strip()

        # Attempt to extract the family name after "belongs to the"
        # This handles cases like "Some prefix; belongs to the X family."
        match_belongs = re.search(r"belongs to the\s+(.+)", family_text, re.IGNORECASE)
        if match_belongs:
            # Take the part after "belongs to the"
            processed_family_text = match_belongs.group(1)
        else:
            # If "belongs to the" is not found, try removing "Belongs to the " from the start
            # This handles cases like "Belongs to the X family."
            processed_family_text = re.sub(
                r"^Belongs to the\s+", "", family_text, flags=re.IGNORECASE
            )

        # Remove trailing period if any, and strip whitespace
        if processed_family_text.endswith("."):
            processed_family_text = processed_family_text[:-1]

        stripped_text = processed_family_text.strip()
        if stripped_text:
            entry_data["Protein families"] = (
                stripped_text[0].upper() + stripped_text[1:]
            )
        else:
            entry_data["Protein families"] = (
                stripped_text  # Assigns empty string if stripped_text is empty
            )

    # Check for "venom" in the complete tissue specificity text
    has_venom_tissue = "venom" in entry_data["Tissue specificity"].lower()

    # Check if entry meets search criteria
    meets_criteria = is_metazoa and (has_venom_tissue or has_toxin_keyword)

    return entry_data, meets_criteria


def process_swissprot_file(input_file, output_file):
    """Process the SwissProt file and extract matching entries to a TSV file."""
    entries_found = []

    try:
        if not input_file.exists():
            print("Error: Input file {} not found.".format(input_file))
            return

        # First count total entries for tqdm
        print("Counting total entries in the file...")
        total_entries = 0
        with open(input_file, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("//"):
                    total_entries += 1
        print("Found {} entries in the file.".format(total_entries))

        # Process the file with tqdm progress bar
        with open(input_file, "r", encoding="utf-8") as f:
            current_entry = []

            # Initialize tqdm progress bar
            pbar = tqdm(total=total_entries, desc="Processing entries")

            for line in f:
                if line.startswith("//"):
                    # End of entry, process it
                    pbar.update(1)

                    if current_entry:
                        entry_data, meets_criteria = parse_entry(current_entry)
                        if meets_criteria:
                            entries_found.append(entry_data)
                        current_entry = []
                else:
                    current_entry.append(line)

            pbar.close()

        # Write results to TSV file
        print(
            "Writing {} matching entries to {}...".format(
                len(entries_found), output_file
            )
        )
        with open(output_file, "w", encoding="utf-8") as out:
            # Write header
            out.write("\t".join(headers) + "\n")

            # Write data
            for entry in entries_found:
                row = [entry.get(field, "") for field in headers]
                out.write("\t".join(row) + "\n")

        print(
            "Processing complete. Found {} matching entries.".format(len(entries_found))
        )
        print("Results saved to {}".format(output_file))

    except Exception as e:
        print("Error processing file: {}".format(str(e)))


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse SwissProt data file to extract animal venom/toxin proteins.",
        epilog="Example: python parse_sprot_dat.py data/raw/201711_sprot.dat data/interim/toxprot_2017.tsv",
    )

    parser.add_argument("input_file", type=str, help="Path to the SwissProt data file")

    parser.add_argument(
        "output_file", type=str, help="Path where the output TSV file will be saved"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    input_path = Path(args.input_file)
    output_path = Path(args.output_file)

    # Print information
    print("Processing SwissProt data file: {}".format(input_path))
    print(
        "Query: (taxonomy_id:33208) AND ((cc_tissue_specificity:venom) OR (keyword:KW-0800))"
    )
    print("Output will be saved to: {}".format(output_path))

    # Process the file
    process_swissprot_file(input_path, output_path)
