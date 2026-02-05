# Data Processing Guide

Detailed reference for ToxProt25 data processing commands. For quick start and analysis, see the [README](../README.md).

## Overview

The data processing pipeline downloads UniProt Swiss-Prot releases (2005–2025), extracts toxin proteins matching ToxProt criteria, and produces enriched CSV + FASTA files for each year.

```
UniProt FTP ──download──▶ XML ──parse──▶ TSV ──clean──▶ CSV + FASTA
     │                      │              │              │
     │                      │              │              └─ data/processed/toxprot/
     │                      │              └─ data/interim/toxprot_parsed/
     │                      └─ data/raw/uniprot_releases/
     └─ ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/
```

## Command Reference

### Quick Reference

| Command    | Description                              | Output                         |
| ---------- | ---------------------------------------- | ------------------------------ |
| `pipeline` | Full workflow (download → parse → clean) | `data/processed/toxprot/`      |
| `download` | Fetch releases from UniProt FTP          | `data/raw/uniprot_releases/`   |
| `parse`    | Extract toxin proteins from XML          | `data/interim/toxprot_parsed/` |
| `clean`    | Add taxonomy, habitat, create CSV/FASTA  | `data/processed/toxprot/`      |

### Global Options

All commands support year selection:

| Option          | Example           | Description             |
| --------------- | ----------------- | ----------------------- |
| `-y`, `--years` | `-y 2020-2025`    | Year range              |
| `-y`, `--years` | `-y 2020 -y 2021` | Multiple specific years |
| (default)       |                   | All years (2005-2025)   |

### Command Details

#### `pipeline`

Runs the complete workflow: download → parse → clean. Processes year-by-year to minimize disk usage (<10GB instead of ~100GB if all files were kept).

```bash
toxprot data pipeline                    # Process all years (2005-2025)
toxprot data pipeline -y 2020-2025       # Process specific range
toxprot data pipeline -y 2025 --force    # Reprocess even if output exists
toxprot data pipeline --keep-raw         # Keep XML files for debugging
```

**Options:**

| Option                    | Default                       | Description                          |
| ------------------------- | ----------------------------- | ------------------------------------ |
| `--raw-dir`               | `data/raw/uniprot_releases`   | Directory for downloaded XML files   |
| `--interim-dir`           | `data/interim/toxprot_parsed` | Directory for intermediate TSV files |
| `--processed-dir`         | `data/processed/toxprot`      | Directory for final output           |
| `--keep-raw/--delete-raw` | `--delete-raw`                | Keep XML files after parsing         |
| `--keep-tsv/--delete-tsv` | `--keep-tsv`                  | Keep TSV files after cleaning        |
| `-f`, `--force`           | `--no-force`                  | Reprocess even if output exists      |
| `-v`, `--verbose`         | `--quiet`                     | Enable verbose logging               |

**Output:** One `toxprot_{year}.csv` and `toxprot_{year}.fasta` per year.

---

#### `download`

Downloads Swiss-Prot releases from the UniProt FTP archive. Uses the first release of each year.

```bash
toxprot data download                    # Download all years
toxprot data download -y 2020-2025       # Download specific range
toxprot data download --list-only        # List available releases
toxprot data download --keep-archives    # Keep .tar.gz files
```

**Options:**

| Option               | Default                     | Description                       |
| -------------------- | --------------------------- | --------------------------------- |
| `-o`, `--output-dir` | `data/raw/uniprot_releases` | Download destination              |
| `--keep-archives`    | False                       | Keep .tar.gz after extraction     |
| `--list-only`        | False                       | List releases without downloading |

**Source:** https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/

**Output:** `{year}_sprot.xml` files (~5GB each compressed, ~100GB total uncompressed).

---

#### `parse`

Parses Swiss-Prot XML files to extract toxin proteins matching ToxProt criteria.

```bash
toxprot data parse -i data/raw/uniprot_releases           # Parse all XML files
toxprot data parse -i data/raw/uniprot_releases -y 2025   # Parse specific year
toxprot data parse data/raw/uniprot_releases/2025_sprot.xml  # Parse single file
toxprot data parse -i data/raw --delete-input             # Delete XML after parsing
```

**Options:**

| Option                        | Default                       | Description                 |
| ----------------------------- | ----------------------------- | --------------------------- |
| `-i`, `--input-dir`           | (required or use positional)  | Directory with XML files    |
| `-o`, `--output-dir`          | `data/interim/toxprot_parsed` | Output directory            |
| `--delete-input/--keep-input` | `--keep-input`                | Delete XML after processing |

**Extraction criteria** (replicates UniProt ToxProt query):

- Taxonomy: Metazoa (ID: 33208)
- Reviewed: Swiss-Prot entries only
- Either:
  - Tissue specificity contains "venom" (free-text match)
  - Has keyword KW-0800 (Toxin)

**Output:** `toxprot_{year}.tsv` with columns for accession, protein name, organism, sequence, PTMs, GO terms, etc.

---

#### `clean`

Cleans parsed TSV files: adds taxonomy hierarchy, standardizes protein families, classifies habitats, and creates FASTA files.

```bash
toxprot data clean                       # Clean all parsed files
toxprot data clean -y 2020-2025          # Clean specific years
toxprot data clean -i custom/parsed      # Custom input directory
```

**Options:**

| Option               | Default                       | Description                     |
| -------------------- | ----------------------------- | ------------------------------- |
| `-i`, `--input-dir`  | `data/interim/toxprot_parsed` | Directory with TSV files        |
| `-o`, `--output-dir` | `data/processed/toxprot`      | Output directory                |
| `--data-dir`         | `data`                        | Base directory for habitat maps |

**Processing steps:**

1. Parse TSV and standardize protein family names
2. Add NCBI taxonomy hierarchy (Phylum, Class, Order, Family, Genus)
3. Classify habitats (terrestrial/marine) using curated mappings
4. Create FASTA with signal peptides removed
5. Add ToxProt definition column (venom_tissue, kw_toxin, both)

**Output:**

- `toxprot_{year}.csv` — Enriched data with taxonomy and habitat
- `toxprot_{year}.fasta` — Sequences with signal peptides removed

## Output Directory Structure

```
data/
├── raw/
│   ├── uniprot_releases/      # Downloaded XML files (deleted by default)
│   │   └── {year}_sprot.xml
│   ├── marine_terrestrial.json # Habitat classification (manual curation)
│   └── habitat_detailed.json   # Detailed habitat mappings
├── interim/
│   └── toxprot_parsed/         # Parsed TSV files
│       └── toxprot_{year}.tsv
└── processed/
    └── toxprot/                # Final output
        ├── toxprot_2005.csv
        ├── toxprot_2005.fasta
        ├── ...
        ├── toxprot_2025.csv
        └── toxprot_2025.fasta
```

## Output File Format

### CSV Columns

| Column                                | Description                                       |
| ------------------------------------- | ------------------------------------------------- |
| `Entry`                               | UniProt accession (e.g., P0DL84)                  |
| `Organism (ID)`                       | NCBI taxonomy ID                                  |
| `Protein families`                    | Standardized protein family name                  |
| `Length`                              | Sequence length (amino acids)                     |
| `Fragment`                            | Whether sequence is a fragment                    |
| `Source tissues`                      | Tissue expression annotations                     |
| `Toxic dose`                          | Toxicity information (LD50, etc.)                 |
| `PTM Keywords`                        | PTM-related UniProt keywords                      |
| `PTM Summary`                         | Parsed PTM annotations                            |
| `Gene Ontology (GO)`                  | All GO term IDs                                   |
| `Gene Ontology (biological process)`  | GO biological process terms                       |
| `Gene Ontology (cellular component)`  | GO cellular component terms                       |
| `Gene Ontology (molecular function)`  | GO molecular function terms                       |
| `Protein existence`                   | Evidence level (1-5)                              |
| `ToxProt definition`                  | Selection criterion: venom_tissue, kw_toxin, both |
| `Scientific_Name`                     | Species scientific name                           |
| `Phylum`                              | Taxonomic phylum                                  |
| `Class`                               | Taxonomic class                                   |
| `Order`                               | Taxonomic order                                   |
| `Family`                              | Taxonomic family                                  |
| `Genus`                               | Taxonomic genus                                   |
| `Species`                             | Species name                                      |
| `Habitat`                             | Terrestrial, Marine, or Freshwater                |
| `Habitat_Detailed`                    | More specific habitat classification              |

### FASTA Format

```
>P01387
TKCYVTPDVKSETCPAGQDICYTETWCDAWCTSRGKRVDLGCAATCPIVKPGVEIKCCSTDNCNPFPTWRKRP
```

Headers contain the UniProt accession. Signal peptides are removed based on UniProt annotations.

## Disk Space Requirements

| Stage       | Per Year | Total (21 years) | Notes                    |
| ----------- | -------- | ---------------- | ------------------------ |
| Archives    | ~200MB   | ~4GB             | Deleted after extraction |
| XML files   | ~2.5GB   | ~50GB            | Deleted by default       |
| TSV files   | ~1MB     | ~20MB            | Kept by default          |
| CSV + FASTA | ~1MB     | ~20MB            | Final output             |

The pipeline processes year-by-year to avoid storing all XML files simultaneously, reducing peak disk usage from ~50GB to ~2.5GB.

## Troubleshooting

### Download failures

```bash
# Retry a specific year
toxprot data download -y 2025

# Check if release exists
toxprot data download --list-only
```

### Parse errors

```bash
# Verify XML file integrity
head -1 data/raw/uniprot_releases/2025_sprot.xml

# Should show: <?xml version='1.0' encoding='UTF-8'?>
```

### Missing taxonomy data

The clean step requires internet access to query NCBI taxonomy. If taxonomy lookup fails:

```bash
# Check network connectivity
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id=9606" | head
```

## See Also

- [UniProt Release History](uniprot_releases.md) — Swiss-Prot release versions
- [Analysis Guide](analysis_summary.md) — Figure generation commands
