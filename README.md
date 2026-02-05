# ToxProt25

Analysis of [ToxProt](https://www.uniprot.org/help/Toxins) — UniProt's curated collection of animal toxin proteins — across 21 years of Swiss-Prot releases (2005-2025).

This project provides a reproducible pipeline that downloads historical Swiss-Prot releases, extracts toxin proteins, and produces per-year datasets (CSV + FASTA) enriched with taxonomy and habitat data. It enables decade-based comparisons (2005 → 2015 → 2025) and full timeline trend analysis.

## Data Flow

```
UniProt FTP ──download──▶ XML ──parse──▶ TSV ──clean──▶ CSV + FASTA
                                                              │
                                                              ▼
                                                    toxprot analysis
                                                              │
                                                              ▼
                                                          Figures
```

## Quick Start

```bash
git clone https://github.com/tsenoner/toxprot25.git
cd toxprot25
uv sync

toxprot data pipeline           # Download and process data (2005-2025)
toxprot analysis pipeline       # Generate all figures
```

Output: `data/processed/toxprot/` (CSV + FASTA per year), `figures/` (all analysis plots).

## CLI Overview

### Data Processing (`toxprot data`)

Process UniProt Swiss-Prot releases to extract toxin proteins.

| Command    | Description                              |
| ---------- | ---------------------------------------- |
| `pipeline` | Full workflow (download → parse → clean) |
| `download` | Fetch releases from UniProt FTP          |
| `parse`    | Extract toxin proteins from XML          |
| `clean`    | Add taxonomy, habitat, create CSV/FASTA  |

**Common options:** `-y` years (e.g., `2020-2025`), `-f` force reprocess, `--keep-raw` retain XML files, `-v` verbose.

See the [Data Processing Guide](docs/data_processing.md) for detailed command options and output descriptions.

### Analysis (`toxprot analysis`)

Generate figures and statistics from processed data.

| Command            | Description                             |
| ------------------ | --------------------------------------- |
| `pipeline`         | Run all analyses with defaults          |
| `summary`          | Dataset statistics (entries, species)   |
| `taxa`             | Taxonomic distribution trends           |
| `families`         | Protein family distributions            |
| `length`           | Sequence length distributions           |
| `habitat`          | Terrestrial vs. marine distribution     |
| `source-tissue`    | Source tissue annotations               |
| `ptm`              | Post-translational modification trends  |
| `go`               | GO term distributions                   |
| `protein-evidence` | Protein existence level evolution       |
| `definitions`      | Compare selection criteria (2025 only)  |

**Default filter:** `--definition venom_tissue`. Use `--definition all` for all entries.

See the [Analysis Guide](docs/analysis_summary.md) for detailed command options and output descriptions.

## ToxProt Selection Criteria

The parser replicates UniProt's ToxProt query, extracting reviewed Metazoa proteins matching **either** criterion:

```
(taxonomy_id:33208) AND (cc_tissue_specificity:venom OR keyword:KW-0800) AND (reviewed:true)
```

| Criterion         | UniProt Field                  | Match Type                              |
| ----------------- | ------------------------------ | --------------------------------------- |
| **Venom tissue**  | `cc_tissue_specificity:venom`  | Free-text: `"venom" in CC_TISSUE_SPECIFICITY` |
| **Toxin keyword** | `keyword:KW-0800`              | Exact match on keyword element          |

The `--definition` flag filters outputs: `all`, `venom_tissue` (default), `kw_toxin`, `both_only`.

## Documentation

- [Data Processing Guide](docs/data_processing.md) — Download, parse, and clean pipeline
- [Analysis Guide](docs/analysis_summary.md) — Figure generation commands
- [UniProt Release History](docs/uniprot_releases.md) — Swiss-Prot release versions (2005-2025)
- [ProtSpace README](src/protspace/README.md) — Protein embedding analysis

## Data Sources

- [ToxProt](https://www.uniprot.org/help/Toxins) — UniProt animal toxin annotation program
- Habitat classification — manual taxonomic curation
