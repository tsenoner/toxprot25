# CLAUDE.md — ToxProt 2025

## Project overview

ToxProt is a curated database of animal venom proteins built from UniProtKB/Swiss-Prot. This repository contains the data processing pipeline, analysis scripts, and manuscript figures for the 2025 publication.

- **Language**: Python 3.12+
- **Package manager**: uv
- **CLI entry point**: `toxprot` (defined in `src/cli.py`)
- **Build system**: hatchling
- **Linting**: ruff (line-length 100, pycodestyle + pyflakes + isort + bugbear + pyupgrade)

## Commands

```bash
# Install dependencies
uv sync

# Run tests
uv run pytest tests/ -q

# Lint
uv run ruff check src/
uv run ruff format src/

# Full ProtSpace pipeline (after H5 embeddings exist)
toxprot analysis protspace pipeline --skip-fasta

# Individual ProtSpace steps
toxprot analysis protspace generate-fasta    # Creates full/mature/active FASTAs
toxprot analysis protspace prepare           # Metadata + H5 filtering
toxprot analysis protspace run-umap          # UMAP projection
toxprot analysis protspace silhouette        # Clustering quality

# Data processing
toxprot data parse --years 2025              # Parse UniProt XML to TSV
toxprot data clean --years 2025              # Clean TSV to CSV + FASTA
```

## Repository structure

```
src/
  cli.py                          # Top-level CLI (click groups: data, analysis)
  config.py                       # Global config constants
  data_processing/
    parse_sprot_xml.py             # UniProt XML parser (primary)
    parse_sprot_dat.py             # Legacy flat-file parser
    clean_data.py                  # TSV -> CSV cleaning, FASTA generation
    pipeline.py                    # Full data pipeline orchestration
    cli.py                         # Data processing CLI commands
  analysis/
    analyze_protein_families.py    # Protein family analysis (top N families)
    analyze_sequence_length.py     # Sequence length analysis
    colors.py                      # Color scheme for protein families/phyla
    helpers.py                     # Shared analysis utilities
    protspace/
      config.py                    # ProtSpace variant configs, filename helpers
      fasta_generator.py           # FASTA generation (full/mature/active)
      metadata_preparer.py         # Metadata CSV + H5 filtering per variant
      umap_runner.py               # UMAP dimensionality reduction
      clustering_analyzer.py       # Silhouette score analysis
      style_generator.py           # ProtSpace style.json generation
      cli.py                       # ProtSpace CLI commands

data/
  raw/                             # Raw UniProt XML/TSV downloads
  interim/toxprot_parsed/          # Parsed TSVs (one per year)
  processed/
    toxprot/                       # Cleaned CSVs (one per year)
    protspace/
      colab/                       # FASTA + H5 embedding files
      intermediates/               # Auto-generated metadata + filtered H5s (gitignored)

figures/                           # Output figures and analysis results
tests/                             # pytest tests
```

## Data pipeline

1. **Parse**: UniProt XML -> `data/interim/toxprot_parsed/toxprot_{year}.tsv`
   - Filters for Metazoa organisms with venom-related annotations
   - Extracts: sequences, signal peptide ranges, propeptide ranges (multi-region, semicolon-separated), chain ranges, PTM features, GO terms, taxonomy
   - Propeptide/chain range extraction requires re-parsing XML (`toxprot data pipeline -y YEAR --force`)
2. **Clean**: TSV -> `data/processed/toxprot/toxprot_{year}.csv` + FASTA
   - Standardizes protein families, merges GO columns
   - Computes `Mature_length` (SP removed) and `Active_length` (SP + propeptide removed)
   - Retains `Signal peptide` and `Propeptide` (yes/no) columns in CSV; drops range columns
   - Adds taxonomy via taxopy (phylum, class, order, family, genus, species, habitat)
3. **Analysis**: Various analysis modules produce figures from the processed CSVs

## ProtSpace analysis

ProtSpace visualizes protein embeddings (ProtT5) in 2D using UMAP, colored by protein family.

### Sequence variants

| Variant | seq_type | Description | Silhouette |
|---------|----------|-------------|------------|
| `full` | full | Complete precursor (with SP) | 0.323 |
| `mature` | mature | Signal peptide removed | 0.474 |
| `mature_clean` | mature | SP removed, no fragments | 0.572 |
| `active` | active | SP + propeptide removed | 0.412 |
| `active_clean` | active | SP + PP removed, no fragments | 0.564 |

**Key finding**: Removing signal peptides significantly improves clustering (0.323 -> 0.474). Additional propeptide removal slightly *decreases* clustering quality (0.474 -> 0.412), suggesting propeptide regions carry family-discriminative information.

Earlier drafts of this table reported lower values (0.262/0.397/0.474/0.382/0.401). Those were computed with a `clustering_analyzer.py` exclude list that didn't cover pandas' `<NA>` literal, so the ~400 proteins without a family annotation were silently counted as their own cluster. The fix (commit adding `<NA>`/`None` to the exclude list) restores the intended semantics; UMAP coordinates and family labels are unchanged.

### Embedding generation

Embeddings can be generated via:
- **Biocentral API**: `biocentral-api` package, server at `https://biocentral.rostlab.org`
- **Google Colab**: Upload FASTAs from `data/processed/protspace/colab/`, run ProtSpace_Embeddings.ipynb

### Config architecture

Variant configs use `seq_type` field (`"full"`, `"mature"`, `"active"`) that maps to:
- FASTA file: `toxprot_{year}_{seq_type}.fasta`
- H5 file: `toxprot_{year}_{seq_type}.h5`

## Sequence length analysis

`toxprot analysis length` generates 6 figures comparing 2005/2015/2025:

| # | File | Length column | Filter | 2005 | 2015 | 2025 |
|---|------|-------------|--------|-----:|-----:|-----:|
| 1 | `sequence_length_full.png` | Length | no fragments | 1,250 | 4,916 | 6,183 |
| 2 | `sequence_length_full_sp.png` | Length | no fragments + has SP | 484 | 3,314 | 4,184 |
| 3 | `sequence_length_mature.png` | Mature_length | no fragments | 1,250 | 4,916 | 6,183 |
| 4 | `sequence_length_active.png` | Active_length | no fragments + has propeptide | 198 | 2,067 | 2,595 |
| 5 | `sequence_length_mature_peptide.png` | Active_length | no fragments | 1,250 | 4,916 | 6,183 |
| 6 | `sequence_length_mature_peptide_with_fragments.png` | Active_length | all entries | 1,379 | 6,012 | 7,418 |

- **Fig 6** matches ProtSpace `active` variant (7,418 proteins, includes fragments)
- **Fig 5** matches ProtSpace `active_clean` variant (6,183 proteins, no fragments)
- Reviewer requested "mature peptide" = bioactive form (SP + propeptide removed) = figs 5/6
- Entries without SP/propeptide annotation are assumed already mature (length unchanged)

## Key data stats (2025 dataset)

- Total entries: 8055 (7418 after venom_tissue filter)
- Has signal peptide: 4799 (4410 venom_tissue)
- Has propeptide: 3076
- Has both SP and propeptide: 2906
- Top 10 protein families cover the main analysis; rest grouped as "Other"

## Conventions

- Protein family names are normalized via `normalize_family_name()` in `analyze_protein_families.py`
- Colors for families/phyla defined in `analysis/colors.py` — shared across all figures
- `ToxProt definition` column: "venom_tissue", "both", or other — most analyses filter to venom_tissue + both
- Reference families come from `get_reference_families()` to ensure consistent ordering across analyses
