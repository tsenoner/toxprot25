# ToxProt25: Analysis of ToxProt (2005-2025)

Comparative analysis of toxin proteins from [ToxProt](https://www.uniprot.org/help/Toxins) across 20 years of UniProtKB releases (2005-2025).

## 🎯 Project Overview

This project analyzes changes in toxin-related proteins across different UniProtKB snapshots:

- **Taxonomic changes**: Species representation and new taxa emergence
- **Habitat patterns**: Marine vs terrestrial protein family distributions
- **Protein families**: Classification and abundance changes
- **GO-term analysis**: Functional annotation comparisons
- **PTM analysis**: Post-translational modification patterns
- **Curation insights**: Protein family renamings and annotation improvements
- **Protein space**: 2D embedding visualization using ProtSpace

## 🚀 Setup

### Installation

```bash
# Install dependencies
uv sync
```

### Data Files

UniProtKB SwissProt DAT files (`.dat`) are not included due to size. Download them:

```bash
# Download all years (2005-2025)
uv run python src/data_processing/download_uniprot_releases.py
```

Files are saved to `data/raw/` as `{year}_sprot.dat`.

## 🔬 Analysis Pipeline

### 1. Data Processing

**Parse SwissProt DAT files** (`src/data_processing/parse_sprot_dat.py`)

Extracts toxin entries matching: `(taxonomy_id:33208) AND ((cc_tissue_specificity:venom) OR (keyword:KW-0800))`

Features:

- Signal peptide processing for mature sequences
- Protein metadata extraction (names, families, length, mass)
- Functional annotations (tissue specificity, toxic dose)
- Post-translational modifications (MOD_RES, CARBOHYD, DISULFID, CROSSLNK, LIPID)
- Automatic PTM vocabulary download (ptmlist.txt)

**Clean and enrich data** (`src/data_processing/clean_data.py`)

- Standardize protein family names
- Generate FASTA files (full and mature sequences)
- Taxonomic enrichment (phylum, class, order, family, genus, species)
- Habitat classification (marine/terrestrial)

### 2. Comparative Analyses

All analysis scripts are in `src/analysis/`:

| Analysis              | Script                               | Outputs                                                    |
| --------------------- | ------------------------------------ | ---------------------------------------------------------- |
| **Protein Families**  | `analyze_protein_families.py`        | Distribution charts, length histograms                     |
| **Taxonomic Changes** | `analyze_taxa.py`                    | Taxa distribution, newcomers by order/family               |
| **Habitat Patterns**  | `analyze_habitat.py`                 | Marine vs terrestrial comparisons, Venn diagrams, heatmaps |
| **PTMs**              | `analyze_ptm.py`                     | Modification type distributions and statistics             |
| **GO Terms**          | `analyze_go_terms.py`                | Functional annotation comparisons                          |
| **Protein Evidence**  | `plot_protein_evidence_sankey.py`    | Evidence flow Sankey diagrams                              |
| **Curation Tracking** | `generate_family_renaming_report.py` | Family name change reports                                 |

### 3. ProtSpace Analysis

Generate protein embeddings and 2D visualizations (`src/protspace/`):

1. `generate_fasta_for_embeddings.py` - Prepare sequences
2. `process_protspace.py` - Generate embeddings (requires ProtSpace)
3. `generate_plots.py` - Create 2D visualizations
4. `analyze_clustering.py` - Clustering quality metrics

See `src/protspace/README.md` for details.

## 📊 Output Structure

```
data/
  ├── raw/           # Source data (habitat classifications, PTM vocabulary)
  ├── interim/       # Parsed SwissProt data
  └── processed/     # Cleaned datasets and ProtSpace embeddings
figures/             # All generated visualizations
notebook/            # Jupyter analysis notebooks
docs/                # Analysis summaries and notes
```

## 📝 Data Sources

- **UniProtKB/SwissProt**: Historical and current releases (.dat files)
- **ToxProt**: March 2025 export (`202503_ToxProt.tsv`)
- **Habitat classification**: Manual taxonomic curation (`marine_terrestrial.json`, `habitat_detailed.json`)
