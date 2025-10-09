# ToxProt25: Analysis of ToxProt (2017-2025)

Comparative analysis of toxin proteins from [ToxProt](https://www.uniprot.org/help/Toxins) between UniProtKB 2017 and 2025 releases.

## 🎯 Project Overview

This project analyzes changes in toxin-related proteins across two UniProtKB snapshots (2017 and 2025):

- **Taxonomic changes**: Species representation and new taxa emergence
- **Habitat patterns**: Marine vs terrestrial protein family distributions
- **Protein families**: Classification and abundance changes
- **GO-term analysis**: Functional annotation comparisons
- **PTM analysis**: Post-translational modification patterns
- **Curation insights**: Protein family renamings and annotation improvements
- **Protein space**: 2D embedding visualization using ProtSpace

## 🚀 Quick Start

### Prerequisites

- Python 3.12+
- UV package manager (recommended) or pip

### Installation

```bash
# Clone the repository
git clone <repository-url>
cd toxprot25

# Install dependencies with UV (recommended)
uv sync

# OR install with pip
pip install -e .
```

## 🔬 Key Analyses

### 1. Data Processing

**Step 1: Parse SwissProt DAT files** (`src/data_processing/parse_sprot_dat.py`)

Extracts entries matching: `(taxonomy_id:33208) AND ((cc_tissue_specificity:venom) OR (keyword:KW-0800))`

Key features:

- Process signal peptides for mature protein sequences
- Extract protein metadata (names, families, length, mass)
- Capture functional annotations (tissue specificity, toxic dose)
- Extract post-translational modifications (MOD_RES, CARBOHYD, DISULFID, CROSSLNK, LIPID)
- Auto-download PTM controlled vocabulary (ptmlist.txt)

**Step 2: Clean and enrich data** (`src/data_processing/clean_data.py`)

- Standardize protein family names
- Create FASTA files with signal peptide removal
- Add taxonomic information (phylum, class, order, family, genus, species)
- Add habitat classification (marine/terrestrial)

**Step 3: Remove fragments** (`src/data_processing/remove_fragments.py`)

- Filter out fragment sequences for specific analyses

### 2. Comparative Analyses (2017 vs 2025)

**Protein Families** (`src/analysis/analyze_protein_families.py`)

- Distribution comparisons with stacked bar charts
- Sequence length histograms
- Summary statistics tables

**Taxonomic Changes** (`src/analysis/analyze_taxa.py`)

- Top taxa distribution comparison
- Newcomer orders and families identification
- Species-level counting

**Habitat Patterns** (`src/analysis/analyze_habitat.py`)

- Marine vs terrestrial protein family distribution
- Percentage and absolute changes
- Dual-habitat family analysis
- Venn diagrams and heatmaps

**Post-Translational Modifications** (`src/analysis/analyze_ptm.py`)

- PTM type distributions
- Modification count analysis
- Comparative statistics

**GO Terms** (`src/analysis/analyze_go_terms.py`)

- Functional annotation comparisons
- GO term enrichment changes

**Curation Tracking** (`src/analysis/generate_family_renaming_report.py`)

- Protein family name changes tracking
- Systematic renaming identification

**Protein Evidence** (`src/analysis/plot_protein_evidence_sankey.py`)

- Sankey diagrams showing protein evidence flow

## 📝 Data Sources

- **UniProtKB/SwissProt**: Nov 2017 (`201711_sprot.dat`) and Jan 2025 (`202501_sprot.dat`)
- **Habitat classification**: Manual curation from taxonomic families
