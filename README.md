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

### 1. Data Processing Pipeline

**SwissProt Parsing** (`src/parse_sprot_dat.py`)

Extracts entries matching: `(taxonomy_id:33208) AND ((cc_tissue_specificity:venom) OR (keyword:KW-0800))`

Key steps:

- Process signal peptides for mature protein sequences
- Extract protein metadata (names, families, length, mass)
- Capture functional annotations (tissue specificity, toxic dose)
- Extract post-translational modifications

### 2. Comparative Analyses (2017 vs 2025)

**Taxonomic Analysis** (`src/taxa_plots_script.py`)

- Top taxa distribution comparison
- Newcomer orders and families identification
- Species-level counting

**Habitat Analysis** (`src/visualize_habitat_protein_changes.py`)

- Marine vs terrestrial protein family distribution
- Percentage and absolute changes
- Dual-habitat family analysis

**GO-term Analysis** (`src/analyze_go_terms.py`)

- Functional annotation comparisons
- GO term enrichment changes

**Curation Analysis** (`src/generate_family_renaming_report.py`)

- Protein family name changes tracking
- Systematic renaming identification

### 3. Visualization Capabilities

- **Stacked bar charts**: Taxonomic and protein family distributions
- **Heatmaps**: Habitat-specific changes
- **Venn diagrams**: Marine/terrestrial overlap analysis
- **Diverging bar charts**: Directional change visualization
- **ProtSpace plots**: 2D protein embedding
- **Sankey diagrams**: Protein evidence flow

## 📝 Data Sources

- **UniProtKB/SwissProt**: Nov 2017 (`201711_sprot.dat`) and Jan 2025 (`202501_sprot.dat`)
- **Habitat classification**: Manual curation from taxonomic families
