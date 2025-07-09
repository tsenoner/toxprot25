# ToxProt25: Comparative Analysis of Venom/Toxin Proteins (2017-2025)

A comprehensive bioinformatics project analyzing the evolution and changes in venom/toxin protein datasets between 2017 and 2025 UniProtKB/SwissProt releases.

## 🎯 Project Overview

This project provides a systematic comparative analysis of toxin-related proteins from two major UniProtKB/SwissProt database snapshots (2017 and 2025). The analysis focuses on:

- **Taxonomic changes**: Evolution of species representation and emergence of new taxa
- **Habitat-specific patterns**: Marine vs terrestrial protein family distributions
- **Protein family dynamics**: Changes in family classifications and abundance
- **Data curation insights**: Identification of protein family renamings and annotations
- **Protein space visualization**: 2D embedding analysis using ProtSpace

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

- Extracts entries matching: `(taxonomy_id:33208) AND ((cc_tissue_specificity:venom) OR (keyword:KW-0800))`
- Processes signal peptides to extract mature protein sequences
- Generates comprehensive metadata including:
  - Protein characteristics (names, families, length, mass)
  - Functional annotations (tissue specificity, toxic dose)
  - Post-translational modifications

### 2. Comparative Analyses (2017 vs 2025)

**Taxonomic Analysis** (`src/taxa_plots_script.py`)

- Top taxa distribution comparison
- Identification of newcomer orders and families
- Species-level counting improvements

**Habitat-Specific Analysis** (`src/visualize_habitat_protein_changes.py`)

- Marine vs terrestrial protein family distribution
- Percentage and absolute change calculations
- Dual-habitat family analysis with heatmaps and diverging bar charts

**Data Curation Analysis** (`src/generate_family_renaming_report.py`)

- Tracks protein family name changes between releases
- Identifies systematic renamings and improvements

### 3. Visualization Capabilities

- **Stacked bar charts**: Taxonomic and protein family distributions
- **Heatmaps**: Habitat-specific percentage changes
- **Venn diagrams**: Overlap analysis between marine/terrestrial habitats
- **Diverging bar charts**: Change visualization with clear directionality
- **ProtSpace plots**: 2D protein embedding visualization

## 🛠️ Technical Stack

- **Core**: Python 3.12+, pandas, numpy
- **Visualization**: matplotlib, seaborn, matplotlib-venn
- **Taxonomic analysis**: taxopy
- **Protein analysis**: ProtSpace
- **Notebooks**: Jupyter Lab
- **Data formats**: CSV, TSV, FASTA, JSON

## 📝 Data Sources

- **UniProtKB/SwissProt**: November 2017 release (`201711_sprot.dat`)
- **UniProtKB/SwissProt**: January 2025 release (`202501_sprot.dat`)
- **Habitat classification**: Manual curation based on taxonomic families
