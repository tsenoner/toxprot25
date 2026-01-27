# ToxProt25: Analysis of ToxProt (2005-2025)

Comparative analysis of toxin proteins from [ToxProt](https://www.uniprot.org/help/Toxins) across 20 years of UniProtKB releases (2005-2025).

## 🎯 Project Overview

This project analyzes changes in toxin-related proteins across different UniProtKB snapshots:

- **Taxonomic changes**: Species representation and new taxa emergence
- **Habitat patterns**: Marine vs terrestrial protein distributions
- **Protein families**: Classification and abundance changes
- **GO-term analysis**: Functional annotation comparisons
- **PTM analysis**: Post-translational modification patterns
- **Protein space**: 2D embedding visualization using ProtSpace

## 🚀 Quick Start

### Installation

```bash
# Clone and install
git clone <repository-url>
cd toxprot25
uv sync
```

### Data Processing Pipeline

The unified pipeline processes UniProt releases year-by-year to minimize disk usage (~2.5GB instead of ~50GB):

```bash
# Process all years (2005-2025): download → parse → clean
toxprot data pipeline

# Process specific years
toxprot data pipeline -y 2020-2025

# Individual stages (if needed)
toxprot data download              # Download releases only
toxprot data parse                 # Parse .dat files only
toxprot data clean                 # Clean parsed data only
```

**Output**: Final datasets in `data/processed/toxprot/` (CSV + FASTA for each year)

## 🔍 Data Filtering Criteria

```
(taxonomy_id:33208) AND (cc_tissue_specificity:venom) AND (reviewed:true)
```

Extracts reviewed Swiss-Prot Metazoa proteins with documented tissue expression in venom/venom glands.

## 📊 Analysis Scripts

Run analyses on processed data:

```bash
# Protein family analysis
python src/analysis/analyze_protein_families.py

# Taxonomic changes
python src/analysis/analyze_taxa.py

# Habitat patterns (marine vs terrestrial)
python src/analysis/analyze_habitat.py

# Post-translational modifications
python src/analysis/analyze_ptm.py

# GO term comparisons
python src/analysis/analyze_go_terms.py

# Protein evidence flow
python src/analysis/plot_protein_evidence_sankey.py
```

## 🧬 ProtSpace Analysis

Generate protein embeddings and 2D visualizations:

```bash
# 1. Prepare sequences
python src/protspace/generate_fasta_for_embeddings.py

# 2. Generate embeddings (requires ProtSpace installation)
python src/protspace/process_protspace.py

# 3. Create visualizations
python src/protspace/generate_plots.py

# 4. Clustering analysis
python src/protspace/analyze_clustering.py
```

See `src/protspace/README.md` for details.

## 📁 Project Structure

```
toxprot25/
├── src/
│   ├── data_processing/     # Pipeline: download, parse, clean
│   ├── analysis/            # Comparative analysis scripts
│   └── protspace/           # Embedding and visualization
├── data/
│   ├── raw/                 # Habitat mappings, PTM vocabulary
│   ├── interim/             # Intermediate parsed files
│   └── processed/           # Final datasets (CSV + FASTA)
└── tests/                   # Unit tests (pytest)
```

## 🛠️ CLI Reference

Use `-h` or `--help` with any command for detailed options:

```bash
toxprot -h
toxprot data -h
toxprot data pipeline -h
```

**Common options:**

- `-y, --years`: Specify years (e.g., `-y 2020-2025`)
- `-f, --force`: Reprocess existing files
- `--keep-dat`: Keep intermediate .dat files
- `-v, --verbose`: Enable verbose logging

## 🧪 Testing

```bash
# Run all tests
uv run pytest tests/ -v

# Run specific test module
uv run pytest tests/test_pipeline.py -v
```

## 📝 Data Sources

- **UniProtKB/SwissProt**: Historical releases (2005-2025)
- **ToxProt**: March 2025 export
- **Habitat classification**: Manual taxonomic curation
