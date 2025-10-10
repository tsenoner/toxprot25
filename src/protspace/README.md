# Protspace Analysis Pipeline

This directory contains scripts for generating protein embedding visualizations using ProtSpace to analyze ToxProt 2025 clustering quality.

## Quick Start

### Prerequisites (One-time setup)

```bash
# Step 0: Generate FASTA files for embedding
uv run python src/protspace/generate_fasta_for_embeddings.py

# Step 0.1: Generate ProtT5 embeddings using Colab
# Open: https://colab.research.google.com/github/tsenoner/protspace/blob/main/examples/notebook/ClickThrough_GenerateEmbeddings.ipynb
# Upload: toxprot_2025_full.fasta and toxprot_2025_mature.fasta
# Download: toxprot_2025_full.h5 (~34MB) and toxprot_2025_mature.h5 (~19MB)
# Place H5 files in: data/processed/protspace/
```

### Main Pipeline

```bash
# Step 1: Prepare metadata variants
uv run python src/protspace/prepare_protspace_metadata.py

# Step 2: Create filtered H5 variants
uv run python src/protspace/create_h5_variants.py

# Step 3: Generate ProtSpace JSON files with UMAP embeddings
uv run python src/protspace/process_protspace.py

# Step 4: Generate visualization plots
uv run python src/protspace/generate_plots.py

# Step 5: Analyze clustering quality
uv run python src/protspace/analyze_clustering.py
```

## Overview

The pipeline processes ToxProt 2025 data through four variants:

1. **All data** - Full dataset with top 15 protein families + NaN + Other (8,055 proteins)
2. **Top 15 full** - Top 15 families, full sequences with signal peptides (4,254 proteins)
3. **Top 15 mature** - Top 15 families, mature sequences, signal peptides removed (4,254 proteins)
4. **Top 15 mature no fragments** - Top 15 families, mature, fragments removed (3,399 proteins)

### Key Results

| Variant                    | Proteins | Silhouette Score | Improvement |
| -------------------------- | -------- | ---------------- | ----------- |
| Top 15 full                | 4,254    | 0.25             | Baseline    |
| Top 15 mature              | 4,254    | 0.45             | +80%        |
| Top 15 mature no fragments | 3,399    | 0.49             | +96%        |

**Conclusion:** Removing signal peptides + fragments improves protein family clustering in pLM embdding space quality.

---

## Detailed Documentation

### Pipeline Scripts

### 0. FASTA Generation for Embeddings (Prerequisite)

```bash
uv run python src/protspace/generate_fasta_for_embeddings.py
```

**Purpose:** Generate FASTA files for ProtT5 embedding generation.

**Inputs:**

- `data/interim/toxprot_2025.tsv` - Interim data with sequences and signal peptide information

**Outputs:**

- `data/processed/protspace/toxprot_2025_full.fasta` - Full sequences (8,055 sequences)
- `data/processed/protspace/toxprot_2025_mature.fasta` - Mature sequences (4,779 with signal peptides removed)

**Statistics:**

- Total sequences: 8,055
- With signal peptides: 4,779 (59.4%)
- Without signal peptides: 3,276 (40.6%)

**Note:** For sequences without signal peptides, the "mature" version is identical to the "full" version.

### 0.1 ProtT5 Embedding Generation

**Purpose:** Generate ProtT5 embeddings for both full and mature sequences.

**Method:** Use the ProtSpace Colab notebook: [ClickThrough_GenerateEmbeddings.ipynb](https://colab.research.google.com/github/tsenoner/protspace/blob/main/examples/notebook/ClickThrough_GenerateEmbeddings.ipynb)

**Process:**

1. **Upload FASTA files to Colab:**

   - `toxprot_2025_full.fasta`
   - `toxprot_2025_mature.fasta`

2. **Run embedding generation in Colab:**

   ```python
   # In the Colab notebook, set your FASTA paths
   fasta_file_full = "toxprot_2025_full.fasta"
   fasta_file_mature = "toxprot_2025_mature.fasta"

   # Generate embeddings (uses ProtT5-XL-UniRef50)
   # Output will be .h5 files containing per-protein embeddings
   ```

3. **Download generated H5 files:**

   - `toxprot_2025_full.h5` (~34MB) - Full sequence embeddings
   - `toxprot_2025_mature.h5` (~19MB) - Mature sequence embeddings

4. **Place H5 files in project:**
   ```bash
   mv toxprot_2025_full.h5 data/processed/protspace/
   mv toxprot_2025_mature.h5 data/processed/protspace/
   ```

**Embedding Details:**

- **Model:** ProtT5-XL-UniRef50 (pre-trained protein language model)
- **Embedding dimension:** 1024 per protein
- **Processing:** Mean pooling of per-residue embeddings
- **Format:** HDF5 files with protein identifiers as keys

**H5 File Usage:**
| H5 File | Used For | Description |
|---------|----------|-------------|
| `toxprot_2025_full.h5` | Variants 1 & 2 | Embeddings from full sequences (with signal peptides) |
| `toxprot_2025_mature.h5` | Variants 3 & 4 | Embeddings from mature sequences (signal peptides cut off) |

**Important:** All `.h5` files are gitignored due to size. Upload to Zenodo/Figshare for data sharing.

### 1. Metadata Preparation

```bash
uv run python src/protspace/prepare_protspace_metadata.py
```

**Purpose:** Generate four metadata CSV files from the base ToxProt 2025 dataset.

**Inputs:**

- `data/processed/toxprot_2025.csv` - Processed ToxProt data
- `data/interim/toxprot_2025.tsv` - Interim data with signal peptide information

**Outputs:**

- `metadata_2025_all.csv` (8,055 entries) - All proteins with top 15 + Other + NaN
- `metadata_2025_top15.csv` (4,254 entries) - Top 15 families, full sequences
- `metadata_2025_top15_mature.csv` (4,254 entries) - Top 15 families, mature sequences
- `metadata_2025_top15_mature_no_fragments.csv` (3,399 entries) - Top 15 families, mature, no fragments

**Column Structure:**

- Variants 1-3: `identifier`, `Protein families`, `Phylum`, `has_fragment`, `has_signal_peptide`
- Variant 4: `identifier`, `Protein families`, `Phylum`, `has_signal_peptide`

### 2. H5 Variant Creation

```bash
uv run python src/protspace/create_h5_variants.py
```

**Purpose:** Create filtered H5 embedding files for each metadata variant.

**Inputs:**

- `toxprot_2025_full.h5` - Full sequence embeddings (from Colab)
- `toxprot_2025_mature.h5` - Mature sequence embeddings (from Colab)
- Metadata CSV files from step 1

**Outputs:**

- `data/processed/protspace/toxprot_2025_all.h5`
- `data/processed/protspace/toxprot_2025_top15.h5`
- `data/processed/protspace/toxprot_2025_top15_mature.h5`
- `data/processed/protspace/toxprot_2025_top15_mature_no_fragments.h5`

**Note:** H5 files are gitignored and should be uploaded to Zenodo/Figshare.

### 3. ProtSpace JSON Generation

```bash
uv run python src/protspace/process_protspace.py
```

**Purpose:** Generate ProtSpace JSON files with UMAP embeddings and apply styling.

**Inputs:**

- H5 files from step 2
- Metadata CSV files from step 1
- `data/processed/protspace/style.json` - Color/shape styling configuration

**Outputs:**

- `protspace_2025_all.json` + `protspace_2025_all_style.json`
- `protspace_2025_top15.json` + `protspace_2025_top15_style.json`
- `protspace_2025_top15_mature.json` + `protspace_2025_top15_mature_style.json`
- `protspace_2025_top15_mature_no_fragments.json` + `protspace_2025_top15_mature_no_fragments_style.json`

**Parameters:**

- UMAP: `n_neighbors=50`, `min_dist=0.5`
- Methods: `pca2` and `umap2`

### 4. Plot Generation

```bash
uv run python src/protspace/generate_plots.py
```

**Purpose:** Generate UMAP visualization plots for all four variants.

**Requirements:** Updated `protspace` package with `ProtSpace` class.

**Inputs:**

- Styled JSON files from step 3

**Outputs:**

- `figures/protspace/2025_all_data.png`
- `figures/protspace/2025_top15_full_sequences.png`
- `figures/protspace/2025_top15_mature_sequences.png`
- `figures/protspace/2025_top15_mature_no_fragments.png`

### 5. Silhouette Score Analysis

```bash
uv run python src/protspace/analyze_clustering.py
```

**Purpose:** Calculate and compare clustering quality across variants using silhouette scores.

**Inputs:**

- JSON files from step 3

**Outputs:**

- `figures/protspace/silhouette_comparison.csv` - Numerical results
- `figures/protspace/silhouette_comparison.png` - Visualization

**Results:** Higher silhouette scores indicate better-separated clusters. Our analysis shows:

- Top 15 full: 0.25
- Top 15 mature: 0.45
- Top 15 mature no fragments: 0.49

This demonstrates that using mature sequences and removing fragments significantly improves clustering quality:

- **Mature sequences**: 80% improvement over full sequences
- **Fragment removal**: Additional 9% improvement

---

## Final Notes

**Data Availability:**

- All `.h5` files are gitignored due to file size constraints
- H5 embedding files should be uploaded to Zenodo/Figshare for data sharing
- FASTA files can be regenerated from the interim TSV if needed

**Top 15 Protein Families:**

1. Phospholipase A2 family (543)
2. Snake three-finger toxin family (536)
3. Long (4 C-C) scorpion toxin superfamily (426)
4. Neurotoxin 10 (Hwtx-1) family (290)
5. Venom metalloproteinase (M12B) family (286)
6. Conotoxin O1 superfamily (281)
7. Short scorpion toxin superfamily (274)
8. Venom Kunitz-type family (261)
9. Peptidase S1 family (232)
10. Arthropod phospholipase D family (220)
11. Conotoxin A superfamily (186)
12. Neurotoxin 19 (CSTX) family (185)
13. Conotoxin M superfamily (183)
14. Snaclec family (183)
15. Non-disulfide-bridged peptide (NDBP) superfamily (168)

## Dependencies

- `protspace` - ProtSpace visualization library
- `h5py` - HDF5 file handling
- `pandas` - Data manipulation
- `numpy` - Numerical operations
- `scikit-learn` - Silhouette score calculation
- `matplotlib` - Plotting (for silhouette analysis)
