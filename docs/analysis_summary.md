# Analysis Guide

Detailed reference for ToxProt25 analysis commands. For quick start and data processing, see the [README](../README.md).

## Overview

The analysis pipeline generates figures and statistics from processed ToxProt data, enabling decade-based comparisons (2005 → 2015 → 2025) and full timeline trend analysis.

```
data/processed/toxprot/
        │
        ├── toxprot_2005.csv
        ├── toxprot_2015.csv
        └── toxprot_2025.csv
                │
                ▼
    toxprot analysis <command>
                │
                ▼
           figures/
           ├── summary statistics
           ├── taxonomic trends
           ├── protein families
           ├── sequence lengths
           ├── habitat distributions
           ├── tissue sources
           ├── PTM frequencies
           ├── GO term coverage
           ├── protein evidence
           └── definition comparisons
```

## Global Options

All analysis commands inherit these options from the `toxprot analysis` group:

### `--definition` Filter

Controls which ToxProt entries are included in the analysis:

| Value          | Description                                             |
| -------------- | ------------------------------------------------------- |
| `venom_tissue` | Entries with venom tissue annotation **(default)**      |
| `kw_toxin`     | Entries with toxin keyword (KW-0800)                    |
| `both_only`    | Only entries matching both criteria                     |
| `all`          | All entries (venom tissue OR toxin keyword)             |

```bash
toxprot analysis summary                    # Uses venom_tissue (default)
toxprot analysis -d all summary             # Include all ToxProt entries
toxprot analysis -d kw_toxin taxa           # Only toxin keyword entries
```

### Directory Options

| Option         | Default                   | Description                   |
| -------------- | ------------------------- | ----------------------------- |
| `--data-dir`   | `data/processed/toxprot/` | Input directory with CSVs     |
| `--output-dir` | `figures/`                | Output directory for figures  |

## Command Reference

### Quick Reference

All commands output to `figures/`.

| Command            | Years            | Description                      |
| ------------------ | ---------------- | -------------------------------- |
| `pipeline`         | All              | Run all analyses sequentially    |
| `summary`          | 2005, 2015, 2025 | Dataset statistics table         |
| `taxa`             | 2005–2025        | Taxonomic distribution trends    |
| `families`         | 2005, 2015, 2025 | Protein family distributions     |
| `length`           | 2005, 2015, 2025 | Sequence length histograms       |
| `habitat`          | 2005, 2015, 2025 | Terrestrial vs marine habitats   |
| `source-tissue`    | 2005, 2015, 2025 | Tissue annotation evolution      |
| `ptm`              | 2005, 2015, 2025 | PTM frequency analysis           |
| `go`               | 2005–2025        | GO term coverage and depth       |
| `protein-evidence` | 2008, 2015, 2025 | Evidence level transitions       |
| `definitions`      | 2025             | Selection criteria comparison    |
| `protspace`        | 2025             | Protein embedding analysis       |

### Command Details

#### `pipeline`

Runs all analysis commands sequentially with default parameters. Useful for regenerating all figures after data updates.

```bash
toxprot analysis pipeline                    # Run all analyses
toxprot analysis -d all pipeline             # Use all entries (not just venom_tissue)
toxprot analysis pipeline --data-dir custom/ # Custom input directory
```

**Output:** All figure directories populated under `figures/`.

---

#### `summary`

Generates a statistics table comparing key metrics across decade snapshots: entry counts, species diversity, taxonomic coverage, and annotation completeness.

```bash
toxprot analysis summary
toxprot analysis -d all summary
toxprot analysis summary -o figures/custom
```

**Output:** `figures/dataset_summary_statistics.png`

**Metrics included:**
- Total entries
- Unique protein families
- Missing protein family annotations (count & %)
- Fragment entries (count & %)
- PTM annotations (count & %)
- Toxic dose annotations (count & %)
- Species count
- Order count

---

#### `taxa`

Analyzes taxonomic composition across all 21 years. Generates trend plots showing top taxa evolution and alluvial diagrams identifying "newcomers" — taxonomic groups appearing in later releases.

```bash
toxprot analysis taxa                        # Default: Order level
toxprot analysis taxa --level Family         # Family-level analysis
toxprot analysis taxa -l all                 # Generate all taxonomic levels
toxprot analysis taxa --skip-trend           # Skip trend plot, only alluvial
```

**Options:**

| Option         | Default | Description                                           |
| -------------- | ------- | ----------------------------------------------------- |
| `--level`      | `Order` | Taxonomic level: Phylum, Class, Order, Family, or all |
| `--skip-trend` | False   | Skip generating trend plot                            |

**Output:** `figures/`
- `top_taxa_trend.png` — Top 5 orders over all 21 years with silhouettes (Squamata/Cobra, Araneae/Spider, Neogastropoda/Conus, Scorpiones/Scorpion, Hymenoptera)
- `taxa_newcomers_alluvial_{level}.png` — Decade-step flow diagram (2005→2015→2025)

**Notes:**
- "Newcomers" are taxa present in a later year but absent in an earlier year
- Trend plot uses all years 2005-2025; alluvial uses decade snapshots

---

#### `families`

Compares protein family distributions across decade snapshots. Shows which toxin families dominate and how their relative abundance changes over time.

```bash
toxprot analysis families
toxprot analysis families --top-n 15         # Show top 15 families
toxprot analysis -d all families
```

**Options:**

| Option    | Default | Description                        |
| --------- | ------- | ---------------------------------- |
| `--top-n` | `10`    | Number of top families to display  |

**Output:** `figures/top_families_alluvial.png` — Alluvial plot showing rank changes

**Notes:**
- Family names are normalized across years using 80+ mappings (e.g., "Snake toxin family" → "Snake three-finger toxin family")
- This enables accurate cross-year comparison despite naming changes

---

#### `length`

Generates sequence length distribution histograms comparing 2005, 2015, and 2025 datasets. Useful for understanding how the size profile of characterized toxins has evolved.

```bash
toxprot analysis length
toxprot analysis -d all length
toxprot analysis length -o figures/custom
```

**Output:** `figures/sequence_length_distribution.png`

**Visualization:**
- Overlaid histograms (not stacked) with 25 AA bins
- Bins: 1-25, 26-50, 51-75, ..., 276-300, 301+
- Latest year plotted at back for visibility

---

#### `habitat`

Analyzes terrestrial vs marine habitat distributions. Creates a two-panel figure showing taxa by habitat and protein family evolution across environments.

```bash
toxprot analysis habitat
toxprot analysis habitat --top-n 20          # Show top 20 families per habitat
toxprot analysis -d all habitat
```

**Options:**

| Option    | Default | Description                           |
| --------- | ------- | ------------------------------------- |
| `--top-n` | `15`    | Number of top families per habitat    |

**Output:** `figures/habitat.png`

**Panel A (Taxa by Habitat):**
- Left: Stacked bars for Entries, Species, Protein Families (terrestrial/marine/shared)
- Right: Breakdown of exclusive vs shared family contributions
- Flow connections showing data distribution

**Panel B (Dual-Habitat Families):**
- Diverging horizontal bar chart (marine LEFT, terrestrial RIGHT)
- Top 15 families appearing in BOTH habitats
- Three color shades per habitat showing 2005→2015→2025 evolution

---

#### `source-tissue`

Tracks source tissue annotations over time. Shows which tissues (venom gland, skin secretion, etc.) are most commonly annotated and how this changes.

```bash
toxprot analysis source-tissue
toxprot analysis source-tissue --top-n 5     # Show top 5 tissues
toxprot analysis -d all source-tissue
```

**Options:**

| Option    | Default | Description                          |
| --------- | ------- | ------------------------------------ |
| `--top-n` | `10`    | Number of top tissues to display     |

**Output:** `figures/source_tissue_alluvial.png`

**Notes:**
- Tissues are exploded from semicolon-separated values in the source data

---

#### `ptm`

Analyzes post-translational modification frequencies across decade snapshots. Shows which PTMs are most common in toxins and how annotation completeness has improved.

```bash
toxprot analysis ptm
toxprot analysis ptm --years 2010,2020       # Custom year comparison
toxprot analysis -d all ptm
```

**Options:**

| Option    | Default          | Description                      |
| --------- | ---------------- | -------------------------------- |
| `--years` | `2005,2015,2025` | Comma-separated years to compare |

**Data source:** PTM annotations extracted from UniProt XML `<feature>` elements:
- Feature types: `modified residue`, `glycosylation site`, `disulfide bond`, `cross-link`, `lipid moiety-binding region`
- Descriptions resolved using UniProt's [ptmlist.txt](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/ptmlist.txt)

**Output:** `figures/ptm_overview.png`

**Panel A:** Horizontal bar chart showing PTM type frequency, sorted by 2025 count

**Panel B:** 2×3 grid of histograms showing per-protein PTM count distributions (1-9, 10+) for top 6 types

**PTM types tracked:** Disulfide bond, Amidation, Glycosylation, Hydroxylation, Pyrrolidone carboxylic acid, Gamma-carboxyglutamic acid, D-amino acid, Bromination, Sulfation, Lipidation

---

#### `go`

Analyzes Gene Ontology term distributions across the full timeline. Shows GO term coverage, annotation depth, and category trends (molecular function, biological process, cellular component).

```bash
toxprot analysis go
toxprot analysis go --top-n 10               # Show top 10 GO terms
toxprot analysis -d all go
```

**Options:**

| Option    | Default | Description                       |
| --------- | ------- | --------------------------------- |
| `--top-n` | `5`     | Number of top GO terms to display |

**Output:** `figures/go_terms_overview.png`

**Panels:**
- **Panel A:** Total GO term annotations over time (3 lines: MF, BP, CC)
- **Panel B:** GO category coverage percentages (entries with ≥1 annotation)
- **Panels C-E:** Top 5 GO terms per category with evolution trends

**Notes:**
- Uses GO hierarchy (`go-basic.obo`) for term relationships
- Analyzes all 21 years (2005-2025), not just decade snapshots

---

#### `protein-evidence`

Tracks protein existence (PE) level transitions over time. Shows how evidence quality has improved as more toxins receive experimental validation.

```bash
toxprot analysis protein-evidence
toxprot analysis -d all protein-evidence
toxprot analysis protein-evidence -o figures/custom
```

**Years:** Uses 2008, 2015, 2025 (not 2005) because PE levels were introduced in UniProt 10.0 (March 2007).

**PE Levels:**
1. Evidence at protein level (experimental)
2. Evidence at transcript level
3. Inferred from homology
4. Predicted
5. Uncertain

**Output:** `figures/protein_evidence_sankey.png`

Alluvial diagram showing:
- PE category bars at each time point
- Flows between categories showing transitions
- "Removed" intermediate nodes showing proteins dropped from dataset
- Which PE categories lost the most proteins over time

---

#### `definitions`

Compares the two ToxProt selection criteria (venom tissue vs toxin keyword) for the 2025 dataset. Shows overlap and entries unique to each criterion.

```bash
toxprot analysis definitions
toxprot analysis definitions --year 2024     # Analyze different year
toxprot analysis definitions -o figures/custom
```

**Options:**

| Option   | Default | Description              |
| -------- | ------- | ------------------------ |
| `--year` | `2025`  | Year of dataset to use   |

**Output:** `figures/definition_comparison.png`

Single-panel Venn→phylum→order flow figure showing entry-level criteria (left), phyla (centre), and orders (right). The Venn diagram shows overlap between venom tissue specificity and toxin keyword criteria, with flows through taxonomic levels.

## Output Directory Structure

```
figures/
├── dataset_summary_statistics.png
├── definition_comparison.png
├── go_terms_overview.png
├── habitat.png
├── protein_evidence_sankey.png
├── protspace_silhouette_comparison.csv
├── protspace_silhouette_comparison.png
├── ptm_overview.png
├── sequence_length_distribution.png
├── source_tissue_alluvial.png
├── taxa_newcomers_alluvial_order.png
├── top_families_alluvial.png
└── top_taxa_trend.png
```

## ProtSpace Embedding Analysis

Protein language model (ProtT5) embeddings enable 2D visualization of toxin relationships. Three variants are compared:

| Variant        | Description                            | Silhouette Score |
| -------------- | -------------------------------------- | ---------------- |
| `full`         | Full-length sequences                  | 0.262            |
| `mature`       | Signal peptides removed                | 0.397            |
| `mature_clean` | Mature sequences, fragments excluded   | 0.474            |

Removing signal peptides and fragments improves clustering quality (silhouette: 0.262 → 0.474).

| Subcommand      | Description                                    |
| ---------------- | ---------------------------------------------- |
| `generate-fasta` | Create FASTA files for embedding generation    |
| `prepare`        | Prepare metadata and filter H5 files           |
| `run-umap`       | Run UMAP dimensionality reduction              |
| `silhouette`     | Analyze clustering quality                     |
| `pipeline`       | Run all steps (except Colab)                   |

UMAP visualizations are created manually via [protspace.app](https://protspace.app/explore).

See the [ProtSpace Guide](protspace.md) for the full pipeline and workflow.

## See Also

- [Data Processing Guide](data_processing.md) — Download, parse, and clean pipeline
- [UniProt Release History](uniprot_releases.md) — Swiss-Prot release versions
- [ProtSpace Guide](protspace.md) — Protein embedding analysis pipeline
