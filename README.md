# ToxProt25

This repository tracks how the [ToxProt](https://www.uniprot.org/help/Toxins) dataset -- UniProt's curated collection of animal toxin proteins -- has evolved across 20 years of Swiss-Prot releases (2005-2025).

It provides a reproducible pipeline that downloads historical Swiss-Prot releases, extracts toxin proteins, and produces per-year datasets (CSV + FASTA) enriched with taxonomy and habitat data. On top of that, it includes analyses of taxonomic distribution, protein families, habitat patterns, PTMs, GO terms, and protein space embeddings.

## Quick Start

```bash
git clone https://github.com/tsenoner/toxprot25.git
cd toxprot25
uv sync

# Run the full pipeline (2005-2025): download -> parse -> clean
toxprot data pipeline

# Or process specific years only
toxprot data pipeline -y 2020-2025
```

Output: `data/processed/toxprot/` (one CSV + FASTA per year).

Run any command with `-h` for full options (e.g. `-f` to force reprocess, `--keep-dat` to retain raw files, `-v` for verbose logging).

## How Proteins Are Selected

The parser extracts reviewed Metazoa proteins from each Swiss-Prot release that match **either** of two criteria (union):

```
(taxonomy_id:33208) AND ((cc_tissue_specificity:venom OR keyword:KW-0800)) AND (reviewed:true)
```

- **Venom tissue** -- the protein's tissue specificity annotation mentions venom or venom glands
- **Toxin keyword** -- the protein carries UniProt keyword KW-0800 ("Toxin")

Every output row includes a `ToxProt definition` column recording which criterion matched: `venom_tissue`, `kw_toxin`, or `both`.

## Analysis

```bash
toxprot analysis taxa                            # Taxonomic distribution plots
toxprot analysis taxa --definition venom_tissue  # Only venom-tissue entries
toxprot analysis families                        # Protein family plots
toxprot analysis families --top-n 20             # Custom top-N families
```

The `--definition` flag filters by selection criterion: `all` (default), `venom_tissue`, `kw_toxin`, or `both_only`.

Additional standalone scripts in `src/analysis/` cover habitat, PTM, and GO-term analyses. ProtSpace embedding analysis lives in `src/protspace/` (see its own README).

## Data Sources

- [UniProtKB/Swiss-Prot](https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases) -- historical releases (2005-2025)
- [ToxProt](https://www.uniprot.org/help/Toxins) -- UniProt animal toxin annotation program
- Habitat classification -- manual taxonomic curation
