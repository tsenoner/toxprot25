#!/usr/bin/env python3
"""Analyze and visualize GO term distributions in ToxProt datasets.

Generates multiple visualizations showing GO term growth patterns,
including hierarchy-aware count propagation to reveal annotation
refinement over time.
"""

import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import obonet
import pandas as pd

from ..config import ALL_YEARS, DATA_DIR, FIGURES_DIR
from .helpers import load_datasets

# Configuration
TOP_N_GO_TERMS = 15
GO_COLUMN = "Gene Ontology (molecular function)"
GO_BP_COLUMN = "Gene Ontology (biological process)"
GO_CC_COLUMN = "Gene Ontology (cellular component)"
DEFAULT_OBO_PATH = Path("data/external/go-basic.obo")

# GO root terms (Molecular Function, Biological Process, Cellular Component)
GO_ROOTS = {"GO:0003674", "GO:0008150", "GO:0005575"}


def load_go_hierarchy(obo_path: Path = DEFAULT_OBO_PATH) -> nx.DiGraph:
    """Load GO hierarchy as NetworkX graph."""
    return obonet.read_obo(str(obo_path))


def extract_go_id(term_str: str) -> str | None:
    """Extract GO ID from a term string like 'GO:0008200 (name)'."""
    match = re.match(r"(GO:\d+)", term_str.strip())
    return match.group(1) if match else None


def get_ancestors(graph: nx.DiGraph, go_id: str) -> set[str]:
    """Get all ancestor terms (parents, grandparents, etc.)."""
    try:
        return nx.descendants(graph, go_id)
    except nx.NetworkXError:
        return set()


def get_term_depth(graph: nx.DiGraph, go_id: str) -> int | None:
    """Calculate depth of a GO term (distance from root)."""
    if go_id not in graph:
        return None

    min_depth = None
    for root in GO_ROOTS:
        if root in graph:
            try:
                path_length = nx.shortest_path_length(graph, go_id, root)
                if min_depth is None or path_length < min_depth:
                    min_depth = path_length
            except nx.NetworkXNoPath:
                continue

    return min_depth


def get_go_counts(df: pd.DataFrame, column: str = GO_COLUMN) -> dict[str, int]:
    """Extract and count GO terms from a DataFrame."""
    counts: dict[str, int] = defaultdict(int)
    if column not in df.columns:
        return dict(counts)
    for entry in df[column].dropna().astype(str):
        for term_str in entry.split(";"):
            go_id = extract_go_id(term_str)
            if go_id:
                counts[go_id] += 1
    return dict(counts)


def get_go_counts_with_propagation(df: pd.DataFrame, go_graph: nx.DiGraph) -> dict[str, int]:
    """Count GO terms with ancestor propagation."""
    counts: dict[str, int] = defaultdict(int)
    for entry in df[GO_COLUMN].dropna().astype(str):
        for term_str in entry.split(";"):
            go_id = extract_go_id(term_str)
            if go_id and go_id in go_graph:
                counts[go_id] += 1
                for ancestor in get_ancestors(go_graph, go_id):
                    counts[ancestor] += 1
    return dict(counts)


def get_go_term_name(go_graph: nx.DiGraph, go_id: str) -> str:
    """Get human-readable name for a GO term."""
    if go_id in go_graph:
        return go_graph.nodes[go_id].get("name", go_id)
    return go_id


def get_go_statistics(df: pd.DataFrame, go_graph: nx.DiGraph) -> dict:
    """Calculate comprehensive GO statistics for a dataset."""
    entries_with_go = df[GO_COLUMN].notna().sum() if GO_COLUMN in df.columns else 0
    total_entries = len(df)

    depths = []
    all_terms = []

    if GO_COLUMN in df.columns:
        for entry in df[GO_COLUMN].dropna().astype(str):
            for term_str in entry.split(";"):
                go_id = extract_go_id(term_str)
                if go_id:
                    all_terms.append(go_id)
                    depth = get_term_depth(go_graph, go_id)
                    if depth is not None:
                        depths.append(depth)

    return {
        "entries_with_go": entries_with_go,
        "total_entries": total_entries,
        "coverage_pct": (entries_with_go / total_entries * 100) if total_entries > 0 else 0,
        "total_annotations": len(all_terms),
        "unique_terms": len(set(all_terms)),
        "average_depth": np.mean(depths) if depths else None,
        "depths": depths,
    }


def load_datasets_with_all_go(
    years: list[int],
    interim_dir: Path = Path("data/interim/toxprot_parsed"),
) -> dict[int, pd.DataFrame]:
    """Load ToxProt datasets from interim TSV files (have all GO categories)."""
    datasets = {}
    for year in years:
        filepath = interim_dir / f"toxprot_{year}.tsv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath, sep="\t")
    return datasets


# =============================================================================
# Visualization Functions
# =============================================================================


def plot_go_overview(
    datasets: dict[int, pd.DataFrame],
    go_graph: nx.DiGraph,
    output_path: Path,
    top_n: int = 5,
):
    """Create multipanel overview figure with key GO insights.

    Panels:
    A. GO category coverage (%)
    B. Depth sum + average depth
    C. Top MF terms
    D. Top BP terms
    E. Top CC terms
    """
    years = sorted(datasets.keys())

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 6, hspace=0.35, wspace=0.6)

    # Tick years for x-axis
    tick_years = [2005, 2010, 2015, 2020, 2025]

    # ==========================================================================
    # Panel A: Total GO Term Counts
    # ==========================================================================
    ax_a = fig.add_subplot(gs[0, 0:3])  # First half of top row

    mf_coverage = []
    bp_coverage = []
    cc_coverage = []
    mf_term_counts = []  # Total number of GO terms (not entries)
    bp_term_counts = []
    cc_term_counts = []

    for year in years:
        if year in datasets:
            df = datasets[year]
            total = len(df)
            # Coverage: entries with at least one annotation
            mf_entries = df[GO_COLUMN].notna().sum() if GO_COLUMN in df.columns else 0
            bp_entries = df[GO_BP_COLUMN].notna().sum() if GO_BP_COLUMN in df.columns else 0
            cc_entries = df[GO_CC_COLUMN].notna().sum() if GO_CC_COLUMN in df.columns else 0
            mf_coverage.append(mf_entries / total * 100 if total > 0 else 0)
            bp_coverage.append(bp_entries / total * 100 if total > 0 else 0)
            cc_coverage.append(cc_entries / total * 100 if total > 0 else 0)

            # Total GO term counts (sum of all terms across all entries)
            mf_terms = (
                sum(len(str(x).split(";")) for x in df[GO_COLUMN].dropna())
                if GO_COLUMN in df.columns
                else 0
            )
            bp_terms = (
                sum(len(str(x).split(";")) for x in df[GO_BP_COLUMN].dropna())
                if GO_BP_COLUMN in df.columns
                else 0
            )
            cc_terms = (
                sum(len(str(x).split(";")) for x in df[GO_CC_COLUMN].dropna())
                if GO_CC_COLUMN in df.columns
                else 0
            )
            mf_term_counts.append(mf_terms)
            bp_term_counts.append(bp_terms)
            cc_term_counts.append(cc_terms)
        else:
            mf_coverage.append(0)
            bp_coverage.append(0)
            cc_coverage.append(0)
            mf_term_counts.append(0)
            bp_term_counts.append(0)
            cc_term_counts.append(0)

    ax_a.plot(
        years,
        mf_term_counts,
        marker="o",
        linewidth=2.5,
        color="#e74c3c",
        label="Molecular Function (MF)",
        markersize=6,
    )
    ax_a.plot(
        years,
        bp_term_counts,
        marker="s",
        linewidth=2.5,
        color="#3498db",
        label="Biological Process (BP)",
        markersize=6,
    )
    ax_a.plot(
        years,
        cc_term_counts,
        marker="^",
        linewidth=2.5,
        color="#2ecc71",
        label="Cellular Component (CC)",
        markersize=6,
    )

    ax_a.set_xlabel("Year", fontsize=12)
    ax_a.set_ylabel("Total GO Term Annotations", fontsize=12)
    ax_a.set_xticks(tick_years)
    ax_a.grid(True, ls="--", alpha=0.4)
    ax_a.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x:,.0f}"))
    ax_a.legend(loc="upper left", fontsize=10)

    ax_a.set_title("A. Total GO Term Annotations", fontsize=13, fontweight="bold")

    # ==========================================================================
    # Panel B: GO Category Coverage (%)
    # ==========================================================================
    ax_b = fig.add_subplot(gs[0, 3:6])  # Second half of top row

    ax_b.plot(
        years,
        mf_coverage,
        marker="o",
        linewidth=2.5,
        color="#e74c3c",
        label="Molecular Function (MF)",
        markersize=6,
    )
    ax_b.plot(
        years,
        bp_coverage,
        marker="s",
        linewidth=2.5,
        color="#3498db",
        label="Biological Process (BP)",
        markersize=6,
    )
    ax_b.plot(
        years,
        cc_coverage,
        marker="^",
        linewidth=2.5,
        color="#2ecc71",
        label="Cellular Component (CC)",
        markersize=6,
    )

    ax_b.set_xlabel("Year", fontsize=12)
    ax_b.set_ylabel("Coverage (%)", fontsize=12)
    ax_b.set_ylim(0, 105)
    ax_b.set_xticks(tick_years)
    ax_b.grid(True, ls="--", alpha=0.4)
    ax_b.legend(loc="lower right", fontsize=10)

    ax_b.set_title("B. GO Category Coverage", fontsize=13, fontweight="bold")

    # ==========================================================================
    # Helper function for top terms panels
    # ==========================================================================
    # Distinct colors for top terms (colorblind-friendly)
    term_colors = [
        "#1f77b4",  # blue
        "#ff7f0e",  # orange
        "#2ca02c",  # green
        "#d62728",  # red
        "#9467bd",  # purple
        "#8c564b",  # brown
        "#e377c2",  # pink
        "#7f7f7f",  # gray
        "#bcbd22",  # olive
        "#17becf",  # cyan
    ]

    def plot_top_terms(ax, datasets, go_column, go_graph, years, title, top_n=5):
        """Plot top GO terms for a category."""
        # Get years with data
        years_with_data = [y for y in years if y in datasets]

        # Get GO counts for each year
        go_counts = {}
        for year in years_with_data:
            counts = {}
            df = datasets[year]
            if go_column in df.columns:
                for entry in df[go_column].dropna().astype(str):
                    for term_str in entry.split(";"):
                        go_id = extract_go_id(term_str)
                        if go_id:
                            counts[go_id] = counts.get(go_id, 0) + 1
            go_counts[year] = counts

        if not years_with_data:
            return

        # Get top terms from latest year
        latest = max(years_with_data)
        sorted_terms = sorted(go_counts[latest].items(), key=lambda x: x[1], reverse=True)
        top_terms = [t[0] for t in sorted_terms[:top_n] if t[0] in go_graph]

        for i, go_id in enumerate(top_terms):
            counts = [go_counts[y].get(go_id, 0) for y in years_with_data]
            name = get_go_term_name(go_graph, go_id)
            ax.plot(
                years_with_data,
                counts,
                marker="o",
                markersize=4,
                linewidth=2,
                label=name,
                color=term_colors[i % len(term_colors)],
            )

        ax.set_xlabel("Year", fontsize=11)
        ax.set_ylabel("Annotations", fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.legend(loc="upper left", fontsize=7, framealpha=0.95)
        ax.set_xticks(tick_years)
        ax.grid(True, ls="--", alpha=0.4)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x:,.0f}"))

    # ==========================================================================
    # Panel C: Top MF Terms
    # ==========================================================================
    ax_c = fig.add_subplot(gs[1, 0:2])
    plot_top_terms(ax_c, datasets, GO_COLUMN, go_graph, years, "C. Top MF Terms", top_n)

    # ==========================================================================
    # Panel D: Top BP Terms
    # ==========================================================================
    ax_d = fig.add_subplot(gs[1, 2:4])
    plot_top_terms(ax_d, datasets, GO_BP_COLUMN, go_graph, years, "D. Top BP Terms", top_n)

    # ==========================================================================
    # Panel E: Top CC Terms
    # ==========================================================================
    ax_e = fig.add_subplot(gs[1, 4:6])
    plot_top_terms(ax_e, datasets, GO_CC_COLUMN, go_graph, years, "E. Top CC Terms", top_n)

    plt.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.2)
    plt.close()


def generate_all_figures(
    datasets: dict[int, pd.DataFrame],
    output_dir: Path,
    go_graph: nx.DiGraph | None = None,
    top_n: int = TOP_N_GO_TERMS,
):
    """Generate all GO term visualization figures."""
    output_dir.mkdir(parents=True, exist_ok=True)

    if go_graph is None:
        go_graph = load_go_hierarchy()

    # Main overview figure (multipanel)
    print("  Generating overview figure...")
    plot_go_overview(datasets, go_graph, output_dir / "go_terms_overview.png", top_n)


def main():
    """Main function for standalone execution."""
    from .helpers import filter_by_definition

    data_dir = DATA_DIR
    output_dir = FIGURES_DIR
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading GO hierarchy...")
    go_graph = load_go_hierarchy()

    print("Loading ToxProt datasets...")
    datasets = load_datasets(ALL_YEARS, data_dir)
    datasets = {year: filter_by_definition(df, "venom_tissue") for year, df in datasets.items()}
    print(f"  Loaded {len(datasets)} years")

    print("Generating figures...")
    generate_all_figures(datasets, output_dir, go_graph)

    print(f"Done. Figures saved to {output_dir}")


if __name__ == "__main__":
    main()
