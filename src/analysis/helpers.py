"""Shared helper functions for ToxProt analysis modules."""

from pathlib import Path

import pandas as pd

DATA_DIR = Path("data/processed/toxprot")


def load_datasets(years: list[int], data_dir: Path = DATA_DIR) -> dict[int, pd.DataFrame]:
    """Load ToxProt datasets for specified years.

    Args:
        years: List of years to load.
        data_dir: Directory containing CSV files named toxprot_{year}.csv.

    Returns:
        Dictionary mapping year to DataFrame.
    """
    datasets = {}
    for year in years:
        filepath = data_dir / f"toxprot_{year}.csv"
        if filepath.exists():
            datasets[year] = pd.read_csv(filepath)
    return datasets


def filter_by_definition(df: pd.DataFrame, definition: str = "all") -> pd.DataFrame:
    """Filter DataFrame by ToxProt definition column.

    Args:
        df: DataFrame with a 'ToxProt definition' column.
        definition: One of 'all', 'venom_tissue', 'kw_toxin', 'both_only'.

    Returns:
        Filtered DataFrame.
    """
    if definition == "all":
        return df
    if "ToxProt definition" not in df.columns:
        return df

    if definition == "venom_tissue":
        return df[df["ToxProt definition"].isin(["venom_tissue", "both"])]
    elif definition == "kw_toxin":
        return df[df["ToxProt definition"].isin(["kw_toxin", "both"])]
    elif definition == "both_only":
        return df[df["ToxProt definition"] == "both"]
    else:
        return df
