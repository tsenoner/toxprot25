import pandas as pd
from pathlib import Path


# --- Path Definitions ---
try:
    BASE_PATH = Path(__file__).resolve().parent.parent
except NameError:  # Fallback for interactive environments
    BASE_PATH = Path.cwd() if Path.cwd().name == "toxprot25" else Path.cwd() / ".."

DATA_PATH = BASE_PATH / "data" / "processed"
OUTPUT_DIR = BASE_PATH / "out"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)  # Ensure 'out' directory exists


# --- Helper Function ---
def get_all_protein_families_global(df, protein_family_col_name):
    """Extract global protein family counts from a given DataFrame."""
    if protein_family_col_name not in df.columns:
        print(
            f"Warning: Column '{protein_family_col_name}' not found in DataFrame for global family counting. Returning empty Series."
        )
        return pd.Series(dtype="int")

    temp_df = df.copy()
    temp_df[protein_family_col_name] = (
        temp_df[protein_family_col_name].fillna("").astype(str)
    )
    valid_families = temp_df[protein_family_col_name][
        temp_df[protein_family_col_name] != ""
    ]
    family_counts = valid_families.value_counts()
    return family_counts.sort_index()


# --- Main Analysis Function ---
def generate_renaming_report(df_2017_orig, df_2025_orig):
    """
    Analyzes potential protein family renamings by tracking changes in the primary protein family
    assigned to the same "Entry" identifier between the 2017 and 2025 datasets.
    Saves the findings to a TSV file.
    """
    print(
        "\nAnalyzing potential protein family renamings by tracking 'Entry' identifiers..."
    )

    protein_col_name = "Protein families"
    entry_col_name = "Entry"

    # Prepare 2017 data: Entry and its primary Protein Family
    df_2017_entry_family = df_2017_orig[[entry_col_name, protein_col_name]].copy()
    df_2017_entry_family[protein_col_name] = (
        df_2017_entry_family[protein_col_name].fillna("").astype(str)
    )
    df_2017_entry_family = df_2017_entry_family[
        df_2017_entry_family[protein_col_name] != ""
    ].rename(columns={protein_col_name: "Family_2017"})

    # Prepare 2025 data: Entry and its primary Protein Family
    df_2025_entry_family = df_2025_orig[[entry_col_name, protein_col_name]].copy()
    df_2025_entry_family[protein_col_name] = (
        df_2025_entry_family[protein_col_name].fillna("").astype(str)
    )
    df_2025_entry_family = df_2025_entry_family[
        df_2025_entry_family[protein_col_name] != ""
    ].rename(columns={protein_col_name: "Family_2025"})

    # Merge based on Entry to find common proteins
    merged_df = pd.merge(
        df_2017_entry_family, df_2025_entry_family, on=entry_col_name, how="inner"
    )

    # Identify entries where the family name has changed
    changed_families_df = merged_df[
        merged_df["Family_2017"] != merged_df["Family_2025"]
    ].copy()

    if changed_families_df.empty:
        print(
            "No protein entries found with changed family names between 2017 and 2025 datasets."
        )
        print(
            "This could mean no families were renamed for proteins present in both datasets, or 'Entry' IDs are not consistent."
        )
        return

    print(
        f"Found {len(changed_families_df)} instances of 'Entry' identifiers with changed protein family names."
    )

    # Count how many entries support each specific renaming
    renaming_counts = (
        changed_families_df.groupby(["Family_2017", "Family_2025"])
        .size()
        .reset_index(name="Entries_Supporting_Rename")
    )

    # Get global counts for context (using the original full dataframes)
    counts_2017_global = get_all_protein_families_global(df_2017_orig, protein_col_name)
    counts_2025_global = get_all_protein_families_global(df_2025_orig, protein_col_name)

    renaming_counts["Old_Family_Global_Count_2017"] = renaming_counts[
        "Family_2017"
    ].apply(lambda x: counts_2017_global.get(x, 0))
    renaming_counts["New_Family_Global_Count_2025"] = renaming_counts[
        "Family_2025"
    ].apply(lambda x: counts_2025_global.get(x, 0))

    renaming_counts = renaming_counts.rename(
        columns={
            "Family_2017": "Old_Family_Name_2017",
            "Family_2025": "New_Family_Name_2025",
        }
    )

    renaming_counts = renaming_counts.sort_values(
        by=["Entries_Supporting_Rename", "Old_Family_Name_2017"],
        ascending=[False, True],
    ).reset_index(drop=True)

    output_filename = "protein_family_renaming.csv"
    output_filepath = OUTPUT_DIR / output_filename
    renaming_counts.to_csv(output_filepath, index=False)

    print(
        f"Generated potential family renamings (by Entry tracking) file: {output_filepath}"
    )
    print(
        "This file lists protein families that changed for the same 'Entry' ID between datasets."
    )
    print(
        "Review this file to understand actual remappings. 'Entries_Supporting_Rename' indicates how many unique protein entries support each listed change.\n"
    )


# --- Main Execution ---
if __name__ == "__main__":
    print(f"Project BASE_PATH resolved to: {BASE_PATH}")
    print(f"Data will be loaded from: {DATA_PATH}")
    print(f"Renaming report will be saved to: {OUTPUT_DIR}")

    try:
        toxprot_2017_df_orig = pd.read_csv(DATA_PATH / "toxprot_2017.csv")
        toxprot_2025_df_orig = pd.read_csv(DATA_PATH / "toxprot_2025.csv")
    except FileNotFoundError as e:
        print(
            f"Error loading original data: {e}. Please ensure CSV files are in {DATA_PATH}"
        )
        exit()

    generate_renaming_report(
        toxprot_2017_df_orig,
        toxprot_2025_df_orig,
    )
    print("\nFamily renaming report generation complete.")
