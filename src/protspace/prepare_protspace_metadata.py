import pandas as pd
from pathlib import Path


def get_top_families_per_dataset(df, column, top_n=15):
    counts = df[column].value_counts()
    top_counts = counts.nlargest(top_n)
    other_count = counts.iloc[top_n:].sum()
    nan_count = df[column].isna().sum()

    result_series = top_counts.copy()
    if other_count > 0:
        result_series["Other"] = other_count
    if nan_count > 0:
        result_series["NaN"] = nan_count
    return result_series


def process_column(df, column, top_n=15):
    # Get value counts
    counts = df[column].value_counts()
    top_values = set(counts.nlargest(top_n).index)

    # Create a new column with processed values
    df[f"{column}_processed"] = df[column].apply(
        lambda x: x if pd.isna(x) or x in top_values else "Other"
    )
    return df


def main():
    # Define paths using pathlib
    input_file = Path("data/processed/toxprot_2025.csv")
    output_dir = Path("data/processed/protspace")
    output_file = output_dir / "metadata_2025.csv"

    # Read the input file
    df = pd.read_csv(input_file)

    # Rename Entry to identifier
    df = df.rename(columns={"Entry": "identifier"})

    # Process Protein families
    df = process_column(df, "Protein families")

    # Process Class and Order
    df = process_column(df, "Class")
    df = process_column(df, "Order")

    # Select and rename columns
    columns_to_keep = [
        "identifier",
        "Protein families_processed",
        "Fragment",
        "Phylum",
        "Habitat",
        "Class_processed",
        "Order_processed",
    ]

    # Rename processed columns
    df = df[columns_to_keep].rename(
        columns={
            "Protein families_processed": "Protein families",
            "Class_processed": "Class",
            "Order_processed": "Order",
        }
    )

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save the processed metadata
    df.to_csv(output_file, index=False)

    # Print some statistics
    print("\nTop 15 Protein families:")
    print(get_top_families_per_dataset(df, "Protein families"))
    print("\nTop 15 Classes:")
    print(get_top_families_per_dataset(df, "Class"))
    print("\nTop 15 Orders:")
    print(get_top_families_per_dataset(df, "Order"))


if __name__ == "__main__":
    main()
