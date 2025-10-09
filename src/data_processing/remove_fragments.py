import pandas as pd
from pathlib import Path

# --- Path Definitions ---
try:
    BASE_PATH = Path(__file__).resolve().parent.parent
except NameError:  # Fallback for interactive environments
    BASE_PATH = Path.cwd()

DATA_PATH = BASE_PATH / "data" / "processed" / "protspace"
INPUT_CSV = DATA_PATH / "metadata_2025.csv"
OUTPUT_CSV = DATA_PATH / "metadata_2025_no_fragments.csv"

# Specific column to check for the exact value "fragment"
TARGET_COLUMN_FOR_FRAGMENT_CHECK = "Fragment"


def remove_fragments_from_metadata(input_file: Path, output_file: Path):
    """
    Reads a protein metadata CSV, removes rows where the 'Fragment' column
    has the exact value "fragment", and saves the result to a new CSV file.

    Args:
        input_file (Path): Path to the input CSV file.
        output_file (Path): Path to save the cleaned CSV file.
    """
    if not input_file.exists():
        print(f"Error: Input file not found at {input_file}")
        return

    print(f"Reading data from {input_file}...")
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    original_row_count = len(df)
    print(f"Original number of rows: {original_row_count}")

    df_cleaned = df.copy()  # Start with a copy to modify

    if TARGET_COLUMN_FOR_FRAGMENT_CHECK in df.columns:
        # Ensure the target column is treated as string for exact matching, handling potential NaNs or other types gracefully.
        # Rows where the column is NA will not match "fragment" and will be kept.
        df_cleaned = df[df[TARGET_COLUMN_FOR_FRAGMENT_CHECK].astype(str).ne("fragment")]
        print(
            f"Checked column '{TARGET_COLUMN_FOR_FRAGMENT_CHECK}' for exact value 'fragment'."
        )
    else:
        print(
            f"Warning: Column '{TARGET_COLUMN_FOR_FRAGMENT_CHECK}' not found in the CSV. No rows removed based on this criterion."
        )

    cleaned_row_count = len(df_cleaned)
    rows_removed = original_row_count - cleaned_row_count

    print(f"Number of rows removed: {rows_removed}")
    print(f"Number of rows after cleaning: {cleaned_row_count}")

    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        df_cleaned.to_csv(output_file, index=False)
        print(f"Cleaned data saved to {output_file}")
    except Exception as e:
        print(f"Error writing cleaned data to CSV: {e}")


if __name__ == "__main__":
    remove_fragments_from_metadata(INPUT_CSV, OUTPUT_CSV)
