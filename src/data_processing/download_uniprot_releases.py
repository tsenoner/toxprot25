#!/usr/bin/env python3
"""
Download UniProt Swiss-Prot releases for each year from 2005 to 2025.

Downloads the first release of each year from:
https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/

The script downloads uniprot_sprot-only*.tar.gz files and extracts the .dat files.
"""

import argparse
import gzip
import logging
import shutil
import subprocess
import tarfile
from pathlib import Path

logger = logging.getLogger(__name__)

# Base URL for UniProt previous releases
BASE_URL = "https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases"

# Mapping of years to release directories and file patterns
# Source: Directory dates from https://ftp.uniprot.org/pub/databases/uniprot/previous_major_releases/
# Pre-2010 releases used version numbers (release4.0 = Feb 2005, release7.0 = Feb 2006, etc.)
# From 2010 onwards, releases use YYYY_MM format
# Format: year -> (release_dir, archive_name_pattern)
RELEASES = {
    2005: ("release4.0", "uniprot_sprot-only4.0.tar.gz"),
    2006: ("release7.0", "uniprot_sprot-only7.0.tar.gz"),
    2007: ("release10.0", "uniprot_sprot-only10.0.tar.gz"),
    2008: ("release13.0", "uniprot_sprot-only13.0.tar.gz"),
    2009: ("release15.0", "uniprot_sprot-only15.0.tar.gz"),
    2010: ("release-2010_01", "uniprot_sprot-only2010_01.tar.gz"),
    2011: ("release-2011_01", "uniprot_sprot-only2011_01.tar.gz"),
    2012: ("release-2012_01", "uniprot_sprot-only2012_01.tar.gz"),
    2013: ("release-2013_01", "uniprot_sprot-only2013_01.tar.gz"),
    2014: ("release-2014_01", "uniprot_sprot-only2014_01.tar.gz"),
    2015: ("release-2015_01", "uniprot_sprot-only2015_01.tar.gz"),
    2016: ("release-2016_01", "uniprot_sprot-only2016_01.tar.gz"),
    2017: ("release-2017_01", "uniprot_sprot-only2017_01.tar.gz"),
    2018: ("release-2018_01", "uniprot_sprot-only2018_01.tar.gz"),
    2019: ("release-2019_01", "uniprot_sprot-only2019_01.tar.gz"),
    2020: ("release-2020_01", "uniprot_sprot-only2020_01.tar.gz"),
    2021: ("release-2021_01", "uniprot_sprot-only2021_01.tar.gz"),
    2022: ("release-2022_01", "uniprot_sprot-only2022_01.tar.gz"),
    2023: ("release-2023_01", "uniprot_sprot-only2023_01.tar.gz"),
    2024: ("release-2024_01", "uniprot_sprot-only2024_01.tar.gz"),
    2025: ("release-2025_01", "uniprot_sprot-only2025_01.tar.gz"),
}


def download_file(url: str, output_path: Path) -> bool:
    """Download a file using curl with progress bar."""
    logger.debug(f"  Downloading: {url}")
    try:
        # Show progress bar only if logging at DEBUG level
        show_progress = logger.isEnabledFor(logging.DEBUG)
        curl_args = ["curl", "-L", "-o", str(output_path), url]
        if show_progress:
            curl_args.insert(2, "-#")  # Add progress bar
        else:
            curl_args.insert(2, "-s")  # Silent mode

        subprocess.run(curl_args, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"  Error downloading: {e}")
        return False


def extract_dat_file(archive_path: Path, output_dir: Path, year: int) -> Path | None:
    """Extract the .dat file from the tar.gz archive."""
    logger.debug(f"  Extracting .dat file from {archive_path.name}...")

    output_file = output_dir / f"{year}_sprot.dat"

    try:
        with tarfile.open(archive_path, "r:gz") as tar:
            # Find the uniprot_sprot.dat.gz file in the archive
            dat_gz_files = [m for m in tar.getmembers() if "sprot.dat" in m.name]

            if not dat_gz_files:
                logger.error("  Error: No sprot.dat file found in archive")
                return None

            # Extract the first matching file found
            dat_member = dat_gz_files[0]
            logger.debug(f"  Found: {dat_member.name}")

            # Extract to output directory
            tar.extract(dat_member, output_dir, filter="data")
            extracted_path = output_dir / dat_member.name

            # If it's gzipped, decompress it
            if extracted_path.suffix == ".gz":
                logger.debug(f"  Decompressing {extracted_path.name}...")
                with gzip.open(extracted_path, "rb") as f_in:
                    with open(output_file, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                extracted_path.unlink()  # Remove the .gz file
            else:
                # Just rename if not gzipped
                extracted_path.rename(output_file)

            logger.debug(f"  Saved: {output_file}")
            return output_file

    except Exception as e:
        logger.error(f"  Error extracting: {e}")
        return None


def download_release(year: int, output_dir: Path, keep_archive: bool = False) -> bool:
    """Download and extract a single year's release."""
    if year not in RELEASES:
        logger.error(f"  Error: No release mapping for year {year}")
        return False

    release_dir, archive_name = RELEASES[year]
    url = f"{BASE_URL}/{release_dir}/knowledgebase/{archive_name}"

    output_dat = output_dir / f"{year}_sprot.dat"

    # Check if already downloaded
    if output_dat.exists():
        logger.debug(f"  Already exists: {output_dat}")
        return True

    # Download archive
    archive_path = output_dir / archive_name

    if not archive_path.exists():
        if not download_file(url, archive_path):
            return False
    else:
        logger.debug(f"  Archive exists: {archive_path}")

    # Extract .dat file
    result = extract_dat_file(archive_path, output_dir, year)

    # If extraction failed, the archive might be corrupted - delete and retry once
    if result is None and archive_path.exists():
        logger.warning("  Archive appears corrupted, deleting and retrying...")
        archive_path.unlink()
        if download_file(url, archive_path):
            result = extract_dat_file(archive_path, output_dir, year)

    # Optionally remove archive to save space
    if result and not keep_archive and archive_path.exists():
        logger.debug(f"  Removing archive: {archive_path.name}")
        archive_path.unlink()

    return result is not None


def main():
    parser = argparse.ArgumentParser(
        description="Download UniProt Swiss-Prot releases (2005-2025)"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/raw"),
        help="Directory to save downloaded files (default: data/raw)",
    )
    parser.add_argument(
        "--years",
        type=int,
        nargs="+",
        default=list(range(2005, 2026)),
        help="Years to download (default: 2005-2025)",
    )
    parser.add_argument(
        "--keep-archives",
        action="store_true",
        help="Keep the .tar.gz archives after extraction",
    )
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="List available releases without downloading",
    )

    args = parser.parse_args()

    if args.list_only:
        print("Available releases:")
        for year, (release_dir, archive) in sorted(RELEASES.items()):
            print(f"  {year}: {release_dir}/knowledgebase/{archive}")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Output directory: {args.output_dir}")
    print(f"Years to download: {min(args.years)}-{max(args.years)}")
    print("=" * 60)

    success_count = 0
    for year in sorted(args.years):
        print(f"\n[{year}]")
        if download_release(year, args.output_dir, args.keep_archives):
            success_count += 1

    print("\n" + "=" * 60)
    print(f"Downloaded {success_count}/{len(args.years)} releases")


if __name__ == "__main__":
    main()
