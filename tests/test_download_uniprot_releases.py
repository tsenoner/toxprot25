"""Tests for UniProt release download script."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
import tempfile
import gzip
import tarfile
import io

from src.data_processing.download_uniprot_releases import (
    RELEASES,
    BASE_URL,
    extract_dat_file,
    download_release,
)


class TestReleasesMapping:
    """Test the RELEASES mapping configuration."""

    def test_all_years_present(self):
        """Verify all years from 2005 to 2025 are in RELEASES."""
        expected_years = set(range(2005, 2026))
        actual_years = set(RELEASES.keys())
        assert actual_years == expected_years

    def test_release_tuple_format(self):
        """Verify each release entry has correct format (dir, archive)."""
        for year, (release_dir, archive_name) in RELEASES.items():
            assert isinstance(release_dir, str)
            assert isinstance(archive_name, str)
            assert archive_name.endswith(".tar.gz")

    def test_pre_2010_releases_use_version_numbers(self):
        """Verify pre-2010 releases use version number format."""
        for year in range(2005, 2010):
            release_dir, _ = RELEASES[year]
            assert release_dir.startswith("release") and not release_dir.startswith("release-")

    def test_post_2010_releases_use_year_format(self):
        """Verify 2010+ releases use YYYY_01 format."""
        for year in range(2010, 2026):
            release_dir, archive_name = RELEASES[year]
            assert f"release-{year}_01" == release_dir
            assert f"{year}_01" in archive_name


class TestURLConstruction:
    """Test URL construction for downloads."""

    def test_url_format(self):
        """Verify URLs are constructed correctly."""
        for year, (release_dir, archive_name) in RELEASES.items():
            expected_url = f"{BASE_URL}/{release_dir}/knowledgebase/{archive_name}"
            assert "ftp.uniprot.org" in expected_url
            assert "knowledgebase" in expected_url


class TestExtraction:
    """Test archive extraction functionality."""

    def test_extract_dat_file_from_tar_gz(self):
        """Test extracting a .dat.gz file from a tar archive."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Create a mock tar.gz archive with a dat.gz file inside
            dat_content = b"ID   TEST_ENTRY\nAC   P12345\n//"
            archive_path = tmpdir / "test_archive.tar.gz"

            # Create gzipped dat content
            dat_gz_buffer = io.BytesIO()
            with gzip.GzipFile(fileobj=dat_gz_buffer, mode="wb") as gz:
                gz.write(dat_content)
            dat_gz_content = dat_gz_buffer.getvalue()

            # Create tar.gz with the dat.gz file
            with tarfile.open(archive_path, "w:gz") as tar:
                info = tarfile.TarInfo(name="uniprot_sprot.dat.gz")
                info.size = len(dat_gz_content)
                tar.addfile(info, io.BytesIO(dat_gz_content))

            # Test extraction
            result = extract_dat_file(archive_path, tmpdir, 2020)

            assert result is not None
            assert result.exists()
            assert result.name == "2020_sprot.dat"
            assert result.read_bytes() == dat_content

    def test_extract_returns_none_on_missing_dat(self):
        """Test that extraction returns None if no dat file in archive."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            archive_path = tmpdir / "empty.tar.gz"

            # Create empty tar.gz
            with tarfile.open(archive_path, "w:gz") as tar:
                pass

            result = extract_dat_file(archive_path, tmpdir, 2020)
            assert result is None


class TestDownloadRelease:
    """Test the download_release function."""

    def test_skip_existing_file(self):
        """Test that existing files are skipped."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            existing_file = tmpdir / "2020_sprot.dat"
            existing_file.write_text("existing content")

            result = download_release(2020, tmpdir)

            assert result is True
            assert existing_file.read_text() == "existing content"

    def test_invalid_year_returns_false(self):
        """Test that invalid years return False."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = download_release(1999, Path(tmpdir))
            assert result is False
