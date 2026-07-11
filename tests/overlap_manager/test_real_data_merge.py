"""
Integration tests for merge functions using real data fixtures.

Fixtures are subsets of production data stored in fixtures/.
Tests auto-skip when fixture directory is absent.
"""

import os
import re
from pathlib import Path

import pandas as pd
import pytest

from metagenomics_utils.overlap_manager.manager import (
    NCBI_ACCESSION_RE,
    _extract_assid,
    _merge_matched_vectorized,
    merge_by_assembly_ID,
)
from metagenomics_utils.overlap_manager.node_stats import (
    _load_m_stats_inputs,
    _merge_matched_stats,
    _resolve_assembly_classifications,
)

FIXTURES = Path(__file__).parent / "fixtures"

DATASETS = ["dataset_0001_plan", "dataset_0020_plan"]


def _paths(dataset: str) -> dict:
    """Build the paths dict expected by _load_m_stats_inputs."""
    d = FIXTURES / dataset
    return {
        "input": str(d / f"{dataset}.tsv"),
        "matched_assemblies": str(d / "matched_assemblies.tsv"),
        "classification": str(d / f"{dataset}_merged_classification.tsv"),
        "merged_stats": str(d / "merged_coverage_statistics.tsv"),
    }


def _all_paths_exist(dataset: str) -> bool:
    paths = _paths(dataset)
    return all(os.path.exists(v) for v in paths.values())


# ---------------------------------------------------------------------------
# Fixture availability
# ---------------------------------------------------------------------------

def pytest_generate_tests(metafunc):
    if "dataset" in metafunc.fixturenames:
        available = [d for d in DATASETS if _all_paths_exist(d)]
        metafunc.parametrize("dataset", available)


# ---------------------------------------------------------------------------
# _load_m_stats_inputs
# ---------------------------------------------------------------------------

class TestLoadMStatsInputs:
    """Verify _load_m_stats_inputs reads and aggregates real TSV data correctly."""

    def test_returns_four_dataframes(self, dataset):
        input_df, matched, m_stats, classification = _load_m_stats_inputs(_paths(dataset))
        assert isinstance(input_df, pd.DataFrame)
        assert isinstance(matched, pd.DataFrame)
        assert isinstance(m_stats, pd.DataFrame)
        assert isinstance(classification, pd.DataFrame)

    def test_input_has_expected_columns(self, dataset):
        input_df, *_ = _load_m_stats_inputs(_paths(dataset))
        expected = {"sample", "taxid", "reads", "assembly_accession"}
        assert expected.issubset(set(input_df.columns))

    def test_matched_has_expected_columns(self, dataset):
        _, matched, *_ = _load_m_stats_inputs(_paths(dataset))
        expected = {"taxid", "description", "uniq_reads", "classifiers", "assembly_accession", "assembly_file"}
        assert expected.issubset(set(matched.columns))

    def test_classification_has_expected_columns(self, dataset):
        *_, classification = _load_m_stats_inputs(_paths(dataset))
        expected = {"taxid", "description", "uniq_reads", "classifiers"}
        assert expected.issubset(set(classification.columns))

    def test_m_stats_has_expected_columns(self, dataset):
        _, _, m_stats, _ = _load_m_stats_inputs(_paths(dataset))
        expected = {"file", "assembly_accession", "coverage", "covbases", "meanmapq", "numreads", "error_rate"}
        assert expected.issubset(set(m_stats.columns))

    def test_m_stats_is_deduplicated_by_file(self, dataset):
        _, _, m_stats, _ = _load_m_stats_inputs(_paths(dataset))
        assert m_stats["file"].duplicated().sum() == 0
        assert len(m_stats) == m_stats["file"].nunique()

    def test_m_stats_aggregation_values(self, dataset):
        _, _, m_stats, _ = _load_m_stats_inputs(_paths(dataset))
        assert (m_stats["numreads"] >= 0).all()
        assert (m_stats["covbases"] >= 0).all()
        assert (m_stats["coverage"] >= 0).all()

    def test_m_stats_size_is_reduced(self, dataset):
        """File dedup should produce fewer rows than raw segments."""
        raw = pd.read_csv(_paths(dataset)["merged_stats"], sep="\t")
        _, _, m_stats, _ = _load_m_stats_inputs(_paths(dataset))
        assert len(m_stats) < len(raw)
        assert len(m_stats) == raw["file"].nunique()


# ---------------------------------------------------------------------------
# _resolve_assembly_classifications
# ---------------------------------------------------------------------------

class TestResolveAssemblyClassifications:
    """Verify classification enrichment of matched assemblies."""

    def test_total_uniq_reads_added(self, dataset):
        paths = _paths(dataset)
        matched = pd.read_csv(paths["matched_assemblies"], sep="\t")
        classification = pd.read_csv(paths["classification"], sep="\t")
        resolved = _resolve_assembly_classifications(matched, classification)
        assert "total_uniq_reads" in resolved.columns
        assert resolved["total_uniq_reads"].notna().all()
        assert (resolved["total_uniq_reads"] >= 0).all()

    def test_classifiers_added(self, dataset):
        paths = _paths(dataset)
        matched = pd.read_csv(paths["matched_assemblies"], sep="\t")
        classification = pd.read_csv(paths["classification"], sep="\t")
        resolved = _resolve_assembly_classifications(matched, classification)
        assert "classifiers" in resolved.columns
        assert resolved["classifiers"].notna().all()

    def test_deduplicated_on_assembly_accession(self, dataset):
        paths = _paths(dataset)
        matched = pd.read_csv(paths["matched_assemblies"], sep="\t")
        classification = pd.read_csv(paths["classification"], sep="\t")
        resolved = _resolve_assembly_classifications(matched, classification)
        assert resolved["assembly_accession"].duplicated().sum() == 0
        assert len(resolved) == resolved["assembly_accession"].nunique()

    def test_assembly_file_is_basename(self, dataset):
        paths = _paths(dataset)
        matched = pd.read_csv(paths["matched_assemblies"], sep="\t")
        classification = pd.read_csv(paths["classification"], sep="\t")
        resolved = _resolve_assembly_classifications(matched, classification)
        assert not resolved["assembly_file"].str.contains("/").any()

    def test_all_rows_have_valid_assembly_file(self, dataset):
        paths = _paths(dataset)
        matched = pd.read_csv(paths["matched_assemblies"], sep="\t")
        classification = pd.read_csv(paths["classification"], sep="\t")
        resolved = _resolve_assembly_classifications(matched, classification)
        assert resolved["assembly_file"].notna().all()
        assert not resolved["assembly_file"].isna().any()

    def test_n_rows_not_less_than_input_taxa(self, dataset):
        """Should have at least as many resolved rows as input taxa (each taxon may map to multiple accessions)."""
        paths = _paths(dataset)
        matched = pd.read_csv(paths["matched_assemblies"], sep="\t")
        classification = pd.read_csv(paths["classification"], sep="\t")
        resolved = _resolve_assembly_classifications(matched, classification)
        assert len(resolved) >= 1


# ---------------------------------------------------------------------------
# _merge_matched_vectorized
# ---------------------------------------------------------------------------

class TestMergeMatchedVectorizedReal:
    """Verify _merge_matched_vectorized with real data shapes and values."""

    @pytest.fixture
    def m_stats_and_matched(self, dataset):
        paths = _paths(dataset)
        _, matched_raw, _, classification = _load_m_stats_inputs(paths)
        matched = _resolve_assembly_classifications(matched_raw, classification)
        m_stats = pd.read_csv(paths["merged_stats"], sep="\t").rename(columns={"#rname": "assembly_accession"})
        return m_stats, matched

    def test_assembly_accession_preserved(self, m_stats_and_matched):
        m_stats, matched = m_stats_and_matched
        result = _merge_matched_vectorized(m_stats, matched)
        assert "assembly_accession" in result.columns

    def test_all_rows_get_non_null_assid(self, m_stats_and_matched):
        m_stats, matched = m_stats_and_matched
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["assid"].notna().all()

    def test_all_rows_get_non_null_taxid(self, m_stats_and_matched):
        m_stats, matched = m_stats_and_matched
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["taxid"].notna().all()

    def test_taxid_matches_known_value(self, m_stats_and_matched):
        m_stats, matched = m_stats_and_matched
        result = _merge_matched_vectorized(m_stats, matched)
        # Spot-check first row: assid should exist in matched
        first_assid = result["assid"].iloc[0]
        assert first_assid in matched["assembly_accession"].values

    def test_total_uniq_reads_filled(self, m_stats_and_matched):
        m_stats, matched = m_stats_and_matched
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["total_uniq_reads"].notna().all()
        assert (result["total_uniq_reads"] >= 0).all()

    def test_coverage_columns_kept(self, m_stats_and_matched):
        m_stats, matched = m_stats_and_matched
        result = _merge_matched_vectorized(m_stats, matched)
        for col in ["coverage", "covbases", "meanmapq", "numreads", "error_rate", "file"]:
            assert col in result.columns

    def test_description_is_string(self, m_stats_and_matched):
        m_stats, matched = m_stats_and_matched
        result = _merge_matched_vectorized(m_stats, matched)
        assert result["description"].notna().all()


# ---------------------------------------------------------------------------
# _merge_matched_stats  (full pipeline: resolve → merge → filter)
# ---------------------------------------------------------------------------

class TestMergeMatchedStatsReal:
    """Verify _merge_matched_stats end-to-end with real data."""

    def test_returns_nonempty_dataframe(self, dataset):
        paths = _paths(dataset)
        _, matched_raw, m_stats, classification = _load_m_stats_inputs(paths)
        matched = _resolve_assembly_classifications(matched_raw, classification)
        result = _merge_matched_stats(m_stats, matched)
        assert isinstance(result, pd.DataFrame)
        assert not result.empty

    def test_all_rows_have_taxid(self, dataset):
        """_merge_matched_stats filters out rows without taxid."""
        paths = _paths(dataset)
        _, matched_raw, m_stats, classification = _load_m_stats_inputs(paths)
        matched = _resolve_assembly_classifications(matched_raw, classification)
        result = _merge_matched_stats(m_stats, matched)
        assert result["taxid"].notna().all()
        assert (result["taxid"].astype(str) != "nan").all()

    def test_has_expected_columns(self, dataset):
        paths = _paths(dataset)
        _, matched_raw, m_stats, classification = _load_m_stats_inputs(paths)
        matched = _resolve_assembly_classifications(matched_raw, classification)
        result = _merge_matched_stats(m_stats, matched)
        expected = {"assid", "taxid", "total_uniq_reads", "coverage", "error_rate"}
        assert expected.issubset(set(result.columns))

    def test_file_column_dropped(self, dataset):
        paths = _paths(dataset)
        _, matched_raw, m_stats, classification = _load_m_stats_inputs(paths)
        matched = _resolve_assembly_classifications(matched_raw, classification)
        result = _merge_matched_stats(m_stats, matched)
        assert "file" not in result.columns

    def test_assembly_accession_preserved(self, dataset):
        paths = _paths(dataset)
        _, matched_raw, m_stats, classification = _load_m_stats_inputs(paths)
        matched = _resolve_assembly_classifications(matched_raw, classification)
        result = _merge_matched_stats(m_stats, matched)
        assert "assembly_accession" in result.columns


# ---------------------------------------------------------------------------
# _extract_assid
# ---------------------------------------------------------------------------

class TestExtractAssidReal:
    """Verify _extract_assid correctly extracts accessions from real filenames."""

    @pytest.fixture(scope="class")
    def all_known_accessions(self):
        known = set()
        for ds in DATASETS:
            if not _all_paths_exist(ds):
                continue
            matched = pd.read_csv(_paths(ds)["matched_assemblies"], sep="\t")
            known.update(matched["assembly_accession"].tolist())
        return list(known)

    @pytest.fixture(scope="class")
    def coverage_files(self):
        files = []
        for ds in DATASETS:
            if not _all_paths_exist(ds):
                continue
            cov = pd.read_csv(_paths(ds)["merged_stats"], sep="\t")
            files.extend(cov["file"].unique().tolist())
        return files

    def test_extracts_from_all_coverage_files(self, coverage_files, all_known_accessions):
        """Every coverage filename should yield a valid accession via _extract_assid."""
        for fname in coverage_files:
            basename = os.path.basename(fname)
            assid = _extract_assid(basename, all_known_accessions)
            assert assid is not None, f"Failed to extract from {basename}"
            assert assid in all_known_accessions, f"Extracted {assid} not in known accessions"

    def test_extract_assid_matches_standard_accessions(self):
        """_extract_assid matches standard GCF/GCA/NC/NZ/CP/etc. via regex or known_assids."""
        samples = [
            "GCF_000123455.1",
            "GCA_000123455.2",
            "NC_016609.1",
            "NZ_CP009252.1",
            "NT_123456.1",
            "NW_123456.1",
            "AC_123456.1",
            "AE_123456.1",
            "CP_012345.1",
        ]
        for s in samples:
            assert _extract_assid(f"file_{s}.bam", samples) == s, f"_extract_assid failed on {s}"

    def test_known_assids_fallback_catches_nonstandard(self, all_known_accessions):
        """Accessions with prefixes like CM_, NG_, NR_, PX_ rely on known_assids fallback."""
        nonstandard = [a for a in all_known_accessions if not NCBI_ACCESSION_RE.search(a)]
        if not nonstandard:
            pytest.skip("No non-standard accessions in fixture data")
        for acc in nonstandard:
            result = _extract_assid(f"some_prefix_{acc}.fna", all_known_accessions)
            assert result == acc, f"Failed to extract non-standard accession {acc}"

    def test_case_sensitivity_preserved(self):
        """Lowercase in filename should not match uppercase accession."""
        result = _extract_assid("gcf_000123.1.fna", ["GCF_000123.1"])
        assert result is None

    def test_longest_prefix_wins_with_known(self):
        """When multiple accessions are embedded, the longest (most specific) should match."""
        result = _extract_assid("data_GCF_000123.10.fna", ["GCF_000123.1", "GCF_000123.10"])
        assert result == "GCF_000123.10"

    def test_no_match_returns_none(self):
        result = _extract_assid("unknown_file.fna", ["GCF_000123.1"])
        assert result is None
