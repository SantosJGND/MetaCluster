"""
Tests for reference_management module.
"""

import os
from pathlib import Path

import pandas as pd
import pytest

os.environ["NCBI_EMAIL"] = "test@example.com"

project_root = Path(__file__).parent.parent.parent


class TestPassport:
    """Tests for Passport dataclass."""

    def test_passport_creation(self):
        from metagenomics_utils.ncbi_tools import Passport

        passport = Passport(taxid="562", accession="NC_000001")
        assert passport.taxid == "562"
        assert passport.accession == "NC_000001"

    def test_passport_prefix_with_accession(self):
        from metagenomics_utils.ncbi_tools import Passport

        passport = Passport(taxid="562", accession="NC_000001")
        assert passport.prefix == "562_NC_000001"

    def test_passport_prefix_without_accession(self):
        from metagenomics_utils.ncbi_tools import Passport

        passport = Passport(taxid="562")
        assert passport.prefix == "562"

    def test_passport_taxid_with_version(self):
        from metagenomics_utils.ncbi_tools import Passport

        passport = Passport(taxid="562.1")
        assert passport.taxid == "562"

    def test_passport_str(self):
        from metagenomics_utils.ncbi_tools import Passport

        passport = Passport(taxid="562", accession="NC_000001")
        assert "562" in str(passport)
        assert "NC_000001" in str(passport)


class TestLocalAssembly:
    """Tests for LocalAssembly dataclass."""

    def test_local_assembly_creation(self):
        from metagenomics_utils.ncbi_tools import LocalAssembly

        assembly = LocalAssembly(taxid="562", accession="NC_000001", file_path="/path/to/file.fasta")
        assert assembly.taxid == "562"
        assert assembly.file_path == "/path/to/file.fasta"


class TestReferenceData:
    """Tests for ReferenceData dataclass."""

    def test_reference_data_creation(self):
        from metagenomics_utils.ncbi_tools import ReferenceData

        ref = ReferenceData(taxid="562", accession="NC_000001", nucleotide_id="123456", assembly_id="GCF_000001")
        assert ref.taxid == "562"
        assert ref.nucleotide_id == "123456"


class TestCompareLineages:
    """Tests for compare_lineages function."""

    def test_identical_lineages(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        lineage1 = "A; B; C; D"
        lineage2 = "A; B; C; D"
        score, level = compare_lineages(lineage1, lineage2)
        assert score == 1.0

    def test_different_lineages(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        lineage1 = "A; B; C"
        lineage2 = "X; Y; Z"
        score, level = compare_lineages(lineage1, lineage2)
        assert score == 0.0

    def test_partial_match(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        lineage1 = "A; B; C; D"
        lineage2 = "A; B; X; Y"
        score, level = compare_lineages(lineage1, lineage2)
        assert score == 0.5
        assert level == "phylum"

    def test_none_lineage(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        score, level = compare_lineages(None, "A; B; C")
        assert score == 0.0
        assert level is None

    def test_empty_lineage(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        score, level = compare_lineages("", "A; B; C")
        assert score == 0.0


class TestMissingPctThreshold:
    """Tests for the max_missing_pct skip threshold."""

    def _missing_pct_exceeds(self, n_missing, n_qualified, max_pct=5.0):
        from reference_management.main import missing_pct_exceeds

        return missing_pct_exceeds(n_missing, n_qualified, max_pct)

    def test_zero_qualified_never_exceeds(self):
        assert self._missing_pct_exceeds(0, 0) is False
        assert self._missing_pct_exceeds(5, 0) is False

    def test_no_missing_never_exceeds(self):
        assert self._missing_pct_exceeds(0, 100) is False

    def test_exactly_threshold_proceeds(self):
        assert self._missing_pct_exceeds(1, 20, 5.0) is False

    def test_over_threshold_skips(self):
        assert self._missing_pct_exceeds(2, 20, 5.0) is True

    def test_custom_threshold(self):
        assert self._missing_pct_exceeds(30, 100, 25.0) is True
        assert self._missing_pct_exceeds(26, 100, 25.0) is True


class TestRetrieveSkip:
    """Retrieve must proceed under the threshold, skip (exit 3) above it."""

    @staticmethod
    def _args(tmp_path, *rows, no_fail_on_missing=False, max_missing_pct=5.0):
        from types import SimpleNamespace

        df = pd.DataFrame(rows, columns=["taxid", "uniq_reads", "assembly_accession", "assembly_file"])
        mapping = str(tmp_path / "references_to_map")

        return df, mapping, SimpleNamespace(
            input_table=str(tmp_path / "in.tsv"),
            assembly_store=str(tmp_path / "store"),
            mapping_references_dir=mapping,
            min_uniq_reads=1,
            max_missing_pct=max_missing_pct,
            no_fail_on_missing=no_fail_on_missing,
        )

    def _retrieve(self, tmp_path, monkeypatch, *rows, no_fail_on_missing=False, max_missing_pct=5.0):
        import reference_management.main as ref_main

        df, mapping, args = self._args(
            tmp_path, *rows, no_fail_on_missing=no_fail_on_missing, max_missing_pct=max_missing_pct
        )
        os.makedirs(mapping, exist_ok=True)

        class _FakeStore:
            def match_taxid_to_assembly(self, path):
                return df.copy()

            def setup_mapping_references(self, d, mapping_references_dir="references_to_map"):
                return None

        monkeypatch.setattr(ref_main, "AssemblyStore", lambda *a, **k: _FakeStore())
        try:
            ref_main.retrieve_assemblies(args)
        except SystemExit as e:
            return e.code
        return 0

    def test_below_threshold_proceeds(self, tmp_path, monkeypatch):
        code = self._retrieve(
            tmp_path,
            monkeypatch,
            [1, 5, "NC_1", "/s/1.fa.gz"],
            [2, 5, "NC_2", "/s/2.fa.gz"],
            [3, 5, None, None],  # 1 of 3 qualified unmatched -> 33% > 5% even with max 40
            max_missing_pct=40.0,
        )
        assert code == 0

    def test_no_missing_proceeds(self, tmp_path, monkeypatch):
        code = self._retrieve(
            tmp_path, monkeypatch, [1, 5, "NC_1", "/s/1.fa.gz"], [2, 5, "NC_2", "/s/2.fa.gz"]
        )
        assert code == 0

    def test_over_threshold_skips_with_exit_3(self, tmp_path, monkeypatch):
        code = self._retrieve(
            tmp_path,
            monkeypatch,
            [1, 5, "NC_1", "/s/1.fa.gz"],
            [2, 5, "NC_2", "/s/2.fa.gz"],
            [3, 5, None, None],
            [4, 5, None, None],
        )
        assert code == 3

    def test_no_fail_on_missing_never_skips(self, tmp_path, monkeypatch):
        code = self._retrieve(
            tmp_path,
            monkeypatch,
            [1, 5, "NC_1", "/s/1.fa.gz"],
            [2, 5, None, None],
            no_fail_on_missing=True,
        )
        assert code == 0

    def test_skip_writes_marker(self, tmp_path, monkeypatch):
        from reference_management.main import retrieve_assemblies

        df, mapping, args = self._args(
            tmp_path, [1, 5, "NC_1", "/s/1.fa.gz"], [2, 5, None, None]
        )
        os.makedirs(mapping, exist_ok=True)

        class _FakeStore:
            def match_taxid_to_assembly(self, path):
                return df.copy()

            def setup_mapping_references(self, d, mapping_references_dir="references_to_map"):
                return None

        monkeypatch.setattr("reference_management.main.AssemblyStore", lambda *a, **k: _FakeStore())
        with pytest.raises(SystemExit) as exc:
            retrieve_assemblies(args)
        assert exc.value.code == 3
        marker = Path(mapping) / "DATASET_SKIPPED.txt"
        assert marker.exists()
        assert "50.0%" in marker.read_text()
