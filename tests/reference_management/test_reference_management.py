"""
Tests for reference_management module.
"""

import gzip
import os
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

os.environ["NCBI_EMAIL"] = "test@example.com"

project_root = Path(__file__).parent.parent.parent


def _write_gzip_fasta(path: Path, content: str = ">seq\nACGT\n") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write(content)


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


class TestSequenceFileValidation:
    """Integrity checks for cached and downloaded FASTA files."""

    def test_valid_gzip_fasta(self, tmp_path):
        from metagenomics_utils.ncbi_tools import is_valid_sequence_file

        path = tmp_path / "valid.fna.gz"
        _write_gzip_fasta(path)
        assert is_valid_sequence_file(str(path)) is True

    def test_truncated_gzip_fasta(self, tmp_path):
        from metagenomics_utils.ncbi_tools import is_valid_sequence_file

        path = tmp_path / "truncated.fna.gz"
        _write_gzip_fasta(path, ">seq\n" + "ACGT" * 1000 + "\n")
        path.write_bytes(path.read_bytes()[:-4])
        assert is_valid_sequence_file(str(path)) is False

    def test_empty_file(self, tmp_path):
        from metagenomics_utils.ncbi_tools import is_valid_sequence_file

        path = tmp_path / "empty.fna.gz"
        path.write_bytes(b"")
        assert is_valid_sequence_file(str(path)) is False

    def test_valid_plain_fasta(self, tmp_path):
        from metagenomics_utils.ncbi_tools import is_valid_sequence_file

        path = tmp_path / "valid.fna"
        path.write_text(">seq\nACGT\n")
        assert is_valid_sequence_file(str(path)) is True

    def test_invalid_plain_fasta(self, tmp_path):
        from metagenomics_utils.ncbi_tools import is_valid_sequence_file

        path = tmp_path / "invalid.fna"
        path.write_text("not a fasta\n")
        assert is_valid_sequence_file(str(path)) is False


class TestAssemblyStoreIntegrity:
    """Local cached assemblies must be valid before they are reused."""

    def _store(self, tmp_path):
        from metagenomics_utils.reference_utils import AssemblyStore

        return AssemblyStore(str(tmp_path / "store"))

    def test_valid_local_assembly_is_returned(self, tmp_path):
        from metagenomics_utils.ncbi_tools import Passport

        store = self._store(tmp_path)
        path = tmp_path / "store" / "562" / "562_GCF_000001.1_sequence.fasta.gz"
        _write_gzip_fasta(path)

        assembly = store.retrieve_local_assembly(Passport(taxid="562", accession="GCF_000001.1"))

        assert assembly is not None
        assert assembly.file_path == str(path)

    def test_truncated_local_assembly_is_quarantined(self, tmp_path):
        from metagenomics_utils.ncbi_tools import Passport

        store = self._store(tmp_path)
        path = tmp_path / "store" / "562" / "562_GCF_000001.1_sequence.fasta.gz"
        _write_gzip_fasta(path, ">seq\n" + "ACGT" * 1000 + "\n")
        path.write_bytes(path.read_bytes()[:-4])

        assembly = store.retrieve_local_assembly(Passport(taxid="562", accession="GCF_000001.1"))

        assert assembly is None
        assert not path.exists()
        assert list(path.parent.glob("*.corrupt.*"))

    def test_invalid_local_assembly_triggers_download(self, tmp_path, monkeypatch):
        from metagenomics_utils.ncbi_tools import Passport, ReferenceData, is_valid_sequence_file

        store = self._store(tmp_path)
        path = tmp_path / "store" / "562" / "562_GCF_000001.1_sequence.fasta.gz"
        _write_gzip_fasta(path, ">seq\n" + "ACGT" * 1000 + "\n")
        path.write_bytes(path.read_bytes()[:-4])
        reference = ReferenceData(taxid="562", accession="GCF_000001.1", assembly_id="1")

        def _fake_retrieve(reference_data, output_path, gzipped=True):
            _write_gzip_fasta(Path(output_path))
            return True

        monkeypatch.setattr(store.ncbi, "retrieve_sequence_databases", _fake_retrieve)

        assembly = store.retrieve_assembly(Passport(taxid="562", accession="GCF_000001.1"), reference_data=reference)

        assert assembly is not None
        assert is_valid_sequence_file(assembly.file_path) is True
        assert list(path.parent.glob("*.corrupt.*"))


class TestAtomicDownloads:
    """Downloaded references must only appear after passing validation."""

    def test_reference_sequence_download_is_atomic(self, tmp_path, monkeypatch):
        import metagenomics_utils.ncbi_tools as ncbi

        class _Handle:
            def read(self):
                return ">seq\nACGT\n"

            def close(self):
                return None

        monkeypatch.setattr(ncbi.Entrez, "efetch", lambda **kwargs: _Handle())
        output = tmp_path / "sequence.fasta.gz"

        assert ncbi.retrieve_reference_sequence("1", str(output), gzipped=True) is True
        assert ncbi.is_valid_sequence_file(str(output)) is True
        assert not list(tmp_path.glob("*.tmp.*"))

    def test_reference_sequence_failure_leaves_no_output(self, tmp_path, monkeypatch):
        import metagenomics_utils.ncbi_tools as ncbi

        class _Handle:
            def read(self):
                raise RuntimeError("boom")

            def close(self):
                return None

        monkeypatch.setattr(ncbi.Entrez, "efetch", lambda **kwargs: _Handle())
        output = tmp_path / "sequence.fasta.gz"

        assert ncbi.retrieve_reference_sequence("1", str(output), gzipped=True) is False
        assert not output.exists()
        assert not list(tmp_path.glob("*.tmp.*"))

    def test_assembly_download_is_atomic(self, tmp_path, monkeypatch):
        import metagenomics_utils.ncbi_tools as ncbi

        class _Handle:
            def close(self):
                return None

        def _fake_run(cmd, capture_output, check):
            with gzip.open(cmd[2], "wt") as handle:
                handle.write(">seq\nACGT\n")
            return SimpleNamespace(returncode=0, stderr=b"")

        monkeypatch.setattr(ncbi.Entrez, "esummary", lambda **kwargs: _Handle())
        monkeypatch.setattr(
            ncbi.Entrez,
            "read",
            lambda handle: {
                "DocumentSummarySet": {
                    "DocumentSummary": [{"FtpPath_RefSeq": "ftp://example/GCF_000001.1", "FtpPath_GenBank": ""}]
                }
            },
        )
        monkeypatch.setattr(ncbi.subprocess, "run", _fake_run)
        output = tmp_path / "assembly.fasta.gz"

        assert ncbi.retrieve_assembly_sequence("1", str(output)) is True
        assert ncbi.is_valid_sequence_file(str(output)) is True
        assert not list(tmp_path.glob("*.tmp.*"))

    def test_invalid_assembly_download_leaves_no_output(self, tmp_path, monkeypatch):
        import metagenomics_utils.ncbi_tools as ncbi

        class _Handle:
            def close(self):
                return None

        def _fake_run(cmd, capture_output, check):
            Path(cmd[2]).write_text("not gzip")
            return SimpleNamespace(returncode=0, stderr=b"")

        monkeypatch.setattr(ncbi.Entrez, "esummary", lambda **kwargs: _Handle())
        monkeypatch.setattr(
            ncbi.Entrez,
            "read",
            lambda handle: {
                "DocumentSummarySet": {
                    "DocumentSummary": [{"FtpPath_RefSeq": "ftp://example/GCF_000001.1", "FtpPath_GenBank": ""}]
                }
            },
        )
        monkeypatch.setattr(ncbi.subprocess, "run", _fake_run)
        output = tmp_path / "assembly.fasta.gz"

        assert ncbi.retrieve_assembly_sequence("1", str(output)) is False
        assert not output.exists()
        assert not list(tmp_path.glob("*.tmp.*"))


class TestMappingReferenceCopy:
    """Mapping references must be validated before reaching the mapper."""

    def _store(self, tmp_path):
        from metagenomics_utils.reference_utils import AssemblyStore

        return AssemblyStore(str(tmp_path / "store"))

    def test_corrupt_source_is_rejected(self, tmp_path):
        store = self._store(tmp_path)
        source = tmp_path / "source.fna.gz"
        _write_gzip_fasta(source, ">seq\n" + "ACGT" * 1000 + "\n")
        source.write_bytes(source.read_bytes()[:-4])
        mapping_dir = tmp_path / "references_to_map"
        df = pd.DataFrame(
            [{"taxid": 562, "assembly_accession": "GCF_000001.1", "assembly_file": str(source)}]
        )

        with pytest.raises(RuntimeError, match="Corrupt assembly file"):
            store.setup_mapping_references(df, mapping_references_dir=str(mapping_dir))
        assert not list(mapping_dir.iterdir())

    def test_corrupt_destination_is_replaced(self, tmp_path):
        store = self._store(tmp_path)
        source = tmp_path / "source.fna.gz"
        _write_gzip_fasta(source)
        mapping_dir = tmp_path / "references_to_map"
        mapping_dir.mkdir()
        dest = mapping_dir / "562_GCF_000001.1.fna.gz"
        _write_gzip_fasta(dest, ">seq\n" + "ACGT" * 1000 + "\n")
        dest.write_bytes(dest.read_bytes()[:-4])
        df = pd.DataFrame(
            [{"taxid": 562, "assembly_accession": "GCF_000001.1", "assembly_file": str(source)}]
        )

        store.setup_mapping_references(df, mapping_references_dir=str(mapping_dir))

        from metagenomics_utils.ncbi_tools import is_valid_sequence_file

        assert is_valid_sequence_file(str(dest)) is True
        assert list(mapping_dir.glob("*.corrupt.*"))
