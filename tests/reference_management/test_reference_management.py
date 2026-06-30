"""
Tests for reference_management module.
"""
import os
import sys
import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import Mock, patch

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
        assembly = LocalAssembly(
            taxid="562",
            accession="NC_000001",
            file_path="/path/to/file.fasta"
        )
        assert assembly.taxid == "562"
        assert assembly.file_path == "/path/to/file.fasta"


class TestReferenceData:
    """Tests for ReferenceData dataclass."""

    def test_reference_data_creation(self):
        from metagenomics_utils.ncbi_tools import ReferenceData
        ref = ReferenceData(
            taxid="562",
            accession="NC_000001",
            nucleotide_id="123456",
            assembly_id="GCF_000001"
        )
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
