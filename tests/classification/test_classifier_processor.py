"""
Tests for classifier_processor module.
"""
import pandas as pd
import pytest
from pathlib import Path

from classification.utils.classifier_processor import (
    CentrifugeOutputProcessor,
    DiamondOutputProcessor,
    KrakenOutputProcessor,
    KrakenUniqOutputProcessor,
    count_prefix_spaces,
    protein_accession_to_taxid,
    taxid_to_description,
)


class TestCountPrefixSpaces:
    """Tests for count_prefix_spaces function."""

    def test_no_spaces(self):
        assert count_prefix_spaces("no_spaces") == 0

    def test_single_space(self):
        assert count_prefix_spaces(" single") == 1

    def test_multiple_spaces(self):
        assert count_prefix_spaces("   triple") == 3


class TestTaxidToDescription:
    """Tests for taxid_to_description function."""

    def test_none_taxid(self):
        assert taxid_to_description(None) is None

    def test_invalid_taxid(self):
        assert taxid_to_description(10**8) is None


class TestProteinAccessionToTaxid:
    """Tests for protein_accession_to_taxid function."""

    def test_none_accession(self):
        result = protein_accession_to_taxid("invalid_accession_xyz")
        assert result is None


@pytest.mark.unit
class TestCentrifugeOutputProcessor:
    """Tests for CentrifugeOutputProcessor class."""

    def test_from_file(self, sample_centrifuge_report: Path):
        processor = CentrifugeOutputProcessor.from_file(
            str(sample_centrifuge_report)
        )
        assert processor is not None
        assert isinstance(processor.report, pd.DataFrame)

    def test_process(self, sample_centrifuge_report: Path):
        processor = CentrifugeOutputProcessor.from_file(
            str(sample_centrifuge_report)
        )
        processor.process()
        assert not processor.final_report.empty
        assert 'description' in processor.final_report.columns
        assert 'taxID' in processor.final_report.columns

    def test_prep_final_report(self, sample_centrifuge_report: Path):
        processor = CentrifugeOutputProcessor.from_file(
            str(sample_centrifuge_report)
        )
        processor.process()
        processor.prep_final_report()
        assert processor.final_report['description'].str.contains('_').any()

    def test_save(self, sample_centrifuge_report: Path, temp_dir: Path):
        processor = CentrifugeOutputProcessor.from_file(
            str(sample_centrifuge_report)
        )
        processor.process()
        output_path = temp_dir / "output.tsv"
        processor.save(str(output_path))
        assert output_path.exists()


@pytest.mark.unit
class TestKrakenOutputProcessor:
    """Tests for KrakenOutputProcessor class."""

    def test_from_file(self, sample_kraken_report: Path):
        processor = KrakenOutputProcessor.from_file(str(sample_kraken_report))
        assert processor is not None
        assert len(processor.nodes) > 0

    def test_kraken_report_to_tree(self):
        data = {
            'PercReads': [0.0, 10.5, 50.5],
            'NumReadsRoot': [0, 1000, 500],
            'Nreads': [0, 1000, 500],
            'RankCode': ['-', '-', 'D'],
            'taxID': [1, 2, 2],
            'name': ['root', 'cellular organisms', 'Bacteria']
        }
        df = pd.DataFrame(data)
        nodes, edges = KrakenOutputProcessor.kraken_report_to_tree(df)
        assert len(nodes) > 0

    @pytest.mark.skip(reason="Known issue with Kraken tree processing")
    def test_process(self, sample_kraken_report: Path):
        processor = KrakenOutputProcessor.from_file(str(sample_kraken_report))
        processor.process()
        assert not processor.final_report.empty

    def test_summarize_leaves(self):
        leaves_simple = {
            (562, "Escherichia coli", 0, 50.0, 50, "S"): [],
            (28901, "Salmonella", 0, 30.0, 30, "S"): [],
        }
        result = KrakenOutputProcessor.summarize_leaves(leaves_simple)
        assert isinstance(result, pd.DataFrame)


@pytest.mark.unit
class TestKrakenUniqueOutputProcessor:
    """Tests for KrakenUniqueOutputProcessor class."""

    def test_from_file(self, sample_krakenunique_report: Path):
        processor = KrakenUniqOutputProcessor.from_file(
            str(sample_krakenunique_report)
        )
        assert processor is not None

    def test_process(self, sample_krakenunique_report: Path):
        processor = KrakenUniqOutputProcessor(str(sample_krakenunique_report), min_uniq_reads=1)
        processor = processor.from_file(str(sample_krakenunique_report))
        processor.process()
        assert 'description' in processor.final_report.columns


@pytest.mark.unit
class TestDiamondOutputProcessor:
    """Tests for DiamondOutputProcessor class."""

    def test_from_file(self, sample_diamond_report: Path):
        processor = DiamondOutputProcessor.from_file(str(sample_diamond_report))
        assert processor is not None

    def test_process(self, sample_diamond_report: Path):
        processor = DiamondOutputProcessor(str(sample_diamond_report), min_uniq_reads=1)
        processor = processor.from_file(str(sample_diamond_report))
        processor.process()
        assert 'description' in processor.final_report.columns


@pytest.mark.unit
class TestEmptyFiles:
    """Tests for handling empty or missing files."""

    def test_centrifuge_empty_file(self, temp_dir: Path):
        empty_file = temp_dir / "empty.txt"
        empty_file.write_text("")
        processor = CentrifugeOutputProcessor.from_file(str(empty_file))
        assert processor.report.empty

    def test_kraken_empty_file(self, temp_dir: Path):
        empty_file = temp_dir / "empty.txt"
        empty_file.write_text("")
        processor = KrakenOutputProcessor.from_file(str(empty_file))
        assert processor.nodes == []
