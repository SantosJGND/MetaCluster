"""
Tests for dataframe_utils module.
"""
import pandas as pd
import pytest


from metagenomics_utils.dataframe_utils import (
    detect_column,
    detect_id_columns,
    rename_columns_to_standard,
    safe_get_column,
    setup_logger,
    validate_required_columns,
)


class TestDetectColumn:
    """Tests for detect_column function."""

    def test_detect_existing_column(self):
        df = pd.DataFrame({'taxid': [1, 2], 'name': ['a', 'b']})
        assert detect_column(df, 'taxid', 'TaxID', 'taxon') == 'taxid'

    def test_detect_second_choice(self):
        df = pd.DataFrame({'TaxID': [1, 2], 'name': ['a', 'b']})
        assert detect_column(df, 'taxid', 'TaxID', 'taxon') == 'TaxID'

    def test_detect_nonexistent_column(self):
        df = pd.DataFrame({'name': ['a', 'b']})
        assert detect_column(df, 'taxid', 'TaxID', 'taxon') is None


class TestDetectIdColumns:
    """Tests for detect_id_columns function."""

    def test_detect_taxid_and_accession(self):
        df = pd.DataFrame({
            'taxid': [1, 2],
            'assembly_accession': ['A', 'B']
        })
        result = detect_id_columns(df)
        assert result['taxid_col'] == 'taxid'
        assert result['accid_col'] == 'assembly_accession'

    def test_detect_taxid_only(self):
        df = pd.DataFrame({'taxid': [1, 2]})
        result = detect_id_columns(df)
        assert result['taxid_col'] == 'taxid'
        assert result['accid_col'] is None

    def test_detect_none(self):
        df = pd.DataFrame({'name': ['a', 'b']})
        result = detect_id_columns(df)
        assert result['taxid_col'] is None
        assert result['accid_col'] is None


class TestRenameColumnsToStandard:
    """Tests for rename_columns_to_standard function."""

    def test_rename_taxid(self):
        df = pd.DataFrame({'TaxID': [1, 2], 'name': ['a', 'b']})
        result = rename_columns_to_standard(df, taxid_col='TaxID')
        assert 'taxid' in result.columns
        assert 'TaxID' not in result.columns

    def test_rename_both(self):
        df = pd.DataFrame({'TaxID': [1, 2], 'accession': ['A', 'B']})
        result = rename_columns_to_standard(
            df, taxid_col='TaxID', accid_col='accession'
        )
        assert 'taxid' in result.columns
        assert 'accid' in result.columns


class TestSafeGetColumn:
    """Tests for safe_get_column function."""

    def test_get_existing_column(self):
        df = pd.DataFrame({'taxid': [1, 2]})
        result = safe_get_column(df, 'taxid', default=0)
        assert result is not None

    def test_get_missing_column(self):
        df = pd.DataFrame({'name': ['a']})
        result = safe_get_column(df, 'taxid', default=0)
        assert result == 0


class TestValidateRequiredColumns:
    """Tests for validate_required_columns function."""

    def test_valid_columns(self):
        df = pd.DataFrame({'a': [1], 'b': [2]})
        validate_required_columns(df, ['a', 'b'], "Test")

    def test_missing_columns(self):
        df = pd.DataFrame({'a': [1]})
        with pytest.raises(ValueError, match="Missing required columns"):
            validate_required_columns(df, ['a', 'b'], "Test")


class TestSetupLogger:
    """Tests for setup_logger function."""

    def test_logger_creation(self):
        logger = setup_logger("test_logger")
        assert logger.name == "test_logger"

    def test_logger_format(self):
        logger = setup_logger("test_format", level=10)
        assert len(logger.handlers) > 0
