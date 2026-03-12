"""
Shared utilities for data processing across the pipeline.
"""
import logging
from typing import Optional, List, Dict, Any

import pandas as pd


def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """
    Set up a standardized logger with consistent formatting.
    
    Args:
        name: Logger name (typically __name__)
        level: Logging level (default: INFO)
    
    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(level)
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    logger.propagate = False
    return logger


def detect_column(df: pd.DataFrame, *possible_names: str) -> Optional[str]:
    """
    Detect which column from possible_names exists in the DataFrame.
    
    Args:
        df: Input DataFrame
        *possible_names: Possible column names to look for
    
    Returns:
        Column name if found, None otherwise
    """
    for name in possible_names:
        if name in df.columns:
            return name
    return None


def detect_id_columns(df: pd.DataFrame) -> Dict[str, Optional[str]]:
    """
    Detect taxonomic ID and accession columns in a DataFrame.
    
    Supported taxid columns: 'taxid', 'TaxID', 'taxon'
    Supported accession columns: 'assembly_accession', 'accession', 'accID', 'accid'
    
    Args:
        df: Input DataFrame
    
    Returns:
        Dictionary with 'taxid_col' and 'accid_col' keys
    """
    taxid_col = detect_column(df, 'taxid', 'TaxID', 'taxon', 'taxID')
    accid_col = detect_column(
        df, 'assembly_accession', 'accession', 'accID', 'accid'
    )
    
    return {'taxid_col': taxid_col, 'accid_col': accid_col}


def rename_columns_to_standard(
    df: pd.DataFrame,
    taxid_col: Optional[str] = None,
    accid_col: Optional[str] = None
) -> pd.DataFrame:
    """
    Rename columns to standard names ('taxid', 'accid').
    
    Args:
        df: Input DataFrame
        taxid_col: Actual taxid column name (if detected)
        accid_col: Actual accession column name (if detected)
    
    Returns:
        DataFrame with standardized column names
    """
    df = df.copy()
    
    if taxid_col and taxid_col != 'taxid':
        df.rename(columns={taxid_col: 'taxid'}, inplace=True)
    
    if accid_col and accid_col != 'accid':
        df.rename(columns={accid_col: 'accid'}, inplace=True)
    
    return df


def validate_required_columns(
    df: pd.DataFrame,
    required: List[str],
    context: str = ""
) -> None:
    """
    Validate that required columns exist in DataFrame.
    
    Args:
        df: Input DataFrame
        required: List of required column names
        context: Context string for error message
    
    Raises:
        ValueError: If any required columns are missing
    """
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(
            f"{context}Missing required columns: {missing}. "
            f"Available columns: {list(df.columns)}"
        )


def safe_get_column(
    df: pd.DataFrame,
    column: str,
    default: Any = None
) -> Any:
    """
    Safely get a column value, returning default if column doesn't exist.
    
    Args:
        df: Input DataFrame
        column: Column name
        default: Default value if column not found
    
    Returns:
        Column value or default
    """
    if column in df.columns:
        return df[column]
    return default
