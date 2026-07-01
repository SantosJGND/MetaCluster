"""
Metagenomics utilities shared package.

This package contains shared utilities for:
- NCBI tools (ncbi_tools.py)
- Reference management (reference_utils.py)
- Dataframe utilities (dataframe_utils.py)
- Overlap manager for clustering analysis (overlap_manager/)
"""

__version__ = "1.0.0"

from metagenomics_utils.ncbi_tools import (
    LocalAssembly,
    NCBITaxonomistWrapper,
    NCBITools,
    Passport,
    ReferenceData,
    compare_lineages,
    retry_with_backoff,
)
from metagenomics_utils.overlap_manager.manager import OverlapManager
