"""
Data loader functions for evaluator module.

These functions handle loading and processing of simulation input data,
clade reports, and taxid information.
"""

import os
import logging
from typing import List, Optional, Tuple

import pandas as pd
import numpy as np

from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper


logger = logging.getLogger(__name__)


def retrieve_simulation_input(study_output_filepath: str) -> pd.DataFrame:
    """
    Retrieve simulation input data from study output directories.
    
    Args:
        study_output_filepath: Path to the study output directory containing dataset subfolders.
    
    Returns:
        DataFrame with all input data merged from all datasets.
    """
    input_files = []
    folders = [f for f in os.listdir(study_output_filepath) if os.path.isdir(os.path.join(study_output_filepath, f))]

    for data_set_name in folders:
        output_filepath = os.path.join(study_output_filepath, f"{data_set_name}", "input", f"{data_set_name}.tsv")
        if not os.path.exists(output_filepath):
            logger.debug(f"Input file not found for {data_set_name}, skipping")
            continue
            
        df = pd.read_csv(output_filepath, sep="\t")

        df['data_set'] = data_set_name.replace("_plan", "")
        if df.empty:
            continue
        input_files.append(df)

    if not input_files:
        logger.warning(f"No input files found in {study_output_filepath}")
        return pd.DataFrame()
    
    all_input_data = pd.concat(input_files, ignore_index=True)
    logger.info(f"Loaded {len(all_input_data)} input records from {len(input_files)} datasets")
    return all_input_data


def process_clade_report(data_set_name: str, study_output_filepath: str) -> List[int]:
    """
    Process clade report for a dataset and extract unique taxids.
    
    Args:
        data_set_name: Name of the dataset.
        study_output_filepath: Path to the study output directory.
    
    Returns:
        List of unique taxids found in the clade report.
    """
    clade_report = os.path.join(study_output_filepath, f"{data_set_name}", "output", "clade_report_with_references.tsv")
    if not os.path.exists(clade_report):
        logger.debug(f"Clade report not found for {data_set_name}")
        return []

    clade_df = pd.read_csv(clade_report, sep="\t").rename(columns={"#rname": "assembly_accession"})
    taxids = clade_df['taxid'].dropna().unique().tolist()
    return taxids


def output_parse(
    study_output_filepath: str, 
    ncbi_wrapper: NCBITaxonomistWrapper, 
    input_taxids: List[int]
) -> pd.DataFrame:
    """
    Parse output from study and create taxid table with lineage information.
    
    Args:
        study_output_filepath: Path to the study output directory.
        ncbi_wrapper: NCBI TaxonomistWrapper for lineage resolution.
        input_taxids: List of input taxids to include.
    
    Returns:
        DataFrame with taxid lineage information (order, family, genus).
    """
    folders = [f for f in os.listdir(study_output_filepath) if os.path.isdir(os.path.join(study_output_filepath, f))]
    all_output_taxids = set()
    for data_set_name in folders:
        output_taxids = process_clade_report(data_set_name, study_output_filepath)
        all_output_taxids.update(output_taxids)
    
    taxids = list(all_output_taxids) + input_taxids

    ncbi_wrapper.resolve_lineages(taxids)

    output_taxids_table = pd.DataFrame(list(all_output_taxids), columns=['taxid'])
    output_taxids_table['order'] = output_taxids_table.apply(
        lambda row: ncbi_wrapper.get_level(row['taxid'], 'order'), axis=1
    )
    output_taxids_table['family'] = output_taxids_table.apply(
        lambda row: ncbi_wrapper.get_level(row['taxid'], 'family'), axis=1
    )
    output_taxids_table['genus'] = output_taxids_table.apply(
        lambda row: ncbi_wrapper.get_level(row['taxid'], 'genus'), axis=1
    )

    return output_taxids_table


def establish_taxids_to_use(
    study_output_filepath: str, 
    ncbi_wrapper: NCBITaxonomistWrapper, 
    input_taxids: List[int], 
    tax_level_to_use: str = 'order', 
    min_tax_count: float = 0.02, 
    normalize: bool = True
) -> pd.DataFrame:
    """
    Establish taxids to use based on frequency threshold at a specific taxonomic level.
    
    Args:
        study_output_filepath: Path to the study output directory.
        ncbi_wrapper: NCBI TaxonomistWrapper for lineage resolution.
        input_taxids: List of input taxids.
        tax_level_to_use: Taxonomic level to use for filtering (default: 'order').
        min_tax_count: Minimum frequency threshold (default: 0.02).
        normalize: Whether to normalize counts (default: True).
    
    Returns:
        DataFrame of taxids to use with lineage information.
    """
    output_taxids_table = output_parse(study_output_filepath, ncbi_wrapper, input_taxids)
    
    if tax_level_to_use not in output_taxids_table.columns:
        logger.warning(f"Tax level '{tax_level_to_use}' not found in output taxids table")
        tax_level_value_counts = pd.Series(dtype=float)
    else:
        tax_level_value_counts = output_taxids_table[tax_level_to_use].value_counts(normalize=normalize)
    
    taxids_to_use = output_taxids_table[
        output_taxids_table[tax_level_to_use].map(tax_level_value_counts) > min_tax_count
    ]

    taxids_to_use = taxids_to_use.dropna(subset=[tax_level_to_use]).drop_duplicates(subset=[tax_level_to_use]).reset_index(drop=True)
    
    unclassified_row = pd.DataFrame({
        'taxid': [0], 
        'order': ['unclassified'], 
        'family': ['unclassified'], 
        'genus': ['unclassified']
    })
    taxids_to_use = pd.concat([taxids_to_use, unclassified_row], ignore_index=True)
    
    logger.info(f"Established {len(taxids_to_use)} taxids to use at level '{tax_level_to_use}'")
    return taxids_to_use


def expand_input_data(
    all_input_data: pd.DataFrame, 
    ncbi_wrapper: NCBITaxonomistWrapper, 
    taxid_plan: pd.DataFrame
) -> pd.DataFrame:
    """
    Expand input data with lineage information from NCBI Taxonomist.
    
    Args:
        all_input_data: Input data DataFrame
        ncbi_wrapper: NCBI TaxonomistWrapper for lineage resolution
        taxid_plan: Taxid plan with lineage information
    
    Returns:
        DataFrame with lineage information added
    """
    all_input_data = all_input_data.merge(taxid_plan, on='taxid', how='left')
    input_tax_df = pd.DataFrame(all_input_data['taxid'].dropna().unique(), columns=['taxid'])
    input_tax_df['order'] = input_tax_df.apply(
        lambda row: ncbi_wrapper.get_level(row['taxid'], 'order'), axis=1
    )
    input_tax_df['family'] = input_tax_df.apply(
        lambda row: ncbi_wrapper.get_level(row['taxid'], 'family'), axis=1
    )
    input_tax_df = input_tax_df.drop_duplicates(subset=['taxid'])

    unclassified_row = pd.DataFrame({
        'taxid': [0], 
        'order': ['unclassified'], 
        'family': ['unclassified'], 
        'genus': ['unclassified']
    })
    input_tax_df = pd.concat([input_tax_df, unclassified_row], ignore_index=True)
    input_tax_df = input_tax_df.replace({np.nan: 'unclassified'})
    
    logger.info(f"Expanded input data with {len(input_tax_df)} unique taxids")
    return input_tax_df


def load_taxid_plan(taxid_plan_filepath: str) -> pd.DataFrame:
    """
    Load taxid plan from file.
    
    Args:
        taxid_plan_filepath: Path to taxid plan file
    
    Returns:
        Cleaned taxid plan DataFrame
    """
    taxid_plan = pd.read_csv(taxid_plan_filepath, sep="\t")
    taxid_plan = taxid_plan[['taxid', 'description', 'lineage']].drop_duplicates(subset=['taxid'])
    logger.info(f"Loaded taxid plan with {len(taxid_plan)} entries")
    return taxid_plan


def get_dataset_folders(study_output_filepath: str, max_training: Optional[int] = None) -> List[str]:
    """
    Get list of dataset folders, optionally limited by max_training.
    
    Args:
        study_output_filepath: Path to study output directory
        max_training: Optional maximum number of folders to return
    
    Returns:
        List of dataset folder names
    """
    from random import sample
    
    folders = [f for f in os.listdir(study_output_filepath) if os.path.isdir(os.path.join(study_output_filepath, f))]
    
    if max_training is not None and len(folders) > max_training:
        folders = sample(folders, max_training)
        logger.info(f"Sampled {max_training} folders from {len(folders)} total")
    
    return folders


def split_train_test_folders(
    folders: List[str], 
    proportion_train: float
) -> Tuple[List[str], List[str]]:
    """
    Split folders into training and test sets.
    
    Args:
        folders: List of all folder names
        proportion_train: Proportion of data for training
    
    Returns:
        Tuple of (training_folders, test_folders)
    """
    split_idx = int(len(folders) * proportion_train)
    training_folders = folders[:split_idx]
    test_folders = folders[split_idx:]
    
    logger.info(f"Split {len(folders)} folders into {len(training_folders)} train / {len(test_folders)} test")
    return training_folders, test_folders


class DataLoader:
    """
    Data loader class for managing data retrieval.
    """
    
    def __init__(self, config):
        """
        Initialize data loader.
        
        Args:
            config: EvaluatorConfig instance
        """
        self.config = config
        self.ncbi_wrapper: Optional[NCBITaxonomistWrapper] = None
        self.all_input_data: Optional[pd.DataFrame] = None
        self.input_tax_df: Optional[pd.DataFrame] = None
        self.taxids_to_use: Optional[pd.DataFrame] = None
        self.taxid_plan: Optional[pd.DataFrame] = None
        self.folders: Optional[List[str]] = None
    
    def initialize(self) -> 'DataLoader':
        """
        Initialize all data loading.
        
        Returns:
            Self for method chaining
        """
        self._load_taxid_plan()
        self._load_datasets()
        self._resolve_lineages()
        self._expand_input_data()
        self._establish_taxids()
        return self
    
    def _load_taxid_plan(self):
        """Load taxid plan."""
        self.taxid_plan = load_taxid_plan(self.config.taxid_plan_filepath)
    
    def _load_datasets(self):
        """Load all dataset folders."""
        self.folders = get_dataset_folders(
            self.config.study_output_filepath, 
            self.config.max_training
        )
        self.all_input_data = retrieve_simulation_input(self.config.study_output_filepath)
        
        output_db = os.path.join(self.config.study_output_filepath, "taxa.db")
        self.ncbi_wrapper = NCBITaxonomistWrapper(db=output_db)
    
    def _resolve_lineages(self):
        """Resolve lineages for input taxids."""
        if self.all_input_data is None or self.ncbi_wrapper is None:
            raise RuntimeError("Data not loaded. Call _load_datasets first.")
        
        input_taxids = self.all_input_data['taxid'].dropna().unique().tolist()
        self.ncbi_wrapper.resolve_lineages(input_taxids)
        logger.info(f"Resolved {len(input_taxids)} lineages")
    
    def _expand_input_data(self):
        """Expand input data with lineage information."""
        if self.all_input_data is None or self.ncbi_wrapper is None or self.taxid_plan is None:
            raise RuntimeError("Data not loaded. Call _load_datasets and _load_taxid_plan first.")
        
        self.input_tax_df = expand_input_data(self.all_input_data, self.ncbi_wrapper, self.taxid_plan)
    
    def _establish_taxids(self):
        """Establish taxids to use."""
        if self.ncbi_wrapper is None or self.all_input_data is None:
            raise RuntimeError("Data not loaded.")
        
        input_taxids = self.all_input_data['taxid'].dropna().unique().tolist()
        self.taxids_to_use = establish_taxids_to_use(
            self.config.study_output_filepath,
            self.ncbi_wrapper,
            input_taxids,
            tax_level_to_use=self.config.tax_level,
            min_tax_count=self.config.taxa_threshold,
            normalize=True
        )
    
    def get_training_folders(self) -> List[str]:
        """Get training folder names."""
        if self.folders is None:
            raise RuntimeError("Data not loaded.")
        split_idx = int(len(self.folders) * self.config.proportion_train)
        return self.folders[:split_idx]
    
    def get_test_folders(self) -> List[str]:
        """Get test folder names."""
        if self.folders is None:
            raise RuntimeError("Data not loaded.")
        split_idx = int(len(self.folders) * self.config.proportion_train)
        return self.folders[split_idx:]
    
    def get_ncbi_wrapper(self) -> NCBITaxonomistWrapper:
        """Get NCBI wrapper."""
        if self.ncbi_wrapper is None:
            raise RuntimeError("Data not loaded.")
        return self.ncbi_wrapper
    
    def get_input_tax_df(self) -> pd.DataFrame:
        """Get input tax DataFrame."""
        if self.input_tax_df is None:
            raise RuntimeError("Data not loaded.")
        return self.input_tax_df
    
    def get_taxids_to_use(self) -> pd.DataFrame:
        """Get taxids to use."""
        if self.taxids_to_use is None:
            raise RuntimeError("Data not loaded.")
        return self.taxids_to_use
