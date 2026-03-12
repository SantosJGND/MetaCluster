"""
Classifier output processors for various metagenomic classifiers.
"""
from __future__ import annotations

import logging
import os
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional, Tuple, Union

import dotenv
import networkx as nx
import pandas as pd
from Bio import Entrez

dotenv.load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL", None)
if Entrez.email is None:
    raise ValueError("NCBI_EMAIL environment variable not set.")


def setup_logger(name: str) -> logging.Logger:
    """Set up a standardized logger."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.propagate = False
    return logger


logger = setup_logger(__name__)


def taxid_to_description(taxid: int) -> Optional[str]:
    """
    Given a taxid, return the corresponding scientific name using NCBI Entrez.
    
    Args:
        taxid: NCBI taxonomy ID
        
    Returns:
        Scientific name or None if not found
    """
    if taxid is None:
        return None
    if taxid >= 10**7:
        return None

    try:
        email = os.getenv("NCBI_EMAIL", None)
        if email is None:
            raise ValueError("NCBI_EMAIL environment variable not set.")
        Entrez.email = email

        handle = Entrez.esummary(db="taxonomy", id=taxid)
        record = Entrez.read(handle)
        handle.close()
        return record[0]["ScientificName"]
    except RuntimeError as e:
        logger.error(f"Error retrieving description for taxid {taxid}: {e}")
        return None
    except Exception as e:
        logger.error(f"Error retrieving description for taxid {taxid}: {e}")
        return None


def protein_accession_to_taxid(accession: str) -> Optional[int]:
    """
    Given a protein accession, return the corresponding taxid using NCBI Entrez.
    
    Args:
        accession: Protein accession number
        
    Returns:
        NCBI taxonomy ID or None if not found
    """
    try:
        Entrez.email = os.getenv("NCBI_EMAIL", "your_email@example.com")
        handle = Entrez.esearch(db="protein", term=accession)
        record = Entrez.read(handle)
        handle.close()
        
        if not record.get("IdList"):
            return None
            
        record_handle = Entrez.esummary(db="protein", id=record["IdList"])
        parsed_record = Entrez.parse(record_handle)
        
        for rec in parsed_record:
            taxid = int(rec["TaxId"])
            return taxid
        
        return None
    except Exception as e:
        logger.error(f"Error retrieving taxid for accession {accession}: {e}")
        return None


def count_prefix_spaces(s: str) -> int:
    """Count leading spaces in a string."""
    return len(s) - len(s.lstrip())


class ClassifierOutputProcessor(ABC):
    """Abstract base class for classifier output processors."""



    def __init__(self, output_path: str) -> None:
        self.output_path = output_path
        self.report: pd.DataFrame = pd.DataFrame()
        self.software_name = ""
        self.final_report: pd.DataFrame = pd.DataFrame(
            columns=["description", "taxID", "software_name"]
        )

    @classmethod
    def from_file(cls, output_path: str) -> "ClassifierOutputProcessor":
        """
        Factory method to create processor from file.
        
        Args:
            output_path: Path to classifier output file
            
        Returns:
            Processor instance with loaded data
        """
        raise NotImplementedError("Subclasses must implement from_file")

    def process(self) -> "ClassifierOutputProcessor":
        """Process the loaded report data."""
        raise NotImplementedError("Subclasses must implement process")

    def prep_final_report(self) -> "ClassifierOutputProcessor":
        """
        Prepare the final report for output.
        
        Returns:
            Self for method chaining
        """
        self.final_report = self.final_report.drop_duplicates(
            subset=["description"]
        )
        self.final_report["description"] = (
            self.final_report["description"]
            .str.replace(" ", "_")
            .str.lower()
        )

        self.final_report["software_name"] = self.software_name

        return self

    def save(self, output_path: str) -> None:
        """
        Save the final report to a file.
        
        Args:
            output_path: Path to save the report
        """
        self.final_report.to_csv(output_path, sep="\t", index=False)


class KrakenUniqOutputProcessor(ClassifierOutputProcessor):
    """Processor for KrakenUnique output files."""


    def __init__(
        self,
        output_path: str,
        min_uniq_reads: int = 10
    ) -> None:
        super().__init__(output_path)
        self.min_uniq_reads = min_uniq_reads
        self.software_name = "KrakenUniq"

    @classmethod
    def from_file(cls, output_path: str) -> "KrakenUniqOutputProcessor":
        """Load KrakenUnique report from file."""
        instance = cls(output_path)
        try:
            report = pd.read_csv(output_path, sep="\t")
            report.columns = [
                "status", "sequence_id", "taxID", "query_length", "map_length"
            ]
            instance.report = report[report["status"] == "C"]
        except (pd.errors.EmptyDataError, FileNotFoundError):
            instance.report = pd.DataFrame(columns=[
                "status", "sequence_id", "taxID", "query_length", "map_length"
            ])
        return instance

    def process(self) -> "KrakenUniqOutputProcessor":
        """Process KrakenUnique report."""
        summary = (
            self.report
            .groupby("taxID")
            .size()
            .reset_index(name="Nreads")
        )
        summary = summary[summary["Nreads"] > self.min_uniq_reads]
        summary['description'] = summary['taxID'].apply(taxid_to_description)
        summary = (
            summary
            .dropna(subset=['description'])
            .sort_values("Nreads", ascending=False)
            .reset_index(drop=True)
        )

        self.final_report = summary[
            ["description", "taxID", "Nreads"]
        ].rename(columns={"Nreads": "uniq_reads"})
        return self


class DiamondOutputProcessor(ClassifierOutputProcessor):
    """Processor for Diamond BLAST output files."""

    def __init__(
        self,
        output_path: str,
        min_uniq_reads: int = 5
    ) -> None:
        super().__init__(output_path)
        self.min_uniq_reads = min_uniq_reads
        self.software_name = "Diamond"

    @classmethod
    def from_file(cls, output_path: str) -> "DiamondOutputProcessor":
        """Load Diamond report from file."""
        instance = cls(output_path)
        try:
            report = pd.read_csv(output_path, sep="\t", header=None)
            report.columns = [
                "query_id", "subject_id", "perc_identity", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"
            ]
            instance.report = report
        except (pd.errors.EmptyDataError, FileNotFoundError):
            instance.report = pd.DataFrame(columns=[
                "query_id", "subject_id", "perc_identity", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"
            ])
        return instance

    def process(self) -> "DiamondOutputProcessor":
        """Process Diamond report."""
        summary = (
            self.report
            .groupby("subject_id")
            .size()
            .reset_index(name="Nreads")
        )
        summary = summary[summary["Nreads"] >= self.min_uniq_reads]
        summary['taxID'] = summary['subject_id'].apply(
            protein_accession_to_taxid
        )
        summary['description'] = summary['taxID'].apply(taxid_to_description)
        summary = (
            summary
            .dropna(subset=['description'])
            .sort_values("Nreads", ascending=False)
            .reset_index(drop=True)
        )
        self.final_report = summary[
            ["description", "taxID", "Nreads"]
        ].rename(columns={"Nreads": "uniq_reads"})
        return self


class CentrifugeOutputProcessor(ClassifierOutputProcessor):
    """Processor for Centrifuge output files."""

    def __init__(
        self,
        output_path: str,
        nuniq_threshold: int = 5
    ) -> None:
        super().__init__(output_path)
        self.nuniq_threshold = nuniq_threshold
        self.software_name = "Centrifuge"

    @classmethod
    def from_file(cls, output_path: str) -> "CentrifugeOutputProcessor":
        """Load Centrifuge report from file."""
        instance = cls(output_path)
        try:
            instance.report = pd.read_csv(output_path, sep="\t")
        except (pd.errors.EmptyDataError, FileNotFoundError):
            instance.report = pd.DataFrame(
                columns=["name", "taxID", "numUniqueReads"]
            )
        return instance

    def process(self) -> "CentrifugeOutputProcessor":
        """Process Centrifuge report."""
        self.report.sort_values(
            "numUniqueReads", ascending=False, inplace=True
        )
        self.report = (
            self.report
            [self.report["numUniqueReads"] >= self.nuniq_threshold]
            .reset_index(drop=True)
        )
        self.final_report = self.report[
            ["name", "taxID", "numUniqueReads"]
        ].rename(columns={
            "name": "description",
            "numUniqueReads": "uniq_reads"
        })
        return self


# Type alias for Kraken tree nodes
KrakenNode = Tuple[int, str, int, float, int, str]
KrakenEdge = Tuple[int, int]


class KrakenOutputProcessor(ClassifierOutputProcessor):
    """Processor for Kraken2 output files."""

    def __init__(
        self,
        output_path: str,
        min_uniq_reads: int = 10
    ) -> None:
        super().__init__(output_path)
        self.min_uniq_reads = min_uniq_reads
        self.nodes_dict: Dict[int, KrakenNode] = {}
        self.edges: List[KrakenEdge] = []
        self.nodes: List[KrakenNode] = []
        self.software_name = "Kraken2"

    @classmethod
    def from_file(cls, output_path: str) -> "KrakenOutputProcessor":
        """Load Kraken2 report from file."""
        instance = cls(output_path)
        try:
            kraken_report = pd.read_csv(output_path, sep="\t")
            kraken_report.columns = [
                "PercReads", "NumReadsRoot", "Nreads", "RankCode", "taxID", "name"
            ]
            kraken_report["prefix_spaces"] = kraken_report["name"].apply(
                count_prefix_spaces
            )
        except (pd.errors.EmptyDataError, FileNotFoundError):
            kraken_report = pd.DataFrame(columns=[
                "PercReads", "NumReadsRoot", "Nreads", "RankCode", "taxID", "name"
            ])

        nodes, edges = cls.kraken_report_to_tree(kraken_report)
        instance.nodes_dict = {node[0]: node for node in nodes}
        instance.edges = edges
        instance.report = kraken_report
        instance.nodes = nodes
        return instance

    def get_node_info(self, node: int) -> KrakenNode:
        """Get node information by taxid."""
        return self.nodes_dict[node]

    def process(self) -> "KrakenOutputProcessor":
        """Process Kraken2 report."""
        leaves_simple = self.get_simplified_leaves(self.nodes, self.edges)
        leaves_simple = {
            self.get_node_info(parent): [
                self.get_node_info(leaf) for leaf in leaves
            ]
            for parent, leaves in leaves_simple.items()
        }
        leaves_summary = self.summarize_leaves(leaves_simple)
        leaves_summary = leaves_summary[
            leaves_summary["Nreads"] > self.min_uniq_reads
        ]
        self.final_report = leaves_summary.rename(columns={"name": "description"})
        self.final_report = self.final_report[
            ["description", "taxID", "Nreads"]
        ].rename(columns={"Nreads": "uniq_reads"})
        return self

    @staticmethod
    def kraken_report_to_tree(
        report: pd.DataFrame
    ) -> Tuple[List[KrakenNode], List[KrakenEdge]]:
        """Convert Kraken2 report to tree structure."""
        nodes_list: List[KrakenNode] = []
        edges_dict: Dict[int, List[int]] = {}
        edges: List[KrakenEdge] = []

        for _, row in report.iterrows():
            name = row["name"].strip()
            tax_id = int(row["taxID"])
            prefix_spaces = count_prefix_spaces(row["name"])
            perc_reads = float(row["PercReads"])
            nreads = int(row["Nreads"])
            tax_rank = str(row["RankCode"])

            node = (tax_id, name, prefix_spaces, perc_reads, nreads, tax_rank)
            nodes_list.append(node)

            parent_node = None
            if prefix_spaces > 0:
                i = len(nodes_list) - 1
                while i > 0 and nodes_list[i][2] >= prefix_spaces:
                    i -= 1
                if i >= 0:
                    parent_node = nodes_list[i]

            if parent_node:
                parent_tax_id = parent_node[0]
                if parent_tax_id not in edges_dict:
                    edges_dict[parent_tax_id] = []
                edges_dict[parent_tax_id].append(tax_id)
                edges.append((parent_tax_id, tax_id))

        return nodes_list, edges

    @staticmethod
    def get_leaves(
        nodes: List[KrakenNode],
        edges: List[KrakenEdge]
    ) -> List[int]:
        """Get leaf nodes from tree."""
        G = nx.DiGraph()
        G.add_edges_from(edges)
        return [node[0] for node in nodes if G.out_degree(node[0]) == 0]

    def get_simplified_leaves(
        self,
        nodes: List[KrakenNode],
        edges: List[KrakenEdge]
    ) -> Dict[int, List[int]]:
        """Get simplified leaf structure."""
        G = nx.DiGraph()
        G.add_edges_from(edges)
        simplified_leaves: Dict[int, List[int]] = {}

        for node in nodes:
            if G.out_degree(node[0]) == 0:
                leaf = node[0]
                parent = list(G.predecessors(leaf))
                parent_detail = (
                    self.get_node_info(parent[0])
                    if parent else None
                )
                while parent and "S" in (parent_detail[5] if parent_detail else ""):
                    leaf = parent[0]
                    parent = list(G.predecessors(leaf))
                    parent_detail = (
                        self.get_node_info(parent[0])
                        if parent else None
                    )
                if parent:
                    parent = parent[0]
                    if parent not in simplified_leaves:
                        simplified_leaves[parent] = []
                    simplified_leaves[parent].append(leaf)

        return simplified_leaves

    @staticmethod
    def summarize_leaves(
        leaves_simple: Dict[KrakenNode, List[KrakenNode]]
    ) -> pd.DataFrame:
        """Summarize leaf nodes into final report."""
        summary: List[KrakenNode] = []
        for parent, leaves in leaves_simple.items():
            if "S" in parent[5]:
                summary.append(parent)
            else:
                best_leaf = max(leaves, key=lambda x: x[3])
                summary.append(best_leaf)

        return pd.DataFrame(
            summary,
            columns=[
                "taxID", "description", "prefix_spaces",
                "perc_reads", "Nreads", "rank_code"
            ]
        )
