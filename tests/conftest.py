"""
Pytest configuration and shared fixtures.
"""
import os
import tempfile
from pathlib import Path
from typing import Generator

import pandas as pd
import pytest

os.environ["NCBI_EMAIL"] = "test@example.com"


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_centrifuge_report(temp_dir: Path) -> Path:
    """Create a sample Centrifuge report file."""
    content = """name\ttaxID\trank\taxReads\tkmers\tnumUniqueReads\tdupRate\tabundances
Escherichia coli\t562\tS\t1000\t500\t100\t0.1\t50.5
Salmonella enterica\t28901\tS\t800\t400\t80\t0.1\t40.2
"""
    file_path = temp_dir / "centrifuge_report.txt"
    file_path.write_text(content)
    return file_path


@pytest.fixture
def sample_kraken_report(temp_dir: Path) -> Path:
    """Create a sample Kraken2 report file."""
    content = """PercReads\tNumReadsRoot\tNreads\tRankCode\ttaxID\tname
0.00\t0\t0\t-\t1\troot
10.50\t1000\t0\t-\t2\tcellular organisms
10.50\t1000\t0\tD\t2157\tArchaea
50.50\t500\t0\tD\t2\tBacteria
30.30\t300\t0\tP\t1224\tProteobacteria
20.20\t200\t0\tO\t9131\tEnterobacterales
15.15\t150\t0\tF\t543\tEnterobacteriaceae
10.10\t100\t0\tG\t561\tEscherichia
 5.05\t50\t0\tS\t562\tEscherichia coli
"""
    file_path = temp_dir / "kraken_report.txt"
    file_path.write_text(content)
    return file_path


@pytest.fixture
def sample_krakenunique_report(temp_dir: Path) -> Path:
    """Create a sample KrakenUnique report file."""
    content = """C\tseq1\t562\t100\t95
C\tseq2\t562\t100\t95
C\tseq3\t28901\t100\t90
U\tseq4\t0\t100\t0
"""
    file_path = temp_dir / "krakenunique_report.txt"
    file_path.write_text(content)
    return file_path


@pytest.fixture
def sample_diamond_report(temp_dir: Path) -> Path:
    """Create a sample Diamond BLAST report file."""
    content = """seq1\tNP_001234.1\t95.5\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200
seq2\tNP_001234.1\t95.5\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200
seq3\tNP_005678.1\t85.0\t100\t0\t0\t1\t100\t1\t100\t1e-40\t180
"""
    file_path = temp_dir / "diamond_report.tsv"
    file_path.write_text(content)
    return file_path


@pytest.fixture
def sample_classification_df() -> pd.DataFrame:
    """Create a sample classification DataFrame."""
    return pd.DataFrame({
        'taxid': [562, 28901, 456],
        'description': ['Escherichia coli', 'Salmonella enterica', 'Unknown'],
        'uniq_reads': [100, 80, 10]
    })


@pytest.fixture
def sample_input_table(temp_dir: Path) -> Path:
    """Create a sample input table for testing."""
    content = """taxid\tdescription\treads
562\tEscherichia coli\t1000
28901\tSalmonella enterica\t800
"""
    file_path = temp_dir / "input_table.tsv"
    file_path.write_text(content)
    return file_path


@pytest.fixture
def mock_ncbi_response():
    """Mock NCBI Entrez response structure."""
    return [
        {
            'TaxId': 562,
            'ScientificName': 'Escherichia coli',
            'Lineage': 'cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia'
        }
    ]
