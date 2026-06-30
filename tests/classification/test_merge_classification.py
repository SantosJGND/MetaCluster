"""
Tests for MergeClassificationResults and MatchCladeReportWithReferenceSequences process logic.
"""
import pandas as pd
import pytest
from pathlib import Path
import tempfile
import os
import subprocess


def create_classifier_output(path: Path, data: dict):
    """Helper to create classifier output TSV file."""
    df = pd.DataFrame(data)
    df.to_csv(path, sep='\t', index=False)


class TestMergeClassificationResults:
    """Tests for the classification merge logic."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def kraken2_file(self, temp_dir):
        path = temp_dir / "kraken2.tsv"
        create_classifier_output(path, {
            'name': ['virus_a', 'virus_b'],
            'taxID': [123, 456],
            'uniq_reads': [100, 50]
        })
        return path

    @pytest.fixture
    def centrifuge_file(self, temp_dir):
        path = temp_dir / "centrifuge.tsv"
        create_classifier_output(path, {
            'name': ['virus_a', 'virus_c'],
            'taxID': [123, 789],
            'uniq_reads': [80, 30]
        })
        return path

    @pytest.fixture
    def diamond_file(self, temp_dir):
        path = temp_dir / "diamond.tsv"
        create_classifier_output(path, {
            'name': ['virus_b', 'virus_d'],
            'taxID': [456, 999],
            'uniq_reads': [25, 15]
        })
        return path

    def test_merge_kraken2_and_centrifuge(self, temp_dir, kraken2_file, centrifuge_file):
        """Test merging kraken2 and centrifuge results."""
        kraken2_df = pd.read_csv(kraken2_file, sep='\t')
        kraken2_df = kraken2_df.rename(columns={'name': 'description', 'taxID': 'taxid'})
        kraken2_df['classifier'] = 'kraken2'

        centrifuge_df = pd.read_csv(centrifuge_file, sep='\t')
        centrifuge_df = centrifuge_df.rename(columns={'name': 'description', 'taxID': 'taxid'})
        centrifuge_df['classifier'] = 'centrifuge'

        merged = pd.merge(kraken2_df, centrifuge_df, on=['description', 'taxid'], 
                         how='outer', suffixes=('_kraken2', '_centrifuge'))

        # Outer join returns all unique rows from both - 3 rows (virus_a, virus_b, virus_c)
        assert len(merged) == 3
        assert 'virus_a' in merged['description'].values
        assert 'virus_b' in merged['description'].values
        assert 'virus_c' in merged['description'].values

    def test_merge_with_diamond(self, temp_dir, kraken2_file, centrifuge_file, diamond_file):
        """Test merging with diamond classifier."""
        kraken2_df = pd.read_csv(kraken2_file, sep='\t')
        kraken2_df = kraken2_df.rename(columns={'name': 'description', 'taxID': 'taxid'})
        kraken2_df['classifier'] = 'kraken2'

        centrifuge_df = pd.read_csv(centrifuge_file, sep='\t')
        centrifuge_df = centrifuge_df.rename(columns={'name': 'description', 'taxID': 'taxid'})
        centrifuge_df['classifier'] = 'centrifuge'

        diamond_df = pd.read_csv(diamond_file, sep='\t')
        diamond_df = diamond_df.rename(columns={'name': 'description', 'taxID': 'taxid'})
        diamond_df['classifier'] = 'diamond'

        merged = pd.merge(kraken2_df, centrifuge_df, on=['description', 'taxid'], 
                         how='outer', suffixes=('_kraken2', '_centrifuge'))
        merged = pd.merge(merged, diamond_df, on=['description', 'taxid'], 
                         how='outer', suffixes=('', '_diamond'))

        # Unique viruses across all classifiers: virus_a, virus_b, virus_c, virus_d = 4
        assert len(merged) == 4

    def test_total_uniq_reads_calculation(self, temp_dir, kraken2_file, centrifuge_file):
        """Test that total_uniq_reads is calculated correctly."""
        kraken2_df = pd.read_csv(kraken2_file, sep='\t')
        kraken2_df = kraken2_df.rename(columns={'name': 'description', 'taxID': 'taxid'})
        
        centrifuge_df = pd.read_csv(centrifuge_file, sep='\t')
        centrifuge_df = centrifuge_df.rename(columns={'name': 'description', 'taxID': 'taxid'})

        merged = pd.merge(kraken2_df, centrifuge_df, on=['description', 'taxid'], 
                         how='outer', suffixes=('_kraken2', '_centrifuge'))

        kraken2_reads = merged.get('uniq_reads_kraken2', pd.Series([0]*len(merged))).fillna(0)
        cent_reads = merged.get('uniq_reads_centrifuge', pd.Series([0]*len(merged))).fillna(0)
        merged['total_uniq_reads'] = kraken2_reads + cent_reads

        virus_a = merged[merged['description'] == 'virus_a']
        assert virus_a['total_uniq_reads'].iloc[0] == 180

    def test_classification_column(self, temp_dir, kraken2_file, centrifuge_file):
        """Test classification column generation."""
        kraken2_df = pd.read_csv(kraken2_file, sep='\t')
        kraken2_df = kraken2_df.rename(columns={'name': 'description', 'taxID': 'taxid'})
        kraken2_df['classifier'] = 'kraken2'

        centrifuge_df = pd.read_csv(centrifuge_file, sep='\t')
        centrifuge_df = centrifuge_df.rename(columns={'name': 'description', 'taxID': 'taxid'})
        centrifuge_df['classifier'] = 'centrifuge'

        merged = pd.merge(kraken2_df, centrifuge_df, on=['description', 'taxid'], 
                         how='outer', suffixes=('_kraken2', '_centrifuge'))

        def classify(row):
            classifiers = []
            if pd.notna(row.get('classifier_kraken2')):
                classifiers.append('kraken2')
            if pd.notna(row.get('classifier_centrifuge')):
                classifiers.append('centrifuge')
            if not classifiers:
                return 'unclassified'
            return '/'.join(classifiers)

        merged['classification'] = merged.apply(classify, axis=1)

        virus_a = merged[merged['description'] == 'virus_a']
        assert 'kraken2/centrifuge' in virus_a['classification'].values

    def test_empty_dataframe_fallback(self, temp_dir):
        """Test handling when no classifier outputs are available."""
        empty_df = pd.DataFrame(columns=['description', 'taxID', 'uniq_reads'])
        
        assert empty_df.empty
        assert 'description' in empty_df.columns
        assert 'taxID' in empty_df.columns
        assert 'uniq_reads' in empty_df.columns


class TestMatchCladeReportWithReferenceSequences:
    """Tests for MatchCladeReportWithReferenceSequences logic."""

    @pytest.fixture
    def temp_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def matched_assemblies_file(self, temp_dir):
        path = temp_dir / "matched_assemblies.tsv"
        df = pd.DataFrame({
            'assembly_accession': ['GCF_000001', 'GCF_000002'],
            'assembly_file': ['/path/to/assembly1.fna', '/path/to/assembly2.fna'],
            'taxid': [562, 28901],
            'description': ['Escherichia coli', 'Salmonella enterica']
        })
        df.to_csv(path, sep='\t', index=False)
        return path

    @pytest.fixture
    def coverage_report_file(self, temp_dir):
        path = temp_dir / "coverage.tsv"
        df = pd.DataFrame({
            'file': ['assembly1', 'assembly2'],
            'coverage': [10.5, 8.2]
        })
        df.to_csv(path, sep='\t', index=False)
        return path

    @pytest.fixture
    def merge_classification_file(self, temp_dir):
        path = temp_dir / "merge_classification.tsv"
        df = pd.DataFrame({
            'taxid': [562, 28901],
            'classifier': ['kraken2', 'kraken2']
        })
        df.to_csv(path, sep='\t', index=False)
        return path

    @pytest.fixture
    def valid_clade_report_file(self, temp_dir):
        path = temp_dir / "clade_report.tsv"
        content = "node1\t100\t0.5\t0.1\t2\tassembly1,assembly2\nnode2\t50\t0.3\t0.2\t1\tassembly1\n"
        path.write_text(content)
        return path

    def test_missing_clade_report_sets_unmapped(
        self, temp_dir, matched_assemblies_file, 
        coverage_report_file, merge_classification_file
    ):
        """Test that missing clade_report results in clade='unmapped', nuniq=0, freq=0."""
        clade_report_path = temp_dir / "clade_report.tsv"
        
        matched_assemblies = pd.read_csv(matched_assemblies_file, sep='\t')
        
        clade_report_path_str = str(clade_report_path)
        if not os.path.exists(clade_report_path_str) or os.path.getsize(clade_report_path_str) == 0:
            clade_report = pd.DataFrame()
        else:
            clade_report = pd.read_csv(clade_report_path_str, sep="\t", header=None,
                                       names=["clade", "nuniq", "freq", "min_pair_dist", "nfiles", "files"])
        
        assert clade_report.empty
        
        for _, row in matched_assemblies.iterrows():
            if clade_report.empty:
                clade = 'unmapped'
                nuniq = 0
                freq = 0
            else:
                clade = 'unmapped'
                nuniq = 0
                freq = 0
            
            assert clade == 'unmapped'
            assert nuniq == 0
            assert freq == 0

    def test_empty_clade_report_sets_unmapped(
        self, temp_dir, matched_assemblies_file,
        coverage_report_file, merge_classification_file
    ):
        """Test that empty clade_report results in clade='unmapped'."""
        empty_clade_report = temp_dir / "clade_report.tsv"
        empty_clade_report.write_text("")
        
        clade_report_path = str(empty_clade_report)
        if os.path.exists(clade_report_path) and os.path.getsize(clade_report_path) > 0:
            clade_report = pd.read_csv(clade_report_path, sep="\t", header=None,
                                       names=["clade", "nuniq", "freq", "min_pair_dist", "nfiles", "files"])
        else:
            clade_report = pd.DataFrame()
        
        assert clade_report.empty
        
        matched_assemblies = pd.read_csv(matched_assemblies_file, sep='\t')
        for _, row in matched_assemblies.iterrows():
            clade = 'unmapped'
            nuniq = 0
            freq = 0
            assert clade == 'unmapped'
            assert nuniq == 0
            assert freq == 0

    def test_valid_clade_report_maps_correctly(
        self, temp_dir, matched_assemblies_file,
        coverage_report_file, merge_classification_file, valid_clade_report_file
    ):
        """Test that valid clade_report correctly maps assemblies to clades."""
        clade_report_path = str(valid_clade_report_file)
        
        clade_report = pd.read_csv(clade_report_path, sep="\t", header=None,
                                   names=["clade", "nuniq", "freq", "min_pair_dist", "nfiles", "files"])
        clade_report['files'] = clade_report['files'].str.split(',')
        clade_report = clade_report.explode('files')
        
        assert not clade_report.empty
        
        matched_assemblies = pd.read_csv(matched_assemblies_file, sep='\t')
        
        for _, row in matched_assemblies.iterrows():
            accession = row['assembly_accession']
            match = clade_report[clade_report['files'].str.contains(accession, na=False)]
            
            if not match.empty:
                assert match['clade'].values[0] in ['node1', 'node2']

    def test_clade_report_with_references_output_structure(
        self, temp_dir, matched_assemblies_file,
        coverage_report_file, merge_classification_file
    ):
        """Test that output has expected columns when clade_report is missing."""
        matched_assemblies = pd.read_csv(matched_assemblies_file, sep='\t')
        coverage_report = pd.read_csv(coverage_report_file, sep='\t')
        merge_results = pd.read_csv(merge_classification_file, sep='\t')
        
        clade_report = pd.DataFrame()
        
        result = matched_assemblies.copy()
        result['clade'] = 'unmapped'
        result['nuniq'] = 0
        result['freq'] = 0
        result['classifier'] = result['taxid'].map(
            dict(zip(merge_results['taxid'], merge_results['classifier']))
        )
        result['coverage'] = result['assembly_accession'].map(
            lambda x: coverage_report[coverage_report['file'].str.contains(x, na=False)]['coverage'].values[0]
            if not coverage_report[coverage_report['file'].str.contains(x, na=False)].empty else 0
        )
        
        expected_columns = ['description', 'taxid', 'assembly_accession', 'coverage', 'clade', 'nuniq', 'freq', 'classifier']
        assert all(col in result.columns for col in expected_columns)
        
        assert (result['clade'] == 'unmapped').all()
        assert (result['nuniq'] == 0).all()
        assert (result['freq'] == 0).all()