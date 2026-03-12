"""
Visualization module for evaluator results.

Generates plots and charts from evaluation results.
"""

import os
from typing import Optional

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from .result_models import BatchEvaluationResult


class ResultVisualizer:
    """
    Generates all plots from BatchEvaluationResult.
    """
    
    def __init__(self, output_dir: str):
        """
        Initialize visualizer.
        
        Args:
            output_dir: Directory to save plots
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
    
    def plot_all(self, results: BatchEvaluationResult) -> None:
        """
        Generate all plots.
        
        Args:
            results: BatchEvaluationResult with evaluation data
        """
        self.plot_precision_distribution(results.test_results)
        self.plot_precision_comparison(results.summary_results)
        self.plot_recall_comparison(results.summary_results)
        self.plot_cross_hit_heatmap(results.cross_hit_composition)
        self.plot_trash_heatmap(results.trash_composition)
        self.plot_recall_improvement(results.summary_results)
        self.plot_probability_metrics(results.summary_results)
    
    def plot_precision_distribution(self, test_results: pd.DataFrame) -> None:
        """
        Plot histogram of overall precision distribution.
        
        Args:
            test_results: Test results DataFrame
        """
        if test_results.empty or 'overall_precision' not in test_results.columns:
            return
        
        plt.figure(figsize=(10, 6))
        sns.histplot(test_results['overall_precision'], bins=20, kde=True)
        plt.xlabel('Overall Precision')
        plt.xlim(0, 2)
        plt.ylabel('Frequency')
        plt.title('Distribution of Overall Precision Across Test Datasets')
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, "overall_precision_histogram.png"))
        plt.close()
    
    def plot_precision_comparison(self, summary_results: pd.DataFrame) -> None:
        """
        Plot boxplot comparing precision metrics.
        
        Args:
            summary_results: Summary results DataFrame
        """
        precision_cols = [
            'overall_precision_raw', 'fuzzy_precision_raw', 'fuzzy_precision_cov_filtered',
            'clade_precision_full', 'clade_precision_post'
        ]
        
        available_cols = [c for c in precision_cols if c in summary_results.columns]
        if not available_cols:
            return
        
        if 'sample' not in summary_results.columns:
            return
            
        melted = summary_results.melt(
            id_vars=['sample'], 
            value_vars=available_cols, 
            var_name='Metric', 
            value_name='Value'
        )
        
        plt.figure(figsize=(12, 8))
        sns.boxplot(x='Metric', y='Value', data=melted)
        plt.title('Comparison of Precision Metrics Across Datasets')
        plt.ylabel('Precision')
        plt.xlabel('Metric')
        plt.ylim(0, 3)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, "precision_metrics_boxplot.png"))
        plt.close()
    
    def plot_recall_comparison(self, summary_results: pd.DataFrame) -> None:
        """
        Plot boxplot comparing recall metrics.
        
        Args:
            summary_results: Summary results DataFrame
        """
        recall_cols = ['recall_raw', 'recall_cov_filtered', 'clade_recall', 'recall_filtered_leaves']
        
        available_cols = [c for c in recall_cols if c in summary_results.columns]
        if not available_cols:
            return
        
        if 'sample' not in summary_results.columns:
            return
            
        melted = summary_results.melt(
            id_vars=['sample'], 
            value_vars=available_cols, 
            var_name='Metric', 
            value_name='Value'
        )
        
        plt.figure(figsize=(12, 6))
        sns.boxplot(x='Metric', y='Value', data=melted)
        plt.title('Comparison of Recall Metrics Across Datasets')
        plt.ylabel('Recall')
        plt.xlabel('Metric')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, "recall_metrics_boxplot.png"))
        plt.close()
    
    def plot_cross_hit_heatmap(self, cross_hit_composition: pd.DataFrame) -> None:
        """
        Plot heatmap of cross-hit composition.
        
        Args:
            cross_hit_composition: Cross-hit composition DataFrame
        """
        if cross_hit_composition.empty or 'tax_level' not in cross_hit_composition.columns:
            return
        
        summary_df = self._composition_summary(cross_hit_composition)
        if summary_df.empty:
            return
        
        summary_df = summary_df.sort_index()
        summary_df = summary_df.reindex(sorted(summary_df.columns), axis=1)
        
        plt.figure(figsize=(12, 8))
        sns.heatmap(summary_df, cmap='viridis', annot=True, fmt=".2f")
        plt.title('Average Cross-Hit Composition by Tax Level')
        plt.xlabel('Taxa')
        plt.ylabel('Tax Level')
        plt.savefig(os.path.join(self.output_dir, "cross_hit_composition_heatmap.png"))
        plt.close()
    
    def plot_trash_heatmap(self, trash_composition: pd.DataFrame) -> None:
        """
        Plot heatmap of trash composition.
        
        Args:
            trash_composition: Trash composition DataFrame
        """
        if trash_composition.empty or 'tax_level' not in trash_composition.columns:
            return
        
        summary_df = self._composition_summary(trash_composition)
        if summary_df.empty:
            return
        
        summary_df = summary_df.sort_index()
        summary_df = summary_df.reindex(sorted(summary_df.columns), axis=1)
        
        plt.figure(figsize=(12, 8))
        sns.heatmap(summary_df, cmap='viridis', annot=True, fmt=".4f")
        plt.title('Average Trash Composition by Tax Level')
        plt.xlabel('Taxa')
        plt.ylabel('Tax Level')
        plt.savefig(os.path.join(self.output_dir, "trash_composition_heatmap.png"))
        plt.close()
    
    def plot_recall_improvement(self, summary_results: pd.DataFrame) -> None:
        """
        Plot histogram of recall improvement.
        
        Args:
            summary_results: Summary results DataFrame
        """
        required_cols = ['recall_raw', 'clade_recall', 'recall_filtered_leaves']
        if not all(c in summary_results.columns for c in required_cols):
            return
        
        if 'sample' not in summary_results.columns:
            return
        
        recall_df = summary_results[['sample', 'recall_raw', 'clade_recall', 'recall_filtered_leaves']].copy()
        recall_df['recall_clade_diff'] = recall_df['clade_recall'] - recall_df['recall_raw']
        recall_df['recall_clade_diff_predicted_leaves'] = recall_df['recall_filtered_leaves'] - recall_df['recall_raw']
        
        plt.figure(figsize=(12, 6))
        sns.histplot(recall_df['recall_clade_diff'], bins=20, kde=True)
        sns.histplot(recall_df['recall_clade_diff_predicted_leaves'], bins=20, kde=True, color='orange')
        plt.title('Distribution of Recall Improvement (Clade Recall - Raw Recall)')
        plt.xlabel('Recall Improvement')
        plt.ylabel('Frequency')
        plt.axvline(0, color='red', linestyle='--', label='No Improvement')
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, "recall_improvement_histogram.png"))
        plt.close()
    
    def plot_probability_metrics(self, summary_results: pd.DataFrame) -> None:
        """
        Plot probability metrics boxplot.
        
        Args:
            summary_results: Summary results DataFrame
        """
        required_cols = ['recall_raw', 'fuzzy_precision_raw', 'overall_precision_raw', 'clade_recall', 'clade_precision_full']
        if not all(c in summary_results.columns for c in required_cols):
            return
        
        if 'sample' not in summary_results.columns:
            return
        
        precisions_df = summary_results[['sample', 'recall_raw', 'fuzzy_precision_raw', 
                                         'overall_precision_raw', 'clade_recall', 'clade_precision_full']].copy()
        precisions_df = precisions_df.drop_duplicates()
        
        precisions_df['Prob_Find_any'] = precisions_df['recall_raw'] * precisions_df['fuzzy_precision_raw']
        precisions_df['Prob_Find_true'] = precisions_df['recall_raw'] * precisions_df['overall_precision_raw']
        precisions_df['Prob_Find_true_clade_full'] = precisions_df['clade_recall'] * precisions_df['clade_precision_full']
        
        prob_cols = ['Prob_Find_any', 'Prob_Find_true', 'Prob_Find_true_clade_full']
        
        melted = precisions_df.melt(
            id_vars=['sample'], 
            value_vars=prob_cols, 
            var_name='Metric', 
            value_name='Value'
        )
        
        plt.figure(figsize=(12, 8))
        sns.boxplot(x='Metric', y='Value', data=melted)
        plt.title('Comparison of Probability Metrics Across Datasets')
        plt.ylabel('Probability')
        plt.xlabel('Metric')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, "probability_metrics_boxplot.png"))
        plt.close()
    
    def _composition_summary(self, composition_df: pd.DataFrame) -> pd.DataFrame:
        """
        Generate composition summary.
        
        Args:
            composition_df: Composition DataFrame
            
        Returns:
            Summary DataFrame
        """
        if 'tax_level' not in composition_df.columns:
            return pd.DataFrame()
        
        summary_list = []
        for tax_level in composition_df['tax_level'].unique():
            subset = composition_df[composition_df['tax_level'] == tax_level]
            mean_values = subset.drop(columns=['taxid', 'tax_level', 'data_set'], errors='ignore').mean()
            mean_values = pd.DataFrame(mean_values).T
            mean_values.insert(0, 'tax_level', tax_level)
            summary_list.append(mean_values)
        
        if not summary_list:
            return pd.DataFrame()
        
        summary_list = pd.concat(summary_list, ignore_index=True)
        summary_list.set_index('tax_level', inplace=True)
        return summary_list
    
    def plot_clade_precision_by_taxlevel(self, summary_results: pd.DataFrame, tax_level: str) -> None:
        """
        Plot clade precision by taxonomic level.
        
        Args:
            summary_results: Summary results DataFrame
            tax_level: Taxonomic level column name
        """
        if tax_level not in summary_results.columns or 'raw_pred_accuracy' not in summary_results.columns:
            return
        
        plt.figure(figsize=(10, 8))
        plt.ylabel('Raw Precision (Post)')
        sns.boxplot(x=tax_level, y='raw_pred_accuracy', data=summary_results)
        plt.xticks(rotation=45)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, "clade_precision.png"))
        plt.close()

    def plot_mutation_rate_vs_crosshits(
        self,
        input_data: pd.DataFrame,
        cross_hit_data: pd.DataFrame,
    ) -> None:
        """
        Plot mutation rate vs cross-hit counts.
        
        Args:
            input_data: Input table DataFrame with mutation_rate column
            cross_hit_data: Cross-hit predictions DataFrame
        """
        if input_data.empty or 'mutation_rate' not in input_data.columns:
            return
        
        merged = input_data.merge(
            cross_hit_data,
            left_on='taxid',
            right_on='taxid',
            how='left',
        )
        
        cross_hit_counts = merged.groupby('mutation_rate').size().reset_index(name='cross_hit_count')
        
        plt.figure(figsize=(10, 6))
        sns.scatterplot(
            data=cross_hit_counts,
            x='mutation_rate',
            y='cross_hit_count',
            size='cross_hit_count',
            sizes=(50, 400),
            color='steelblue',
            alpha=0.7,
        )
        plt.xlabel('Mutation Rate')
        plt.ylabel('Cross-Hit Count')
        plt.title('Mutation Rate vs Cross-Hit Counts')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, "mutation_rate_vs_crosshits.png"))
        plt.close()

    def plot_cluster_size_distribution(
        self,
        matched_assemblies: pd.DataFrame,
        input_data: pd.DataFrame,
    ) -> None:
        """
        Plot cluster size distribution for matched taxids from input.
        
        Args:
            matched_assemblies: Matched assemblies DataFrame with cluster info
            input_data: Input table DataFrame with taxid column
        """
        if matched_assemblies.empty or input_data.empty:
            return
        
        input_taxids = set(input_data['taxid'].unique())
        
        matched_filtered = matched_assemblies[
            matched_assemblies['taxid'].isin(input_taxids)
        ]
        
        if matched_filtered.empty:
            return
        
        cluster_sizes = matched_filtered.groupby('taxid').size().reset_index(name='cluster_size')
        
        plt.figure(figsize=(12, 6))
        
        plt.subplot(1, 2, 1)
        sns.histplot(data=cluster_sizes, x='cluster_size', bins=30, kde=True, color='coral')
        plt.xlabel('Cluster Size (leaves per taxid)')
        plt.ylabel('Frequency')
        plt.title('Distribution of Cluster Sizes')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        top_clusters = cluster_sizes.nlargest(20, 'cluster_size')
        sns.barplot(data=top_clusters, x='cluster_size', y='taxid', palette='viridis')
        plt.xlabel('Cluster Size')
        plt.ylabel('Taxid')
        plt.title('Top 20 Largest Clusters')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, "cluster_size_distribution.png"))
        plt.close()

    def generate_html_report(self, results: BatchEvaluationResult) -> str:
        """
        Generate HTML report embedding all plots.
        
        Args:
            results: BatchEvaluationResult with evaluation data
        
        Returns:
            Path to generated HTML file
        """
        from datetime import datetime
        
        self.plot_all(results)
        
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        stats = results.get_summary_stats() if hasattr(results, 'get_summary_stats') else {}
        
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Model Evaluation Results</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        h1 {{
            color: #333;
            border-bottom: 2px solid #333;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #555;
            margin-top: 30px;
        }}
        .metric {{
            display: inline-block;
            margin: 10px 20px;
            padding: 15px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .metric-label {{
            font-size: 12px;
            color: #888;
        }}
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
            color: #333;
        }}
        .plot {{
            margin: 20px 0;
            padding: 15px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .plot img {{
            max-width: 100%;
            height: auto;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            color: #888;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <h1>Model Evaluation Results</h1>
    <p>Generated: {timestamp}</p>
    
    <h2>Summary Metrics</h2>
"""
        
        if stats:
            html_content += '<div class="metrics">\n'
            for metric, values in stats.items():
                if isinstance(values, dict):
                    for stat_name, value in values.items():
                        if isinstance(value, (int, float)):
                            html_content += f"""    <div class="metric">
        <div class="metric-label">{metric} ({stat_name})</div>
        <div class="metric-value">{value:.4f}</div>
    </div>\n"""
            html_content += '</div>\n'
        
        plot_files = [
            "overall_precision_histogram.png",
            "precision_metrics_boxplot.png", 
            "recall_metrics_boxplot.png",
            "recall_improvement_histogram.png",
            "probability_metrics_boxplot.png",
            "cross_hit_composition_heatmap.png",
            "trash_composition_heatmap.png",
            "mutation_rate_vs_crosshits.png",
            "cluster_size_distribution.png",
        ]
        
        html_content += """
    <h2>Visualizations</h2>
"""
        for plot_file in plot_files:
            plot_path = os.path.join(self.output_dir, plot_file)
            if os.path.exists(plot_path):
                plot_title = plot_file.replace(".png", "").replace("_", " ").title()
                html_content += f"""
    <div class="plot">
        <h3>{plot_title}</h3>
        <img src="{plot_file}" alt="{plot_title}"/>
    </div>
"""
        
        html_content += f"""
    <div class="footer">
        <p>Metagenomics Evaluation Pipeline</p>
    </div>
</body>
</html>
"""
        
        output_path = os.path.join(self.output_dir, "evaluation_report.html")
        with open(output_path, 'w') as f:
            f.write(html_content)
        
        return output_path


try:
    import plotly.express as px
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False


class PlotlyVisualizer:
    """
    Interactive visualization using Plotly.
    
    Requires plotly to be installed.
    """
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
    
    def plot_precision_distribution(self, test_results: pd.DataFrame) -> go.Figure:
        """Create interactive precision distribution histogram."""
        if not PLOTLY_AVAILABLE or test_results.empty:
            return None
        
        fig = px.histogram(
            test_results, 
            x='overall_precision',
            nbins=20,
            title='Distribution of Overall Precision',
            labels={'overall_precision': 'Precision'}
        )
        fig.update_layout(bargap=0.1)
        return fig
    
    def plot_recall_comparison(self, summary: pd.DataFrame) -> go.Figure:
        """Create interactive recall comparison box plot."""
        if not PLOTLY_AVAILABLE or summary.empty:
            return None
        
        recall_cols = ['recall_raw', 'recall_cov_filtered', 'clade_recall']
        available = [c for c in recall_cols if c in summary.columns]
        
        if not available:
            return None
        
        melted = summary.melt(id_vars=['sample'], value_vars=available, 
                            var_name='Metric', value_name='Value')
        
        fig = px.box(melted, x='Metric', y='Value', title='Recall Metrics Comparison')
        return fig
    
    def plot_cross_hit_heatmap(self, composition: pd.DataFrame) -> go.Figure:
        """Create interactive cross-hit composition heatmap."""
        if not PLOTLY_AVAILABLE or composition.empty:
            return None
        
        fig = px.imshow(
            composition.select_dtypes(include='number').T,
            title='Cross-Hit Composition',
            labels=dict(x='Taxa', y='Tax Level', color='Value')
        )
        return fig
    
    def save_interactive(self, results: BatchEvaluationResult) -> None:
        """Save all plots as interactive HTML."""
        self.plot_precision_distribution(results.test_results)
        self.plot_recall_comparison(results.summary_results)
        self.plot_cross_hit_heatmap(results.cross_hit_composition)


class ReportGenerator:
    """
    Generates HTML reports from evaluation results.
    
    Combines static plots with evaluation metrics into a
    comprehensive HTML report.
    """
    
    def __init__(self, output_dir: str, template: Optional[str] = None):
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.template = template or self._default_template()
        self.sections: list[dict] = []
    
    def add_plot(self, title: str, plot_path: str, description: str = "") -> 'ReportGenerator':
        """Add a plot section to the report."""
        self.sections.append({
            'type': 'plot',
            'title': title,
            'path': plot_path,
            'description': description,
        })
        return self
    
    def add_metrics(self, title: str, metrics: dict, description: str = "") -> 'ReportGenerator':
        """Add a metrics section to the report."""
        self.sections.append({
            'type': 'metrics',
            'title': title,
            'metrics': metrics,
            'description': description,
        })
        return self
    
    def add_text(self, title: str, content: str) -> 'ReportGenerator':
        """Add a text section to the report."""
        self.sections.append({
            'type': 'text',
            'title': title,
            'content': content,
        })
        return self
    
    def generate(self, output_path: str) -> None:
        """Generate HTML report."""
        html = self.template.format(
            title="Evaluation Report",
            generated_at=pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
            sections=self._render_sections(),
        )
        
        with open(output_path, 'w') as f:
            f.write(html)
    
    def _render_sections(self) -> str:
        """Render all sections as HTML."""
        html_parts = []
        for section in self.sections:
            if section['type'] == 'plot':
                html_parts.append(self._render_plot(section))
            elif section['type'] == 'metrics':
                html_parts.append(self._render_metrics(section))
            elif section['type'] == 'text':
                html_parts.append(self._render_text(section))
        return '\n'.join(html_parts)
    
    def _render_plot(self, section: dict) -> str:
        rel_path = os.path.relpath(section['path'], self.output_dir)
        return f"""
        <div class="section">
            <h2>{section['title']}</h2>
            <p>{section.get('description', '')}</p>
            <img src="{rel_path}" alt="{section['title']}" />
        </div>
        """
    
    def _render_metrics(self, section: dict) -> str:
        rows = []
        for key, value in section['metrics'].items():
            rows.append(f"<tr><td>{key}</td><td>{value}</td></tr>")
        
        return f"""
        <div class="section">
            <h2>{section['title']}</h2>
            <p>{section.get('description', '')}</p>
            <table>
                <thead><tr><th>Metric</th><th>Value</th></tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>
        """
    
    def _render_text(self, section: dict) -> str:
        return f"""
        <div class="section">
            <h2>{section['title']}</h2>
            <p>{section['content']}</p>
        </div>
        """
    
    def _default_template(self) -> str:
        return """<!DOCTYPE html>
<html>
<head>
    <title>{title}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #333; }}
        h2 {{ color: #666; margin-top: 30px; }}
        .section {{ margin-bottom: 40px; }}
        img {{ max-width: 100%; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <h1>{title}</h1>
    <p>Generated: {generated_at}</p>
    {sections}
</body>
</html>"""


def generate_report(
    results: BatchEvaluationResult,
    output_dir: str,
    plot_dir: Optional[str] = None,
) -> str:
    """
    Convenience function to generate a complete evaluation report.
    
    Args:
        results: Evaluation results
        output_dir: Output directory for report
        plot_dir: Directory containing plots (defaults to output_dir)
    
    Returns:
        Path to generated HTML report
    """
    if plot_dir is None:
        plot_dir = output_dir
    
    visualizer = ResultVisualizer(plot_dir)
    visualizer.plot_all(results)
    
    report_gen = ReportGenerator(plot_dir)
    
    stats = results.get_summary_stats()
    if stats:
        flat_stats = {}
        for metric, values in stats.items():
            for stat_name, value in values.items():
                flat_stats[f"{metric}_{stat_name}"] = value
        
        report_gen.add_metrics("Summary Statistics", flat_stats)
    
    report_gen.add_text("Overview", 
        f"Evaluated {results.get_dataset_count()} datasets. "
        f"Success: {results.metadata.get('successful', 0)}, "
        f"Failed: {results.metadata.get('failed', 0)}")
    
    output_path = os.path.join(output_dir, "evaluation_report.html")
    report_gen.generate(output_path)
    
    return output_path
