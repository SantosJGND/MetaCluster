"""
Result models with JSON schema support for evaluator module.

JSON Schema:
{
    "type": "object",
    "properties": {
        "generated_at": {"type": "string", "format": "date-time"},
        "metadata": {"type": "object"},
        "test_results": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "overall_precision": {"type": "number"},
                    "data_set": {"type": "string"}
                }
            }
        },
        "summary_results": {
            "type": "array",
            "items": {"type": "object"}
        },
        "trash_composition": {
            "type": "array",
            "items": {"type": "object"}
        },
        "cross_hit_composition": {
            "type": "array",
            "items": {"type": "object"}
        }
    }
}
"""

from dataclasses import dataclass, field, asdict
from typing import Optional, Any
import pandas as pd
import json
from datetime import datetime


@dataclass
class PrecisionMetrics:
    raw_pred_accuracy: float = 0.0
    fuzzy_precision_raw: float = 0.0
    fuzzy_precision_cov_filtered: float = 0.0
    overall_precision_raw: float = 0.0
    clade_precision_full: float = 0.0
    clade_precision_post: float = 0.0
    
    def to_dict(self) -> dict:
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict) -> 'PrecisionMetrics':
        return cls(**data)


@dataclass  
class RecallMetrics:
    recall_raw: float = 0.0
    recall_cov_filtered: float = 0.0
    clade_recall: float = 0.0
    recall_filtered_leaves: float = 0.0
    
    def to_dict(self) -> dict:
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict) -> 'RecallMetrics':
        return cls(**data)


@dataclass
class CrossHitMetrics:
    predicted_cross_hits: int = 0
    cleanup_accuracy: float = 0.0
    
    def to_dict(self) -> dict:
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: dict) -> 'CrossHitMetrics':
        return cls(**data)


@dataclass
class DatasetResult:
    data_set: str
    sample: str = ''
    input_taxid_count: int = 0
    output_raw: int = 0
    output_cov_filtered: int = 0
    predicted_clades_pre: int = 0
    predicted_clades_post: int = 0
    
    precision: PrecisionMetrics = field(default_factory=PrecisionMetrics)
    recall: RecallMetrics = field(default_factory=RecallMetrics)
    cross_hit: CrossHitMetrics = field(default_factory=CrossHitMetrics)
    
    trash_composition: Optional[dict] = None
    cross_hit_composition: Optional[dict] = None
    
    def to_dict(self) -> dict:
        result = {
            'data_set': self.data_set,
            'sample': self.sample,
            'input_taxid_count': self.input_taxid_count,
            'output_raw': self.output_raw,
            'output_cov_filtered': self.output_cov_filtered,
            'predicted_clades_pre': self.predicted_clades_pre,
            'predicted_clades_post': self.predicted_clades_post,
            'precision': self.precision.to_dict(),
            'recall': self.recall.to_dict(),
            'cross_hit': self.cross_hit.to_dict(),
        }
        if self.trash_composition is not None:
            result['trash_composition'] = self.trash_composition
        if self.cross_hit_composition is not None:
            result['cross_hit_composition'] = self.cross_hit_composition
        return result
    
    @classmethod
    def from_dict(cls, data: dict) -> 'DatasetResult':
        return cls(
            data_set=data['data_set'],
            sample=data.get('sample', ''),
            input_taxid_count=data.get('input_taxid_count', 0),
            output_raw=data.get('output_raw', 0),
            output_cov_filtered=data.get('output_cov_filtered', 0),
            predicted_clades_pre=data.get('predicted_clades_pre', 0),
            predicted_clades_post=data.get('predicted_clades_post', 0),
            precision=PrecisionMetrics.from_dict(data.get('precision', {})),
            recall=RecallMetrics.from_dict(data.get('recall', {})),
            cross_hit=CrossHitMetrics.from_dict(data.get('cross_hit', {})),
            trash_composition=data.get('trash_composition'),
            cross_hit_composition=data.get('cross_hit_composition'),
        )


@dataclass
class BatchEvaluationResult:
    test_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    summary_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    trash_composition: pd.DataFrame = field(default_factory=pd.DataFrame)
    cross_hit_composition: pd.DataFrame = field(default_factory=pd.DataFrame)
    
    metadata: dict = field(default_factory=dict)
    
    def to_json(self, filepath: str) -> None:
        """Save complete results to JSON file."""
        output = {
            'generated_at': datetime.now().isoformat(),
            'metadata': self.metadata,
            'test_results': self.test_results.to_dict(orient='records'),
            'summary_results': self.summary_results.to_dict(orient='records'),
            'trash_composition': self.trash_composition.to_dict(orient='records') if not self.trash_composition.empty else [],
            'cross_hit_composition': self.cross_hit_composition.to_dict(orient='records') if not self.cross_hit_composition.empty else [],
        }
        with open(filepath, 'w') as f:
            json.dump(output, f, indent=2, default=str)
    
    @classmethod
    def from_json(cls, filepath: str) -> 'BatchEvaluationResult':
        """Load results from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        return cls(
            test_results=pd.DataFrame(data.get('test_results', [])),
            summary_results=pd.DataFrame(data.get('summary_results', [])),
            trash_composition=pd.DataFrame(data.get('trash_composition', [])),
            cross_hit_composition=pd.DataFrame(data.get('cross_hit_composition', [])),
            metadata=data.get('metadata', {}),
        )
    
    def save_tsv(self, output_dir: str) -> None:
        """Save DataFrames to TSV files (legacy format)."""
        self.test_results.to_csv(f"{output_dir}/test_datasets_overall_precision.tsv", sep="\t", index=False)
        self.summary_results.to_csv(f"{output_dir}/test_datasets_summary_results.tsv", sep="\t", index=False)
        self.trash_composition.to_csv(f"{output_dir}/test_datasets_trash_composition.tsv", sep="\t", index=False)
        self.cross_hit_composition.to_csv(f"{output_dir}/test_datasets_cross_hit_composition.tsv", sep="\t", index=False)
    
    def get_summary_stats(self) -> dict:
        """Compute summary statistics for precision metrics."""
        precision_cols = [
            'overall_precision_raw', 'fuzzy_precision_raw', 'fuzzy_precision_cov_filtered',
            'clade_precision_full', 'clade_precision_post'
        ]
        
        available_cols = [c for c in precision_cols if c in self.summary_results.columns]
        
        if not available_cols:
            return {}
        
        stats_df = self.summary_results[available_cols].describe()
        return stats_df.to_dict()
    
    def get_dataset_count(self) -> int:
        """Get number of datasets processed."""
        return self.summary_results['data_set'].nunique() if 'data_set' in self.summary_results.columns else 0


def create_empty_result() -> BatchEvaluationResult:
    """Create an empty BatchEvaluationResult."""
    return BatchEvaluationResult(
        test_results=pd.DataFrame(columns=['overall_precision', 'data_set']),
        summary_results=pd.DataFrame(),
        trash_composition=pd.DataFrame(),
        cross_hit_composition=pd.DataFrame(),
    )
