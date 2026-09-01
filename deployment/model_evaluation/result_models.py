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
        "spurious_composition": {
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

import json
import logging
import os
from dataclasses import asdict, dataclass, field
from datetime import datetime

import pandas as pd

logger = logging.getLogger(__name__)


def write_pipeline_metadata(rows, output_dir: str, filename: str = "pipeline_metadata.tsv") -> str:
    """Write pipeline/cohort metadata as a two-column TSV, atomically.

    Args:
        rows: Iterable of (metric, value) tuples.
        output_dir: Directory in which to write the file.
        filename: Output filename (default: pipeline_metadata.tsv).

    Returns:
        Path to the written file.
    """
    os.makedirs(output_dir, exist_ok=True)
    df = pd.DataFrame(rows, columns=["metric", "value"])
    path = os.path.join(output_dir, filename)
    tmp_path = f"{path}.tmp.{os.getpid()}"
    df.to_csv(tmp_path, sep="\t", index=False)
    os.replace(tmp_path, path)
    logger.info(f"Saved pipeline metadata to {path}")
    return path


@dataclass
class PrecisionMetrics:
    purity_raw: float = 0.0
    purity_cov_filtered: float = 0.0
    precision_best_match: float = 0.0
    clade_precision_full: float = 0.0
    clade_precision_post: float = 0.0
    clade_precision_fixed: float = 0.0

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "PrecisionMetrics":
        return cls(**data)


@dataclass
class RecallMetrics:
    recall_raw: float = 0.0
    recall_cov_filtered: float = 0.0
    clade_recall_pre_cleanup: float = 0.0
    clade_recall_post_cleanup: float = 0.0
    recall_filtered_leaves: float = 0.0
    recall_fixed_filter: float = 0.0
    recall_metrics: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "RecallMetrics":
        return cls(**data)


@dataclass
class CrossHitMetrics:
    predicted_cross_hits: int = 0
    cross_hit_specificity: float = 0.0
    cross_hit_precision: float = 0.0
    cross_hit_recall: float = 0.0
    cross_hit_f1: float = 0.0
    total_true_cross_hits: int = 0
    total_cross_hit_reads_mapped: int = 0
    cross_hit_counts_per_class: list[dict] | None = None
    cross_hit_reads_per_class: list[dict] | None = None
    spurious_hit_counts_per_class: list[dict] | None = None
    spurious_hit_reads_per_class: list[dict] | None = None
    reads_simulated_per_class: list[dict] | None = None

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "CrossHitMetrics":
        return cls(**data)


@dataclass
class DatasetResult:
    data_set: str
    input_df: pd.DataFrame
    sample: str = ""
    input_taxid_count: int = 0
    output_raw: int = 0
    output_taxid_count: int = 0
    output_cov_filtered: int = 0
    predicted_clades_pre: int = 0
    predicted_clades_post: int = 0
    predicted_clades_fixed: int = 0

    input_composition: dict | None = None
    input_read_counts: list[dict] | None = None
    reads_simulated_per_class: list[dict] | None = None
    precision: PrecisionMetrics = field(default_factory=PrecisionMetrics)
    recall: RecallMetrics = field(default_factory=RecallMetrics)
    cross_hit: CrossHitMetrics = field(default_factory=CrossHitMetrics)

    spurious_composition: dict | None = None
    cross_hit_composition: dict | None = None

    def to_dict(self) -> dict:
        result = {
            "data_set": self.data_set,
            "input_df": self.input_df.to_dict(orient="records"),
            "sample": self.sample,
            "input_taxid_count": self.input_taxid_count,
            "input_composition": self.input_composition,
            "input_read_counts": self.input_read_counts,
            "reads_simulated_per_class": self.reads_simulated_per_class,
            "output_raw": self.output_raw,
            "output_taxid_count": self.output_taxid_count,
            "output_cov_filtered": self.output_cov_filtered,
            "predicted_clades_pre": self.predicted_clades_pre,
            "predicted_clades_post": self.predicted_clades_post,
            "predicted_clades_fixed": self.predicted_clades_fixed,
            "precision": self.precision.to_dict(),
            "recall": self.recall.to_dict(),
            "cross_hit": self.cross_hit.to_dict(),
        }
        if self.spurious_composition is not None:
            result["spurious_composition"] = self.spurious_composition
        if self.cross_hit_composition is not None:
            result["cross_hit_composition"] = self.cross_hit_composition
        return result

    @classmethod
    def from_dict(cls, data: dict) -> "DatasetResult":
        return cls(
            data_set=data["data_set"],
            input_df=pd.DataFrame(data.get("input_df", [])),
            sample=data.get("sample", ""),
            input_taxid_count=data.get("input_taxid_count", 0),
            input_composition=data.get("input_composition"),
            input_read_counts=data.get("input_read_counts"),
            reads_simulated_per_class=data.get("reads_simulated_per_class"),
            output_raw=data.get("output_raw", 0),
            output_taxid_count=data.get("output_taxid_count", 0),
            output_cov_filtered=data.get("output_cov_filtered", 0),
            predicted_clades_fixed=data.get("predicted_clades_fixed", 0),
            predicted_clades_pre=data.get("predicted_clades_pre", 0),
            predicted_clades_post=data.get("predicted_clades_post", 0),
            precision=PrecisionMetrics.from_dict(data.get("precision", {})),
            recall=RecallMetrics.from_dict(data.get("recall", {})),
            cross_hit=CrossHitMetrics.from_dict(data.get("cross_hit", {})),
            spurious_composition=data.get("spurious_composition"),
            cross_hit_composition=data.get("cross_hit_composition"),
        )


@dataclass
class BatchEvaluationResult:
    input_df: pd.DataFrame = field(default_factory=pd.DataFrame)
    test_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    summary_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    spurious_composition: pd.DataFrame = field(default_factory=pd.DataFrame)
    cross_hit_composition: pd.DataFrame = field(default_factory=pd.DataFrame)

    metadata: dict = field(default_factory=dict)

    def to_json(self, filepath: str) -> None:
        """Save complete results to JSON file."""
        output = {
            "generated_at": datetime.now().isoformat(),
            "metadata": self.metadata,
            "test_results": self.test_results.to_dict(orient="records"),
            "summary_results": self.summary_results.to_dict(orient="records"),
            "spurious_composition": self.spurious_composition.to_dict(orient="records")
            if not self.spurious_composition.empty
            else [],
            "cross_hit_composition": self.cross_hit_composition.to_dict(orient="records")
            if not self.cross_hit_composition.empty
            else [],
        }
        with open(filepath, "w") as f:
            json.dump(output, f, indent=2, default=str)

    @classmethod
    def from_json(cls, filepath: str) -> "BatchEvaluationResult":
        """Load results from JSON file."""
        with open(filepath) as f:
            data = json.load(f)

        return cls(
            test_results=pd.DataFrame(data.get("test_results", [])),
            summary_results=pd.DataFrame(data.get("summary_results", [])),
            spurious_composition=pd.DataFrame(data.get("spurious_composition", [])),
            cross_hit_composition=pd.DataFrame(data.get("cross_hit_composition", [])),
            metadata=data.get("metadata", {}),
        )

    def save_tsv(self, output_dir: str) -> None:
        """Save DataFrames to TSV files (legacy format)."""
        self.input_df.to_csv(f"{output_dir}/test_datasets_input_df.tsv", sep="\t", index=False)
        self.test_results.to_csv(f"{output_dir}/test_datasets_overall_precision.tsv", sep="\t", index=False)
        self.summary_results.to_csv(f"{output_dir}/test_datasets_summary_results.tsv", sep="\t", index=False)
        self.spurious_composition.to_csv(f"{output_dir}/test_datasets_spurious_composition.tsv", sep="\t", index=False)
        self.cross_hit_composition.to_csv(
            f"{output_dir}/test_datasets_cross_hit_composition.tsv", sep="\t", index=False
        )

    def write_agent_output(self, filepath: str) -> None:
        """Save comprehensive evaluation results in agent-parseable JSON format.

        Includes aggregate summary statistics (mean, median, std, quartiles)
        alongside per-dataset records for easy downstream parsing.
        """
        summary = self.summary_results

        def _describe(s):
            if s.empty:
                return {}
            return {
                "mean": float(s.mean()),
                "median": float(s.median()),
                "std": float(s.std()),
                "q1": float(s.quantile(0.25)),
                "q3": float(s.quantile(0.75)),
                "min": float(s.min()),
                "max": float(s.max()),
            }

        precision_cols = [
            "precision_best_match",
            "purity",
            "purity_cov_filtered",
            "precision_clade_full",
            "precision_clade_post_cleanup",
            "precision_clade_fixed",
        ]
        recall_cols = [
            "recall_baseline",
            "recall_baseline_cov_filtered",
            "recall_clade_pre_cleanup",
            "recall_clade_post_cleanup",
            "recall_after_recall_filter",
            "recall_fixed_max_12",
        ]
        cross_hit_cols = [
            "cross_hit_precision",
            "cross_hit_recall",
            "cross_hit_f1",
            "cross_hit_specificity",
            "predicted_cross_hits",
            "total_true_cross_hits",
        ]

        aggregate = {}
        if not summary.empty:
            for col in precision_cols + recall_cols:
                if col in summary.columns:
                    aggregate[col] = _describe(summary[col])
            for col in cross_hit_cols:
                if col in summary.columns:
                    aggregate[col] = _describe(summary[col])

            if all(c in summary.columns for c in ["recall_baseline", "precision_best_match"]):
                prob = summary["recall_baseline"] * summary["precision_best_match"]
                aggregate["prob_find_true"] = _describe(prob)
            if all(c in summary.columns for c in ["recall_baseline", "purity"]):
                prob = summary["recall_baseline"] * summary["purity"]
                aggregate["prob_find_any"] = _describe(prob)
            if all(c in summary.columns for c in ["recall_clade_post_cleanup", "precision_clade_post_cleanup"]):
                prob = summary["recall_clade_post_cleanup"] * summary["precision_clade_post_cleanup"]
                aggregate["prob_find_true_clade_clean"] = _describe(prob)

        output = {
            "generated_at": datetime.now().isoformat(),
            "pipeline_version": self.metadata.get("pipeline_version", "unknown"),
            "metadata": self.metadata,
            "aggregate_summary": aggregate if aggregate else None,
            "test_results": self.test_results.to_dict(orient="records"),
            "summary_results": summary.to_dict(orient="records") if not summary.empty else [],
            "spurious_composition": self.spurious_composition.to_dict(orient="records")
            if not self.spurious_composition.empty
            else [],
            "cross_hit_composition": self.cross_hit_composition.to_dict(orient="records")
            if not self.cross_hit_composition.empty
            else [],
        }

        with open(filepath, "w") as f:
            json.dump(output, f, indent=2, default=str)

    def get_summary_stats(self) -> dict:
        """Compute summary statistics for precision metrics."""
        precision_cols = [
            "precision_best_match",
            "purity",
            "purity_cov_filtered",
            "precision_clade_full",
            "precision_clade_post_cleanup",
            "precision_clade_fixed",
        ]

        available_cols = [c for c in precision_cols if c in self.summary_results.columns]

        if not available_cols:
            return {}

        stats_df = self.summary_results[available_cols].describe()
        return stats_df.to_dict()

    def get_dataset_count(self) -> int:
        """Get number of datasets processed."""
        return self.summary_results["data_set"].nunique() if "data_set" in self.summary_results.columns else 0

    def save_metadata(self, output_dir: str) -> None:
        """Save pipeline metadata to TSV.

        Writes pipeline_metadata.tsv with columns: metric, value.
        Includes dataset counts and names of skipped/failed datasets.
        """
        rows = [
            ("total_datasets", self.metadata.get("total_datasets", 0)),
            ("successful", self.metadata.get("successful", 0)),
            ("failed", self.metadata.get("failed", 0)),
            (
                "skipped_count",
                self.metadata.get("skipped_count", len(self.metadata.get("skipped", []))),
            ),
        ]

        skipped = self.metadata.get("skipped", [])
        if skipped:
            rows.append(("skipped_datasets", ";".join(skipped)))

        errors = self.metadata.get("errors", [])
        if errors or self.metadata.get("failed_datasets"):
            failed_datasets = self.metadata.get("failed_datasets") or ";".join(errors)
            rows.append(("failed_datasets", failed_datasets))

        write_pipeline_metadata(rows, output_dir)


def create_empty_result() -> BatchEvaluationResult:
    """Create an empty BatchEvaluationResult."""
    return BatchEvaluationResult(
        test_results=pd.DataFrame(columns=["recall_baseline", "precision_fixed", "precision_clade_post", "data_set"]),
        summary_results=pd.DataFrame(),
        spurious_composition=pd.DataFrame(),
        cross_hit_composition=pd.DataFrame(),
    )
