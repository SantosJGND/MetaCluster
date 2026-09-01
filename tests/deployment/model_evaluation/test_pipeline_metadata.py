"""Tests for the pipeline/cohort metadata contract (changes_requests.md)."""

import argparse
import json

import pandas as pd
import pytest

from deployment.model_evaluation.batch_evaluator import BatchEvaluator
from deployment.model_evaluation.exceptions import EvaluatorError
from deployment.model_evaluation.result_models import (
    BatchEvaluationResult,
    CrossHitMetrics,
    DatasetResult,
    PrecisionMetrics,
    RecallMetrics,
    write_pipeline_metadata,
)


def _meta_map(df: pd.DataFrame) -> dict:
    """Map metric column to value, coercing numeric-looking values to int."""
    out = {}
    for metric, value in zip(df["metric"], df["value"]):
        try:
            out[metric] = int(value)
        except (TypeError, ValueError):
            out[metric] = str(value)
    return out


def test_write_pipeline_metadata_atomic(tmp_path):
    rows = [("total_attempted", 10), ("extracted", 8), ("dropped", 2)]
    path = write_pipeline_metadata(rows, str(tmp_path))
    df = pd.read_csv(path, sep="\t")
    assert list(df.columns) == ["metric", "value"]
    m = _meta_map(df)
    assert m == {"total_attempted": 10, "extracted": 8, "dropped": 2}
    leftovers = [f for f in tmp_path.iterdir() if f.name.endswith(".tmp.")]
    assert leftovers == []


def test_save_metadata_rows(tmp_path):
    meta = {
        "total_datasets": 5,
        "successful": 2,
        "failed": 1,
        "skipped": ["s1"],
        "skipped_count": 1,
        "errors": ["boom"],
        "failed_datasets": "boom",
    }
    res = BatchEvaluationResult()
    res.metadata = meta
    res.save_metadata(str(tmp_path))

    m = _meta_map(pd.read_csv(tmp_path / "pipeline_metadata.tsv", sep="\t"))
    assert m["total_datasets"] == 5
    assert m["successful"] == 2
    assert m["failed"] == 1
    assert m["skipped_count"] == 1
    assert m["skipped_datasets"] == "s1"
    assert m["failed_datasets"] == "boom"


def test_to_json_tsv_parity(tmp_path):
    meta = {
        "total_datasets": 5,
        "successful": 2,
        "failed": 1,
        "skipped": ["s1"],
        "skipped_count": 1,
        "errors": ["boom"],
        "failed_datasets": "boom",
    }
    res = BatchEvaluationResult(metadata=meta)
    res.to_json(str(tmp_path / "evaluation_results.json"))
    res.save_metadata(str(tmp_path))

    data = json.loads((tmp_path / "evaluation_results.json").read_text())
    jmeta = data["metadata"]
    tmeta = _meta_map(pd.read_csv(tmp_path / "pipeline_metadata.tsv", sep="\t"))

    assert jmeta["skipped_count"] == 1
    assert jmeta["failed_datasets"] == "boom"
    assert jmeta["skipped"] == ["s1"]
    assert jmeta["successful"] == tmeta["successful"]
    assert jmeta["failed"] == tmeta["failed"]
    assert jmeta["skipped_count"] == tmeta["skipped_count"]
    assert jmeta["failed_datasets"] == tmeta["failed_datasets"]


def test_aggregate_results_metadata_no_success():
    evaluator = object.__new__(BatchEvaluator)
    res = BatchEvaluator._aggregate_results(
        evaluator,
        results=[],
        errors=[EvaluatorError("boom")],
        skipped=["skipped_ds"],
    )
    meta = res.metadata
    assert meta["successful"] == 0
    assert meta["failed"] == 1
    assert meta["skipped_count"] == 1
    assert meta["skipped"] == ["skipped_ds"]
    assert meta["failed_datasets"] == "boom"


def test_aggregate_results_metadata_success():
    result = DatasetResult(data_set="d1", input_df=pd.DataFrame())
    result.precision = PrecisionMetrics()
    result.recall = RecallMetrics()
    result.cross_hit = CrossHitMetrics()

    evaluator = object.__new__(BatchEvaluator)
    res = BatchEvaluator._aggregate_results(evaluator, results=[result], errors=[], skipped=[])
    meta = res.metadata
    assert meta["total_datasets"] == 1
    assert meta["successful"] == 1
    assert meta["failed"] == 0
    assert meta["skipped_count"] == 0
    assert meta["skipped"] == []
    assert meta["failed_datasets"] == ""


def test_extractor_emits_cohort_metadata(tmp_path, monkeypatch):
    import deployment.model_evaluation.analysis_data_extractor as adx

    study = tmp_path / "study"
    for name in ["ds_ok", "ds_skip", "ds_fail", "gap_folder"]:
        (study / name).mkdir(parents=True)
    out = tmp_path / "out"
    out.mkdir()

    monkeypatch.setattr(adx, "discover_datasets", lambda p: ["ds_ok", "ds_skip", "ds_fail"])

    def fake_process(ds, *args, **kwargs):
        if ds == "ds_ok":
            record = {
                "data_set": ds,
                "input_taxid_count": 2,
                "output_raw": 5,
                "output_taxid_count": 4,
                "tp_count": 3,
                "cross_hit_count": 1,
                "spurious_count": 0,
                "overall_recall": 0.75,
                "recall_cov_filtered": 0.5,
                "last_best_match_relindex": 0.4,
            }
            return record, [], [], [], []
        if ds == "ds_skip":
            return None, None, None, None, None
        raise ValueError("boom")

    monkeypatch.setattr(adx, "process_dataset", fake_process)

    class FakeNCBI:
        def __init__(self, *args, **kwargs):
            self.lineages = {}

        def resolve_lineages(self, taxids):
            return None

    monkeypatch.setattr(adx, "NCBITaxonomistWrapper", FakeNCBI)

    namespace = argparse.Namespace(
        study_output=str(study),
        ncbi_db="taxa.db",
        output_dir=str(out),
        cross_hit_threshold=0.3,
        min_taxonomic_score=0.7,
        explanatory=False,
    )
    monkeypatch.setattr(adx, "parse_args", lambda: namespace)

    adx.main()

    m = _meta_map(pd.read_csv(out / "pipeline_metadata.tsv", sep="\t"))
    assert m["total_attempted"] == 4
    assert m["extracted"] == 1
    assert m["dropped"] == 3
    assert m["failed"] == 1
    assert m["skipped"] == 1
    assert m["study_gaps"] == 1
    assert m["skipped_datasets"] == "ds_skip"
    assert m["failed_datasets"] == "boom"
    assert m["study_gap_datasets"] == "gap_folder"

    per_dataset = pd.read_csv(out / "per_dataset_metrics.tsv", sep="\t")
    assert len(per_dataset) == m["extracted"]

    assert m["dropped"] == m["skipped"] + m["failed"] + m["study_gaps"]