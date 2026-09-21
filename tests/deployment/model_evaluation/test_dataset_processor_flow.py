"""
Tests for ``DatasetProcessor.process()`` pre/post-cleanup OM chaining.

Verifies the coherence fix (see README § "Pre- vs Post-Cleanup Semantics"):

* ``_apply_crosshit_cleanup`` mutates the **same** OverlapManager that
  ``_predict_clades_postcleanup`` consumes — so cross-hit cleanup feeds the
  reported post-cleanup clade columns.
* ``apply_recall_filter=False`` (``--no-apply-recall-filter``) skips the recall
  filter entirely and passes the original full OverlapManager to post cleanup.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from deployment.model_evaluation.config import EvaluatorConfig
from deployment.model_evaluation.dataset_processor import DatasetProcessor

INPUT_SUMMARY = pd.DataFrame(
    {
        "sample": ["s1"],
        "taxid": [9606],
        "reads": [100],
        "mutation_rate": [0.01],
    }
)


@pytest.fixture(autouse=True)
def _ncbi_email():
    import os

    os.environ.setdefault("NCBI_EMAIL", "test@example.com")


def _make_processor(apply_recall_filter: bool, enable_cross_hit: bool) -> DatasetProcessor:
    base = Path(__file__).parent / "fixtures"
    config = EvaluatorConfig(
        study_output_filepath=base / "minimal_study",
        taxid_plan_filepath=base / "minimal_taxid_plan.tsv",
        analysis_output_filepath=base / "tmp_analysis",
        apply_recall_filter=apply_recall_filter,
        enable_cross_hit=enable_cross_hit,
        use_cache=False,
    )
    processor = DatasetProcessor.__new__(DatasetProcessor)
    processor.config = config
    processor.ncbi = object()
    processor.min_uniq_reads = 1
    processor.max_missing_pct = 5.0
    return processor


def _wire_process(processor: DatasetProcessor, om, filtered_om, result) -> None:
    """Replace every step of ``process()`` with mocks returning sentinels."""
    processor._load_data = MagicMock(return_value=(om, INPUT_SUMMARY))
    processor._compute_baseline_metrics = MagicMock(return_value=result)
    processor._predict_clades_precleanup = MagicMock(return_value=result)
    processor._apply_recall_filter = MagicMock(return_value=(result, filtered_om))
    processor._apply_fixed_filter = MagicMock(return_value=(result, MagicMock()))
    processor._apply_crosshit_cleanup = MagicMock(return_value=result)
    processor._predict_clades_postcleanup = MagicMock(return_value=result)


class TestPostCleanupOMChaining:
    def test_crosshit_cleanup_and_postcleanup_share_the_recall_filtered_om(self):
        om, filtered_om, result = object(), object(), object()
        processor = _make_processor(apply_recall_filter=True, enable_cross_hit=True)
        _wire_process(processor, om, filtered_om, result)

        with patch("deployment.model_evaluation.dataset_processor.get_m_stats_matrix") as g, patch(
            "deployment.model_evaluation.dataset_processor.passes_assembly_completeness"
        ) as p, patch("deployment.model_evaluation.dataset_processor.DatasetResult") as dr:
            g.return_value = MagicMock()
            p.return_value = True
            dr.return_value = result

            assert processor.process("ds1") is result

        processor._apply_crosshit_cleanup.assert_called_once_with("ds1", filtered_om, result)
        processor._predict_clades_postcleanup.assert_called_once_with("ds1", filtered_om, result)

    def test_crosshit_disabled_passes_recall_filtered_om_to_postcleanup_only(self):
        om, filtered_om, result = object(), object(), object()
        processor = _make_processor(apply_recall_filter=True, enable_cross_hit=False)
        _wire_process(processor, om, filtered_om, result)

        with patch("deployment.model_evaluation.dataset_processor.get_m_stats_matrix"), patch(
            "deployment.model_evaluation.dataset_processor.passes_assembly_completeness", return_value=True
        ), patch("deployment.model_evaluation.dataset_processor.DatasetResult", return_value=result):
            processor.process("ds1")

        processor._apply_crosshit_cleanup.assert_not_called()
        processor._predict_clades_postcleanup.assert_called_once_with("ds1", filtered_om, result)

    def test_recall_filter_disabled_uses_the_original_om(self):
        om, filtered_om, result = object(), object(), object()
        processor = _make_processor(apply_recall_filter=False, enable_cross_hit=True)
        _wire_process(processor, om, filtered_om, result)

        with patch("deployment.model_evaluation.dataset_processor.get_m_stats_matrix"), patch(
            "deployment.model_evaluation.dataset_processor.passes_assembly_completeness", return_value=True
        ), patch("deployment.model_evaluation.dataset_processor.DatasetResult", return_value=result):
            processor.process("ds1")

        processor._apply_recall_filter.assert_not_called()
        processor._apply_crosshit_cleanup.assert_called_once_with("ds1", om, result)
        processor._predict_clades_postcleanup.assert_called_once_with("ds1", om, result)


class TestCLIWiring:
    def test_config_defaults(self):
        base = Path(__file__).parent / "fixtures"
        config = EvaluatorConfig(
            study_output_filepath=base / "minimal_study",
            taxid_plan_filepath=base / "minimal_taxid_plan.tsv",
            analysis_output_filepath=base / "tmp_analysis",
        )
        assert config.apply_recall_filter is True

    def test_no_apply_recall_filter_flag(self, monkeypatch):
        from deployment.model_evaluation import evaluate

        monkeypatch.setattr(
            "sys.argv",
            [
                "evaluate.py",
                "--study_output_filepath",
                "tests/deployment/model_evaluation/fixtures/minimal_study",
                "--taxid_plan_filepath",
                "tests/deployment/model_evaluation/fixtures/minimal_taxid_plan.tsv",
                "--analysis_output_filepath",
                "/tmp/opencode/model_eval_test_out",
                "--no-apply-recall-filter",
            ],
        )
        args = evaluate.get_args()
        assert args.apply_recall_filter is False
