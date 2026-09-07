"""
Tests for the recall-loss fixes:

* F1 classification-aware ``compute_recall``
* F3 exact-taxid shortcut in ``find_best_match``
* F4 ``recall_raw`` no longer implicitly coverage-filtered
* F5 symmetric ``compare_lineages`` denominator
* Phase 1 ``diagnose_recall_gap`` bucketing
"""

from __future__ import annotations

import os

import pandas as pd
import pytest

os.environ.setdefault("NCBI_EMAIL", "test@example.com")


# ---------------------------------------------------------------------------
# F5: compare_lineages
# ---------------------------------------------------------------------------


class TestSymmetricCompareLineages:
    def test_deeper_leaf_lineage_scores_full_match(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        leaf = "Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes; Nidovirales; Coronaviridae; Betacoronavirus; Sarbecovirus"
        inp = "Viruses; Riboviria; Orthornavirae; Pisuviricota; Pisoniviricetes; Nidovirales; Coronaviridae; Betacoronavirus"
        score, _ = compare_lineages(leaf, inp)
        assert score == 1.0
        # symmetric
        assert compare_lineages(inp, leaf)[0] == 1.0

    def test_previous_asymmetry_would_have_dropped_match(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        # 7 shared ranks, leaf has 3 extra: old score = 7/10 = 0.7 (borderline); new = 1.0
        shared = "; ".join(f"R{i}" for i in range(7))
        leaf = shared + "; X; Y; Z"
        score, _ = compare_lineages(leaf, shared)
        assert score == 1.0

    def test_partial_match_unchanged(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        assert compare_lineages("A; B; C; D", "A; B; X; Y")[0] == 0.5

    def test_one_empty_lineage_does_not_divide_by_zero(self):
        from metagenomics_utils.ncbi_tools import compare_lineages

        score, level = compare_lineages("A; B", "")
        assert score == 0.0
        assert level is None


# ---------------------------------------------------------------------------
# F3: exact taxid shortcut
# ---------------------------------------------------------------------------


class _FakeWrapper:
    TAX_SPECIES = "species"

    def __init__(self, lineages=None):
        self.lineages = lineages or {}
        self.calls = 0

    def compare_lineages_relative(self, taxid1, taxid2):
        self.calls += 1
        return 0.0, None

    def get_name(self, taxid):
        return f"name_{taxid}"

    @staticmethod
    def level_is_below(level, threshold_level="species"):
        from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper

        return NCBITaxonomistWrapper.level_is_below(level, threshold_level)


class TestFindBestMatchExactTaxid:
    def test_identical_taxid_matches_without_lineage(self):
        from metagenomics_utils.overlap_manager.node_stats import find_best_match

        wrapper = _FakeWrapper()
        taxid, level, score = find_best_match(11676, [10359, 11676, 10245], wrapper)
        assert taxid == 11676
        assert level == "species"
        assert score == 1.0

    def test_string_vs_int_taxids_are_equal(self):
        from metagenomics_utils.overlap_manager.node_stats import find_best_match

        wrapper = _FakeWrapper()
        taxid, _, score = find_best_match("11676", [11676], wrapper)
        assert taxid == 11676
        assert score == 1.0

    def test_nan_taxid_does_not_match(self):
        from metagenomics_utils.overlap_manager.node_stats import find_best_match

        wrapper = _FakeWrapper()
        taxid, level, score = find_best_match(float("nan"), [11676], wrapper)
        assert taxid is None
        assert score == 0.0

    def test_update_df_best_match_uses_shortcut(self):
        from metagenomics_utils.overlap_manager.node_stats import update_df_best_match

        wrapper = _FakeWrapper()
        row = pd.Series({"taxid": 11676})
        out = update_df_best_match(row, [11676], wrapper, min_score=0.7)
        assert out["best_match_taxid"] == 11676
        assert out["best_match_score"] == 1.0
        assert out["name"] == "name_11676"


# ---------------------------------------------------------------------------
# F4 + F1: compute_recall
# ---------------------------------------------------------------------------


def _m_stats(rows):
    return pd.DataFrame(
        rows,
        columns=["best_match_taxid", "best_match_is_best", "is_trash", "coverage"],
    )


class TestComputeRecall:
    def test_zero_coverage_best_counts_in_raw_not_cov(self):
        from deployment.model_evaluation.metrics import compute_recall

        m = _m_stats([[1, True, False, 0.0], [2, True, False, 5.0]])
        inputs = pd.DataFrame({"taxid": [1, 2, 3]})
        raw, cov, kept, kept_cov = compute_recall(m, inputs)
        assert raw == pytest.approx(2 / 3)
        assert cov == pytest.approx(1 / 3)
        assert kept == [1, 2]
        assert kept_cov == [2]

    def test_classification_evidence_extends_raw_only(self):
        from deployment.model_evaluation.metrics import compute_recall

        m = _m_stats([[2, True, False, 5.0]])
        inputs = pd.DataFrame({"taxid": [1, 2, 3]})
        raw, cov, kept, kept_cov = compute_recall(m, inputs, extra_detected_taxids={3, 99})
        assert raw == pytest.approx(2 / 3)
        assert cov == pytest.approx(1 / 3)
        assert kept == [2, 3]
        assert kept_cov == [2]

    def test_extra_taxids_outside_input_are_ignored(self):
        from deployment.model_evaluation.metrics import compute_recall

        m = _m_stats([])
        inputs = pd.DataFrame({"taxid": [1, 2]})
        raw, cov, kept, _ = compute_recall(m, inputs, extra_detected_taxids={5, 6})
        assert raw == 0.0
        assert cov == 0.0
        assert kept == []

    def test_default_behaviour_unchanged(self):
        from deployment.model_evaluation.metrics import compute_recall

        m = _m_stats([[1, True, False, 1.0], [2, False, True, 1.0]])
        inputs = pd.DataFrame({"taxid": [1, 2]})
        raw, cov, kept, _ = compute_recall(m, inputs)
        assert raw == 0.5
        assert cov == 0.5
        assert kept == [1]


class TestBestMatchMarking:
    def test_zero_coverage_row_is_marked_best_when_only_candidate(self):
        """F4: `_compute_best_matches` should no longer skip coverage==0 rows."""
        import numpy as np

        from metagenomics_utils.overlap_manager import node_stats

        m_stats = pd.DataFrame(
            {
                "taxid": [11676, 10359],
                "assid": ["NC_001802.1", "NC_006273.2"],
                "coverage": [0.0, 12.5],
                "error_rate": [0.0, 0.01],
                "numreads": [0, 100],
                "total_uniq_reads": [10, 100],
            }
        )

        class _OM:
            leaves = ["NC_001802.1", "NC_006273.2"]

        monkey = pytest.MonkeyPatch()
        try:
            monkey.setattr(node_stats, "dataframe_update_with_lineage", lambda df, w: df)
            monkey.setattr(
                "metagenomics_utils.overlap_manager.manager.merge_by_assembly_ID", lambda df: df, raising=True
            )
            wrapper = _FakeWrapper(lineages={11676: {}, 10359: {}})
            out = node_stats._compute_best_matches(m_stats, [11676, 10359], wrapper, _OM(), 0.7)
        finally:
            monkey.undo()

        best = out[out["best_match_is_best"] == True]
        assert set(int(t) for t in best["best_match_taxid"]) == {11676, 10359}
        zero_cov_best = out[(out["coverage"] == 0.0) & (out["best_match_is_best"] == True)]
        assert len(zero_cov_best) == 1
        assert np.isclose(zero_cov_best["best_match_score"].iloc[0], 1.0)


# ---------------------------------------------------------------------------
# Result model round-trip
# ---------------------------------------------------------------------------


def test_recall_metrics_from_dict_tolerates_missing_and_extra_keys():
    from deployment.model_evaluation.result_models import RecallMetrics

    rm = RecallMetrics.from_dict({"recall_raw": 0.4, "unknown_key": 1})
    assert rm.recall_raw == 0.4
    assert rm.recall_assembly_raw == 0.0
    assert rm.recall_classification_raw == 0.0
    assert RecallMetrics.from_dict(rm.to_dict()) == rm


# ---------------------------------------------------------------------------
# Classification-aware helpers (dataset_processor / analysis_data_extractor)
# ---------------------------------------------------------------------------


def _write_classification(tmp_path, dataset, rows):
    d = tmp_path / dataset / "classification"
    d.mkdir(parents=True)
    pd.DataFrame(rows, columns=["taxid", "description", "uniq_reads", "classifiers"]).to_csv(
        d / f"{dataset}_merged_classification.tsv", sep="\t", index=False
    )


def test_load_classifier_detected_taxids_filters_zero_reads(tmp_path):
    from deployment.model_evaluation.dataset_processor import load_classifier_detected_taxids

    _write_classification(tmp_path, "ds1", [[1, "a", 10, "k2"], [2, "b", 0, "k2"], [3, "c", 5, "cf"]])
    assert load_classifier_detected_taxids(str(tmp_path), "ds1") == {1, 3}
    assert load_classifier_detected_taxids(str(tmp_path), "missing") == set()


def test_collect_all_matched_taxids_reads_leaf_taxids(tmp_path):
    from deployment.model_evaluation.analysis_data_extractor import collect_all_matched_taxids

    for ds, taxids in {"ds1": [11676, 10359], "ds2": [10359, 2697049.0]}.items():
        d = tmp_path / ds / "output"
        d.mkdir(parents=True)
        pd.DataFrame({"taxid": taxids, "assembly_accession": ["X"] * len(taxids)}).to_csv(
            d / "matched_assemblies.tsv", sep="\t", index=False
        )
    (tmp_path / "ds3").mkdir()  # no output dir -> skipped silently
    assert collect_all_matched_taxids(str(tmp_path), ["ds1", "ds2", "ds3"]) == [10359, 11676, 2697049]


def test_collect_recall_data_reports_assembly_and_classification(tmp_path):
    from deployment.model_evaluation.analysis_data_extractor import collect_recall_data

    m_stats = pd.DataFrame({"best_match_is_best": [True], "best_match_taxid": [1.0]})
    inputs = pd.DataFrame({"taxid": [1, 2, 3], "reads": [10, 20, 30], "mutation_rate": [0.0, 0.1, 0.2]})

    class _W:
        def get_level(self, tid, level):
            return "Ord"

    recs = collect_recall_data(m_stats, inputs, "ds", _W(), classifier_detected_taxids={2})
    by_tid = {r["taxid"]: r for r in recs}
    assert (by_tid[1]["recalled"], by_tid[1]["recalled_assembly"], by_tid[1]["recalled_classification"]) == (1, 1, 0)
    assert (by_tid[2]["recalled"], by_tid[2]["recalled_assembly"], by_tid[2]["recalled_classification"]) == (1, 0, 1)
    assert (by_tid[3]["recalled"], by_tid[3]["recalled_assembly"], by_tid[3]["recalled_classification"]) == (0, 0, 0)


# ---------------------------------------------------------------------------
# Phase 1: diagnose_recall_gap
# ---------------------------------------------------------------------------


class TestDiagnoseRecallGap:
    def _frames(self):
        inputs = pd.DataFrame(
            {
                "taxid": [1, 2, 3, 4, 5],
                "reads": [10, 100, 1000, 50, 500],
                "mutation_rate": [0.0, 0.1, 0.2, 0.05, 0.3],
            }
        )
        classification = pd.DataFrame({"taxid": [2, 3, 4, 5], "uniq_reads": [5, 50, 3, 20]})
        matched = pd.DataFrame(
            {
                "taxid": [3, 4, 5],
                "assembly_accession": ["NC_000003.1", "NC_000004.1", "NC_000005.1"],
                "assembly_file": ["/s/3/3_NC_000003.1_sequence.fasta.gz", "/s/4/x.fa.gz", "/s/5/y.fa.gz"],
            }
        )
        coverage = pd.DataFrame(
            {
                "#rname": ["NC_000003.1", "NC_000004.1"],
                "numreads": [0, 42],
                "coverage": [0.0, 30.0],
                "file": ["ds_plan_3_NC_000003.1.fna.sorted.coverage.txt", "ds_plan_4_NC_000004.1.fna.sorted.coverage.txt"],
            }
        )
        return inputs, classification, matched, coverage

    def test_buckets(self):
        from deployment.model_evaluation.analysis_scripts.diagnose_recall_gap import classify_recall_gap

        out = classify_recall_gap(*self._frames()).set_index("taxid")
        assert out.loc[1, "bucket"] == "absent_classification"
        assert out.loc[2, "bucket"] == "classified_no_assembly"
        assert out.loc[3, "bucket"] == "assembly_zero_reads"
        assert out.loc[4, "bucket"] == "assembly_with_reads"
        assert out.loc[5, "bucket"] == "assembly_no_coverage"
        assert out.loc[4, "mapped_reads"] == 42
        assert out.loc[2, "classifier_uniq_reads"] == 5

    def test_missing_files_fall_back_gracefully(self):
        from deployment.model_evaluation.analysis_scripts.diagnose_recall_gap import classify_recall_gap

        inputs, *_ = self._frames()
        out = classify_recall_gap(inputs, None, None, None)
        assert (out["bucket"] == "absent_classification").all()

    def test_summary_shares_sum_to_one(self):
        from deployment.model_evaluation.analysis_scripts.diagnose_recall_gap import (
            classify_recall_gap,
            summarise_buckets,
        )

        out = classify_recall_gap(*self._frames())
        summary = summarise_buckets(out)
        assert summary["share"].sum() == pytest.approx(1.0)
        assert summary["count"].sum() == 5

    def test_driver_end_to_end(self, tmp_path):
        from deployment.model_evaluation.analysis_scripts.diagnose_recall_gap import main

        inputs, classification, matched, coverage = self._frames()
        study = tmp_path / "study"
        ds = study / "ds1"
        (ds / "input").mkdir(parents=True)
        (ds / "classification").mkdir()
        (ds / "output").mkdir()
        inputs.to_csv(ds / "input" / "ds1.tsv", sep="\t", index=False)
        classification.to_csv(ds / "classification" / "ds1_merged_classification.tsv", sep="\t", index=False)
        matched.to_csv(ds / "output" / "matched_assemblies.tsv", sep="\t", index=False)
        coverage.to_csv(ds / "output" / "merged_coverage_statistics.tsv", sep="\t", index=False)

        out_dir = tmp_path / "diag"
        rc = main(["--study-output", str(study), "--output-dir", str(out_dir), "--no-plot"])
        assert rc == 0
        per_taxid = pd.read_csv(out_dir / "recall_gap_per_taxid.tsv", sep="\t")
        assert len(per_taxid) == 5
        assert set(per_taxid["data_set"]) == {"ds1"}
        assert (out_dir / "recall_gap_summary.tsv").exists()
        assert (out_dir / "recall_gap_per_dataset.tsv").exists()


# ---------------------------------------------------------------------------
# F2: AssemblyStore local lookup
# ---------------------------------------------------------------------------


class TestAssemblyStoreLocalLookup:
    @pytest.fixture
    def store(self, tmp_path, monkeypatch):
        from metagenomics_utils import reference_utils

        monkeypatch.setattr(reference_utils, "NCBITools", lambda: object())
        return reference_utils.AssemblyStore(str(tmp_path / "store"))

    def test_exact_filename(self, store):
        from metagenomics_utils.ncbi_tools import Passport

        d = os.path.join(store.store_path, "562")
        os.makedirs(d)
        open(os.path.join(d, "562_NC_000913.3_sequence.fasta.gz"), "w").close()
        res = store.retrieve_local_assembly(Passport(taxid="562", accession="NC_000913.3"))
        assert res is not None
        assert res.file_path.endswith("562_NC_000913.3_sequence.fasta.gz")
        assert res.accession == "NC_000913.3"

    def test_alternate_filename_in_taxid_dir(self, store):
        from metagenomics_utils.ncbi_tools import Passport

        d = os.path.join(store.store_path, "562")
        os.makedirs(d)
        open(os.path.join(d, "NC_000913.3.fna.gz"), "w").close()
        res = store.retrieve_local_assembly(Passport(taxid="562", accession="NC_000913.3"))
        assert res is not None
        assert res.file_path.endswith("NC_000913.3.fna.gz")

    def test_no_accession_uses_any_sequence_file(self, store):
        from metagenomics_utils.ncbi_tools import Passport

        d = os.path.join(store.store_path, "562")
        os.makedirs(d)
        open(os.path.join(d, "562_GCF_000005845.2_sequence.fasta.gz"), "w").close()
        res = store.retrieve_local_assembly(Passport(taxid="562"))
        assert res is not None
        assert res.accession == "GCF_000005845.2"

    def test_accession_found_under_other_taxid_dir(self, store):
        from metagenomics_utils.ncbi_tools import Passport

        d = os.path.join(store.store_path, "83333")
        os.makedirs(d)
        open(os.path.join(d, "83333_NC_000913.3_sequence.fasta.gz"), "w").close()
        res = store.retrieve_local_assembly(Passport(taxid="562", accession="NC_000913.3"))
        assert res is not None
        assert "83333" in res.file_path

    def test_missing_returns_none(self, store):
        from metagenomics_utils.ncbi_tools import Passport

        assert store.retrieve_local_assembly(Passport(taxid="1", accession="NC_1.1")) is None
