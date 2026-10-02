"""
Tests for ``reference_details`` — the per-reference and per-leaf tracking tables.

Covers the two things downstream analysts rely on:

* ``loss_stage`` is a strict first-match cascade in pipeline order, so the
  stage counts form a proper funnel (the first four stage names must stay
  compatible with ``analysis_scripts/diagnose_recall_gap.py``).
* The ``loss_*`` flags are *independent*, so a reference that fails several
  conditions at once still reports each of them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from deployment.model_evaluation.analysis_data_extractor import load_classifier_detected_taxids
from deployment.model_evaluation.reference_details import (
    LEGACY_INPUT_DF_COLUMNS,
    LOSS_FLAGS,
    LOSS_LOSSES,
    LOSS_STAGES,
    REFERENCE_DETAIL_COLUMNS,
    apply_clade_membership,
    assign_loss_stages,
    build_reference_leaf_table,
    load_reference_evidence,
    mark_clade_leaves,
    mark_filter_survival,
    order_reference_detail,
)

# ---------------------------------------------------------------------------
# fixtures / builders
# ---------------------------------------------------------------------------

# Four leaves exercising: best match with reads, a zero-read assembly, a
# trash-flagged leaf, and a non-best leaf.
M_STATS = pd.DataFrame(
    {
        "taxid": [10, 10, 20, 30],
        "leaf": ["L0", "L1", "L2", "L3"],
        "coverage": [11.0, 0.0, 5.0, 2.0],
        "covbases": [100, 0, 50, 20],
        "meanmapq": [60.0, 0.0, 50.0, 40.0],
        "numreads": [500, 0, 250, 100],
        "error_rate": [0.01, 0.0, 0.02, 0.03],
        "total_uniq_reads": [1000, 10, 900, 800],
        "best_match_taxid": [100, 100, 200, 300],
        "best_match_is_best": [True, False, False, True],
        "best_match_score": [1.0, 0.0, 0.0, 1.0],
        "best_match_level": ["species", "species", "species", "species"],
        "is_trash": [False, True, True, False],
        "is_crosshit": [False, False, False, False],
        "max_shared": [0.0, 0.9, 0.5, 0.0],
        "cross_hit_match": [None, 100, None, None],
        "crosshit_match_score": [0.0, 0.8, 0.0, 0.0],
        "crosshit_match_level": [None, "species", None, None],
        "crossh_hit_dist": [1.0, 0.1, 0.4, 1.0],
    },
    index=["A0", "A1", "A2", "A3"],
)

CLADES = pd.DataFrame(
    {
        "node": ["N0"],
        "n_leaves": [2],
        "leaves": [["L0", "L3"]],
        "node_precision": [0.5],
    }
)


def _write_study(root, data_set="ds1", *, classification: str, coverage: str):
    """Materialise the source tables ``load_reference_evidence`` reads.

    ``matched_assemblies.tsv`` is deliberately *not* written: its ``taxid`` column
    is the matched assembly's taxid, so the reference evidence must come from
    ``m_stats["best_match_taxid"]`` instead. See
    ``test_n_clean_matches_follows_the_recall_criterion``.
    """
    out = root / data_set / "output"
    cls = root / data_set / "classification"
    out.mkdir(parents=True)
    cls.mkdir(parents=True)
    (cls / f"{data_set}_merged_classification.tsv").write_text(classification)
    (out / "merged_coverage_statistics.tsv").write_text(coverage)
    return str(root)


CURRENT_CLASSIFICATION = (
    "taxid\tdescription\tuniq_reads\tclassifiers\n100\tspecies a\t17\tKraken2\n200\tspecies b\t5\tCentrifuge;Kraken2\n"
)
LEGACY_CLASSIFICATION = "description\ttaxid\tclassifier_centrifuge\tclassifier_kraken2\tclassification\nspecies a\t100\t-\tKraken2\t1\nspecies b\t200\tCentrifuge\t-\t1\n"
CURRENT_COVERAGE = "#rname\tcoverage\nA0\t11.0\nA2\t5.0\nA3\t2.0\n"
LEGACY_COVERAGE = "index\tcoverage\nA0\t11.0\nA2\t5.0\nA3\t2.0\n"


# ---------------------------------------------------------------------------
# schema tolerance
# ---------------------------------------------------------------------------


class TestMatchSemantics:
    """The evidence columns must agree with ``metrics.compute_recall``.

    ``compute_recall`` counts a reference as recalled when it has at least one
    leaf that is neither trash nor non-best
    (``is_trash == False and best_match_is_best == True``). The evidence table
    has to use exactly the same criterion, otherwise the loss funnel would
    contradict the reported recall.
    """

    @staticmethod
    def _recall_matched_taxids(m_stats: pd.DataFrame) -> set:
        clean = m_stats.dropna(subset=["best_match_taxid"])
        clean = clean[(clean["is_trash"] == False) & (clean["best_match_is_best"] == True)]  # noqa: E712
        return set(clean["best_match_taxid"].astype(int))

    def test_n_clean_matches_reproduces_the_recall_criterion(self, tmp_path):
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        taxids = [100, 200, 300]
        detail = load_reference_evidence(root, "ds1", M_STATS, taxids)

        assert set(detail.loc[detail["n_clean_matches"] > 0, "taxid"].astype(int)) == self._recall_matched_taxids(
            M_STATS
        )

    def test_only_trash_leaves_count_as_no_clean_match(self, tmp_path):
        """Taxid 100's L1 is trash, so it must not count as a clean match."""
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100])

        assert detail.loc[0, "n_leaves"] == 2
        assert detail.loc[0, "n_clean_matches"] == 1

    def test_n_leaves_counts_every_leaf_regardless_of_trash(self, tmp_path):
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100])

        assert detail.loc[0, "n_leaves"] == 2

    def test_coverage_detection_uses_the_coverage_table_accessions(self, tmp_path):
        root = _write_study(
            tmp_path,
            classification=CURRENT_CLASSIFICATION,
            coverage="#rname\tcoverage\nA0\t11.0\n",
        )
        detail = load_reference_evidence(root, "ds1", M_STATS, [100, 200, 300])

        # Only A0 has a coverage record; A1/A2/A3 do not.
        assert detail["n_covered_leaves"].tolist() == [1, 0, 0]
        assert detail["has_coverage_row"].tolist() == [True, False, False]

    def test_matched_assemblies_file_is_not_the_join_key(self, tmp_path):
        """A present-but-misleading matched_assemblies.tsv must be ignored.

        That file's ``taxid`` column is the matched assembly's taxid, so joining
        input taxids on it reported that no reference was ever matched. Taxids 10
        and 20 appear only in that file and must not acquire a match.
        """
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        (tmp_path / "ds1" / "output" / "matched_assemblies.tsv").write_text(
            "taxid\tdescription\tuniq_reads\tclassifiers\tassembly_accession\n"
            "10\tasm a\t100\tKraken2\tA0\n"
            "20\tasm b\t90\tKraken2\tA2\n"
        )
        detail = load_reference_evidence(root, "ds1", M_STATS, [100, 200, 300, 10, 20])

        # 100 -> one clean leaf; 200 -> its only leaf is trash; 300 -> one clean leaf.
        assert detail["n_clean_matches"].tolist() == [1, 0, 1, 0, 0]

    def test_absent_m_stats_reports_no_matches(self, tmp_path):
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", pd.DataFrame(), [100])

        assert detail.loc[0, "n_leaves"] == 0
        assert detail.loc[0, "n_clean_matches"] == 0
        assert detail.loc[0, "best_accession"] is None


class TestClassificationSchemaTolerance:
    def test_current_schema_yields_reads_and_classifier_names(self, tmp_path):
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100, 200])

        assert detail["in_classification"].tolist() == [True, True]
        assert detail["classifier_uniq_reads"].tolist() == [17.0, 5.0]
        assert detail.loc[0, "classifier_names"] == "Kraken2"
        assert detail.loc[1, "classifier_names"] == "Centrifuge;Kraken2"

    def test_legacy_schema_tolerates_different_column_names(self, tmp_path):
        root = _write_study(tmp_path, classification=LEGACY_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100, 200])

        assert detail["in_classification"].tolist() == [True, True]
        assert detail.loc[0, "classifier_names"] == "Kraken2"
        assert detail.loc[1, "classifier_names"] == "Centrifuge"

    def test_missing_read_column_yields_nan_not_zero(self, tmp_path):
        """No ``uniq_reads`` column must not masquerade as 'zero reads'."""
        root = _write_study(tmp_path, classification=LEGACY_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100])

        assert detail["classifier_uniq_reads"].isna().all()

    @pytest.mark.parametrize("coverage", [CURRENT_COVERAGE, LEGACY_COVERAGE])
    def test_both_coverage_schemas_detect_covered_assemblies(self, tmp_path, coverage):
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=coverage)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100, 200, 300])

        # A0 covered, A1 absent from coverage; A2 and A3 covered.
        assert detail["has_coverage_row"].tolist() == [True, True, True]

    def test_uncovered_assembly_is_flagged(self, tmp_path):
        """A taxid whose only assembly lacks a coverage row is not covered."""
        root = _write_study(
            tmp_path,
            classification=CURRENT_CLASSIFICATION,
            coverage="#rname\tcoverage\nA2\t5.0\n",
        )
        detail = load_reference_evidence(root, "ds1", M_STATS, [100, 200])

        assert detail["has_coverage_row"].tolist() == [False, True]

    def test_taxid_absent_from_classification(self, tmp_path):
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100, 999])

        assert detail["in_classification"].tolist() == [True, False]
        assert pd.isna(detail.loc[1, "classifier_uniq_reads"])

    def test_accepts_a_path_object_study_root(self, tmp_path):
        """``EvaluatorConfig.study_output_filepath`` is a ``Path``, not a ``str``.

        Passing a ``Path`` used to raise ``AttributeError`` inside the guarded
        detail builder, silently dropping the whole feature for every dataset.
        """
        _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)

        detail = load_reference_evidence(tmp_path, "ds1", M_STATS, [100])

        assert detail["in_classification"].tolist() == [True]


# ---------------------------------------------------------------------------
# loss_stage cascade
# ---------------------------------------------------------------------------


def _stage_for(**kwargs) -> str:
    """Build a one-row detail frame with overridden evidence fields."""
    row = {
        "taxid": 1,
        "in_classification": True,
        "classifier_reads_available": True,
        "classifier_uniq_reads": 5.0,
        "n_clean_matches": 1,
        "has_coverage_row": True,
        "numreads": 10.0,
        "n_leaves": 1,
        "best_accession": "A0",
        "is_trash_best": False,
        "leaf_rank": 0.0,
        "found_in_clade_post": True,
    }
    row.update(kwargs)
    return str(assign_loss_stages(pd.DataFrame([row]), keep_index=None)["loss_stage"].iloc[0])


def _flags_for(**kwargs):
    """Return the assigned flags for a one-row detail frame."""
    row = {
        "taxid": 1,
        "in_classification": True,
        "classifier_reads_available": True,
        "classifier_uniq_reads": 5.0,
        "n_clean_matches": 1,
        "has_coverage_row": True,
        "numreads": 10.0,
        "n_leaves": 1,
        "best_accession": "A0",
        "is_trash_best": False,
        "leaf_rank": 0.0,
        "found_in_clade_post": True,
    }
    row.update(kwargs)
    return assign_loss_stages(pd.DataFrame([row]), keep_index=None).iloc[0]


class TestLossStageCascade:
    def test_healthy_reference_survives(self):
        assert _stage_for() == "survived_post_cleanup"

    def test_absent_classification(self):
        assert _stage_for(in_classification=False, n_clean_matches=0, n_leaves=0, best_accession=None) == (
            "absent_classification"
        )

    def test_classified_no_assembly(self):
        """Listed by the classifier but with no reads and no assembly match."""
        assert _stage_for(n_clean_matches=0, n_leaves=0, best_accession=None, classifier_uniq_reads=0.0) == (
            "classified_no_assembly"
        )

    def test_classifier_credit_without_assembly_is_not_a_loss(self):
        """The bug this taxonomy exists for: `loss_stage` must agree with recall."""
        assert _stage_for(n_clean_matches=0, n_leaves=0, best_accession=None, classifier_uniq_reads=5.0) == (
            "survived_classifier_only"
        )

    def test_absent_classification_is_not_rescued_by_a_read_count(self):
        """No read count can credit a reference absent from the classification."""
        assert _stage_for(
            in_classification=False, n_clean_matches=0, n_leaves=0, best_accession=None, classifier_uniq_reads=5.0
        ) == "absent_classification"

    def test_legacy_classification_without_counts_credits_every_listed_taxid(self):
        """`load_classifier_detected_taxids` drops the >0 filter when absent."""
        assert _stage_for(
            n_clean_matches=0, n_leaves=0, best_accession=None, classifier_reads_available=False, classifier_uniq_reads=float("nan")
        ) == "survived_classifier_only"

    def test_assembly_zero_reads(self):
        assert _stage_for(numreads=0.0) == "assembly_zero_reads"

    def test_lost_clade_prediction(self):
        assert _stage_for(found_in_clade_post=False) == "lost_clade_prediction"

    def test_lost_recall_truncation(self):
        # truncation requires the leaf to exist and rank >= keep_index
        frame = pd.DataFrame(
            [
                {
                    "taxid": 1,
                    "in_classification": True,
                    "classifier_reads_available": True,
                    "classifier_uniq_reads": 5.0,
                    "n_clean_matches": 1,
                    "has_coverage_row": True,
                    "numreads": 10.0,
                    "n_leaves": 1,
                    "best_accession": "A0",
                    "is_trash_best": False,
                    "leaf_rank": 30.0,
                    "found_in_clade_post": False,
                }
            ]
        )
        assert str(assign_loss_stages(frame, keep_index=12)["loss_stage"].iloc[0]) == "lost_recall_truncation"

    def test_earliest_stage_wins_over_later_ones(self):
        """A reference can satisfy several conditions; the earliest must win."""
        stage = _stage_for(
            in_classification=False,
            n_clean_matches=0,
            numreads=0.0,
            is_trash_best=True,
            found_in_clade_post=False,
            best_accession=None,
            n_leaves=0,
        )
        assert stage == "absent_classification"

    def test_clade_loss_does_not_mask_an_earlier_failure(self):
        stage = _stage_for(
            n_clean_matches=0, is_trash_best=True, classifier_uniq_reads=0.0, found_in_clade_post=False
        )
        assert stage == "classified_no_assembly"

    def test_trash_best_match_is_reported_by_a_flag_not_a_dead_stage(self):
        """No cascade bucket can win once a clean match exists, so use the flag."""
        row = _flags_for(is_trash_best=True, n_clean_matches=0, classifier_uniq_reads=0.0)
        assert row["lost_best_match_trash"]
        assert row["lost_assembly_lookup"]

    def test_every_assigned_stage_is_in_the_documented_taxonomy(self):
        stages = {
            _stage_for(),
            _stage_for(in_classification=False, n_clean_matches=0, n_leaves=0, best_accession=None),
            _stage_for(n_clean_matches=0, n_leaves=0, best_accession=None, classifier_uniq_reads=0.0),
            _stage_for(n_clean_matches=0, n_leaves=0, best_accession=None),
            _stage_for(numreads=0.0),
            _stage_for(has_coverage_row=False),
            _stage_for(found_in_clade_post=False),
        }
        assert stages <= set(LOSS_STAGES)

    def test_diagnose_recall_gap_bucket_names_are_all_present(self):
        """Names must stay compatible with the diagnose_recall_gap script."""
        assert {
            "absent_classification",
            "classified_no_assembly",
            "assembly_no_coverage",
            "assembly_zero_reads",
        } <= set(LOSS_STAGES)

    def test_no_stage_is_unreachable_by_construction(self):
        """`no_best_match` / `best_match_trash` cannot win; they were removed."""
        assert "no_best_match" not in LOSS_STAGES
        assert "best_match_trash" not in LOSS_STAGES

    def test_only_two_stages_are_actual_recall_losses(self):
        """Everything else describes how far a *recalled* reference got."""
        assert LOSS_LOSSES == ("absent_classification", "classified_no_assembly")
        assert set(LOSS_LOSSES) < set(LOSS_STAGES)


# ---------------------------------------------------------------------------
# independent loss flags
# ---------------------------------------------------------------------------


class TestLossFlags:
    def test_flags_are_present_and_boolean(self, tmp_path):
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100, 200, 300, 999])
        apply_clade_membership(detail, CLADES, M_STATS, "post")
        out = assign_loss_stages(detail, keep_index=1)

        for flag in LOSS_FLAGS:
            assert flag in out.columns
            assert out[flag].dtype == bool

    def test_flags_are_independent_not_a_copy_of_the_cascade(self):
        """A trash leaf that is also truncated reports both flags, one stage."""
        detail = pd.DataFrame(
            [
                {
                    "taxid": 1,
                    "in_classification": True,
                    "classifier_reads_available": True,
                    "classifier_uniq_reads": 5.0,
                    "n_clean_matches": 1,
                    "has_coverage_row": True,
                    "numreads": 10.0,
                    "n_leaves": 1,
                    "best_accession": "A0",
                    "is_trash_best": True,
                    "leaf_rank": 30.0,
                    "found_in_clade_post": False,
                }
            ]
        )
        out = assign_loss_stages(detail, keep_index=12).iloc[0]

        assert out["loss_stage"] == "lost_recall_truncation"
        assert out["lost_recall_truncation"]
        assert out["lost_clade_prediction"]
        assert out["lost_best_match_trash"]
        assert out["recalled_classification_aware"]
        assert not out["survived_post_cleanup"]


# ---------------------------------------------------------------------------
# clade membership
# ---------------------------------------------------------------------------


class TestCladeMembership:
    @pytest.mark.parametrize("stage", ["pre", "post", "fixed"])
    def test_membership_marks_the_matching_taxids(self, stage):
        frame = pd.DataFrame(
            {
                "taxid": [100, 200, 300, 999],
                # L0 is taxid 100's best leaf, L3 is taxid 300's; both are in N0.
                "best_leaf": ["L0", "L1", "L3", None],
            }
        )
        apply_clade_membership(frame, CLADES, M_STATS, stage)

        assert frame[f"found_in_clade_{stage}"].tolist() == [True, False, True, False]
        assert frame.loc[0, f"clade_{stage}_node"] == "N0"
        assert frame.loc[2, f"clade_{stage}_node"] == "N0"
        assert frame.loc[0, f"clade_{stage}_precision"] == 0.5

    def test_membership_without_best_leaf_falls_back_to_taxid_map(self):
        """A bare ``taxid`` frame still gets membership, just no node/precision."""
        frame = pd.DataFrame({"taxid": [100, 200, 300]})
        apply_clade_membership(frame, CLADES, M_STATS, "post")

        assert frame["found_in_clade_post"].tolist() == [True, False, True]
        assert frame["clade_post_node"].isna().all()
        assert frame["clade_post_precision"].isna().all()

    def test_missing_clade_frame_marks_nothing(self):
        frame = pd.DataFrame({"taxid": [100, 300]})
        apply_clade_membership(frame, None, M_STATS, "post")

        assert not frame["found_in_clade_post"].any()

    def test_taxid_in_two_clades_reports_a_count(self):
        two = pd.DataFrame({"node": ["N0", "N1"], "leaves": [["L0"], ["L0"]], "node_precision": [0.9, 0.7]})
        frame = pd.DataFrame({"taxid": [100]})
        apply_clade_membership(frame, two, M_STATS, "post")

        assert frame.loc[0, "found_in_clade_post"]
        assert frame.loc[0, "clade_post_n_nodes"] == 2


# ---------------------------------------------------------------------------
# leaf table
# ---------------------------------------------------------------------------


class TestLeafTable:
    def test_taxid_is_the_input_taxid_not_the_assembly_taxid(self):
        """M_STATS already has its own ``taxid``; the leaf table must not reuse it."""
        leaves = build_reference_leaf_table("ds1", M_STATS, [100, 300])

        # A2 best-matches taxid 200, which is not an input, so it is excluded.
        assert leaves["taxid"].tolist() == [100, 100, 300]
        assert leaves["assembly_taxid"].tolist() == [10, 10, 30]

    def test_only_input_taxids_and_unassigned_leaves_are_kept(self):
        """Taxid 200 (a cross-hit target) is not an input, so it is dropped."""
        leaves = build_reference_leaf_table("ds1", M_STATS, [100])

        assert set(leaves["taxid"].dropna().astype(int)) == {100}

    def test_unassigned_leaves_are_retained_with_null_taxid(self):
        frame = M_STATS.copy()
        frame.loc["A2", "best_match_taxid"] = np.nan
        leaves = build_reference_leaf_table("ds1", frame, [100])

        assert leaves["taxid"].isna().sum() == 1
        assert (~leaves["is_input_taxid"]).sum() == 1

    def test_grain_is_unique_per_taxid_assembly(self):
        leaves = build_reference_leaf_table("ds1", M_STATS, [100, 200, 300])

        assert not leaves.duplicated(["data_set", "taxid", "assembly_accession"]).any()

    def test_leaf_rank_is_the_recall_filter_order(self):
        leaves = build_reference_leaf_table("ds1", M_STATS, [100, 200, 300])

        assert leaves["leaf_rank"].tolist() == sorted(leaves["leaf_rank"].tolist())

    def test_empty_m_stats_returns_schema_only_frame(self):
        leaves = build_reference_leaf_table("ds1", pd.DataFrame(), [100])

        assert leaves.empty
        assert "data_set" in leaves.columns
        assert "taxid" in leaves.columns


class TestFilterSurvival:
    def test_survival_requires_both_rank_and_coverage(self):
        """Both filters truncate by rank *and* drop zero-coverage leaves."""
        leaves = build_reference_leaf_table("ds1", M_STATS, [100])
        out = mark_filter_survival(leaves, keep_index=2, fixed_keep_index=2)

        by_leaf = out.set_index("leaf")["survives_recall_filter"]
        assert by_leaf["L0"]  # rank 0, coverage 11
        assert not by_leaf["L1"]  # rank 1 passes rank but coverage == 0

    def test_rank_beyond_keep_index_is_truncated(self):
        leaves = build_reference_leaf_table("ds1", M_STATS, [100])
        out = mark_filter_survival(leaves, keep_index=1, fixed_keep_index=1)

        assert not out.loc[out["leaf"] == "L1", "survives_recall_filter"].item()

    def test_none_keep_index_means_no_truncation(self):
        """A disabled filter must not report every leaf as lost.

        With ``keep_index=None`` the rank cutoff is inactive, so only the
        coverage condition applies: L0 (covered) survives, L1 (zero coverage)
        does not.
        """
        leaves = build_reference_leaf_table("ds1", M_STATS, [100])
        out = mark_filter_survival(leaves, keep_index=None, fixed_keep_index=None)

        by_leaf = out.set_index("leaf")["survives_recall_filter"]
        assert by_leaf["L0"]
        assert not by_leaf["L1"]


class TestLeafCladeMembership:
    def test_leaf_clade_flags_agree_with_the_wide_table(self):
        leaves = build_reference_leaf_table("ds1", M_STATS, [100, 200, 300])
        mark_clade_leaves(leaves, CLADES, M_STATS, "post")

        assert leaves.loc[leaves["leaf"].isin(["L0", "L3"]), "in_clade_post"].all()

    def test_wide_and_leaf_tables_agree_on_clade_membership(self):
        detail = pd.DataFrame({"taxid": [100, 200, 300]})
        apply_clade_membership(detail, CLADES, M_STATS, "post")
        leaves = build_reference_leaf_table("ds1", M_STATS, [100, 200, 300])
        mark_clade_leaves(leaves, CLADES, M_STATS, "post")

        wide = set(detail.loc[detail["found_in_clade_post"], "taxid"])
        long = set(leaves.loc[leaves["in_clade_post"], "taxid"].dropna().astype(int))
        assert wide == long


# ---------------------------------------------------------------------------
# schema / ordering guarantees
# ---------------------------------------------------------------------------


class TestSchemaGuarantees:
    def test_detail_columns_are_emitted_in_documented_order(self, tmp_path):
        root = _write_study(tmp_path, classification=CURRENT_CLASSIFICATION, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(root, "ds1", M_STATS, [100])
        ordered = order_reference_detail(detail)

        expected = [c for c in REFERENCE_DETAIL_COLUMNS if c in detail.columns]
        assert list(ordered.columns)[: len(expected)] == expected

    def test_legacy_columns_keep_their_original_names_and_order(self):
        """Existing readers select these by name and by leading position."""
        assert LEGACY_INPUT_DF_COLUMNS == (
            "sample",
            "taxid",
            "reads",
            "mutation_rate",
            "order",
            "family",
            "genus",
            "match_index",
            "match_coverage",
            "found_in_recall_filter",
            "found_in_fixed_filter",
            "data_set",
        )

    def test_evidence_load_does_not_mutate_the_input_m_stats(self):
        before = M_STATS.copy()
        build_reference_leaf_table("ds1", M_STATS, [100])

        pd.testing.assert_frame_equal(M_STATS, before)


# ---------------------------------------------------------------------------
# loss_stage <-> classification-aware recall
# ---------------------------------------------------------------------------


class TestLossStageMatchesRecall:
    """`loss_stage` and the headline recall must use one definition.

    Taxid 200 is the reference that used to be mislabelled: the classifier
    found it (``uniq_reads=5``) but its only leaf (A2) is trash, so there is no
    clean assembly match. It is recalled, therefore it must not be a loss stage.
    The invariant under test is
    ``recalled_classification_aware == loss_stage not in LOSS_LOSSES``.
    """

    def _detail(self, tmp_path, classification=CURRENT_CLASSIFICATION, keep_index=1):
        _write_study(tmp_path, classification=classification, coverage=CURRENT_COVERAGE)
        detail = load_reference_evidence(tmp_path, "ds1", M_STATS, [100, 200, 300, 999])
        apply_clade_membership(detail, CLADES, M_STATS, "post")
        return assign_loss_stages(detail, keep_index=keep_index).set_index("taxid")

    def test_classifier_found_reference_is_credited_not_lost(self, tmp_path):
        row = self._detail(tmp_path).loc[200]

        assert row["n_clean_matches"] == 0
        assert row["recalled_classification_aware"]
        assert row["recalled_classifier_only"]
        assert not row["recalled_assembly"]
        assert row["loss_stage"] == "survived_classifier_only"

    def test_reference_in_neither_source_is_still_a_loss(self, tmp_path):
        """Taxid 999 has no leaves and no classification entry."""
        row = self._detail(tmp_path).loc[999]

        assert not row["recalled_classification_aware"]
        assert row["loss_stage"] == "absent_classification"

    def test_invariant_holds_for_every_row(self, tmp_path):
        out = self._detail(tmp_path)

        assert set(out["loss_stage"]) <= set(LOSS_STAGES)
        assert (out["recalled_classification_aware"] != out["loss_stage"].isin(LOSS_LOSSES)).all()

    def test_recalled_set_equals_compute_recall_detected_set(self, tmp_path):
        """Compare against the real helper, not a re-implementation of it."""
        out = self._detail(tmp_path)

        expected = TestMatchSemantics._recall_matched_taxids(M_STATS) | load_classifier_detected_taxids(
            str(tmp_path), "ds1"
        )
        got = set(out.loc[out["recalled_classification_aware"]].index.astype(int))

        assert got == expected & set(out.index.astype(int))

    def test_legacy_classification_shapes_agree_with_the_helper(self, tmp_path):
        """No `uniq_reads` column => every listed taxid counts as detected."""
        out = self._detail(tmp_path, classification=LEGACY_CLASSIFICATION)

        expected = load_classifier_detected_taxids(str(tmp_path), "ds1")
        got = set(out.loc[out["recalled_classification_aware"]].index.astype(int))

        assert expected == {100, 200}
        assert got == expected | TestMatchSemantics._recall_matched_taxids(M_STATS) & set(out.index.astype(int))
        # 100 also has a clean assembly match, so it is assembly-recalled.
        assert out.loc[100, "recalled_assembly"]
        assert not out.loc[100, "recalled_classifier_only"]
        assert out.loc[200, "loss_stage"] == "survived_classifier_only"

    def test_recall_counts_do_not_double_count(self, tmp_path):
        out = self._detail(tmp_path)

        assert out["recalled_classification_aware"].sum() == out["recalled_assembly"].sum() + out[
            "recalled_classifier_only"
        ].sum()
        assert len(out) == 4
