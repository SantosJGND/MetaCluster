"""
Per-reference downstream tracking for the evaluator module.

Builds two complementary tables that make every input reference traceable
through the whole downstream clustering pipeline:

1. **Wide per-reference table** (:func:`load_reference_evidence`) — one row per
   input taxid, carrying classifier evidence, mapping statistics, per-stage
   membership, and an ordered ``loss_stage`` attribution.
2. **Long per-leaf table** (:func:`build_reference_leaf_table`) — one row per
   ``(taxid, assembly_accession)`` leaf, exposing the multiple leaves and
   non-best leaves that the wide table necessarily aggregates away.

The ``loss_stage`` taxonomy intentionally reuses the bucket names of
``analysis_scripts/diagnose_recall_gap.py`` for the four upstream buckets, so
the two producers can be cross-checked against each other.
"""

import logging
import os
from typing import Any

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ``loss_stage`` is a first-match cascade whose terminal success condition is the
# **classification-aware** recall definition used for the headline
# ``recall_baseline``: a reference is recalled if it either has a clean
# assembly match (``n_clean_matches > 0``) or was called by the classifier.
#
# Hence the invariant
#
#     recalled_classification_aware == loss_stage NOT IN
#         {absent_classification, classified_no_assembly}
#
# Only those two buckets are genuine recall losses. Every other bucket is a
# *degradation* state describing how far a recalled reference got.
#
# ``absent_classification``, ``classified_no_assembly``, ``assembly_no_coverage``
# and ``assembly_zero_reads`` keep the names used by
# ``analysis_scripts/diagnose_recall_gap.py``.
#
# ``no_best_match`` / ``best_match_trash`` are deliberately absent: once the
# classifier is credited, a reference whose only leaves are trash is recalled
# (``survived_classifier_only``), and once there is a clean match there is by
# definition a best match — so neither can ever win the cascade. Both conditions
# are still reported, as the ``n_clean_matches`` / ``is_trash_best`` columns and
# the ``lost_no_best_match`` / ``lost_best_match_trash`` flags.
LOSS_STAGES = (
    "absent_classification",
    "survived_classifier_only",
    "classified_no_assembly",
    "assembly_no_coverage",
    "assembly_zero_reads",
    "lost_recall_truncation",
    "lost_clade_prediction",
    "survived_post_cleanup",
)

# Independent boolean flags, one per condition. Their sum per dataset is
# expected to exceed the count of `loss_stage` rows, because the conditions
# overlap (e.g. a reference can be both zero-reads and beyond keep_index).
# The only stages that count against classification-aware recall. Everything
# else in LOSS_STAGES is a degradation state on an already-recalled reference.
LOSS_LOSSES = ("absent_classification", "classified_no_assembly")

LOSS_FLAGS = (
    # Classification-aware recall decomposition (mutually consistent with
    # `loss_stage`, see the invariant above).
    "recalled_classification_aware",
    "recalled_assembly",
    "recalled_classifier_only",
    # Per-condition diagnostics.
    "lost_before_classification",
    "lost_assembly_lookup",
    "lost_no_coverage",
    "lost_zero_reads",
    "lost_no_best_match",
    "lost_best_match_trash",
    "lost_recall_truncation",
    "lost_clade_prediction",
    "survived_post_cleanup",
)

# Mapping/aggregate columns taken from the taxid's best-match leaf.
MAPPING_COLUMNS = (
    "best_accession",
    "best_leaf",
    "leaf_rank",
    "numreads",
    "total_uniq_reads",
    "covbases",
    "meanmapq",
    "error_rate",
    "best_match_score",
    "best_match_level",
    "max_coverage",
    "max_shared",
    "n_leaves",
    "n_clean_matches",
    "n_covered_leaves",
    "is_trash_best",
)

EVIDENCE_COLUMNS = (
    "in_classification",
    "classifier_reads_available",
    "classifier_uniq_reads",
    "classifier_names",
    "has_coverage_row",
)

CLADE_MEMBERSHIP_COLUMNS = (
    "found_in_clade_pre",
    "found_in_clade_post",
    "found_in_clade_fixed",
    "clade_pre_node",
    "clade_post_node",
    "clade_fixed_node",
    "clade_pre_precision",
    "clade_post_precision",
    "clade_fixed_precision",
    "clade_pre_n_nodes",
    "clade_post_n_nodes",
    "clade_fixed_n_nodes",
)

REFERENCE_DETAIL_COLUMNS = EVIDENCE_COLUMNS + MAPPING_COLUMNS + CLADE_MEMBERSHIP_COLUMNS + ("loss_stage",) + LOSS_FLAGS

# The pre-existing ``test_datasets_input_df.tsv`` schema, in its original order.
# New columns are appended *after* this block so downstream readers that select
# by name or by leading position keep working unchanged.
LEGACY_INPUT_DF_COLUMNS = (
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


LEAF_COLUMNS = (
    # `taxid` is the *input* taxid this leaf was matched to (i.e. the wide
    # table's join key), not m_stats' own `taxid` column, which holds the
    # matched assembly's taxid and is kept separately as `assembly_taxid`.
    "taxid",
    "assembly_taxid",
    "assembly_accession",
    "leaf",
    "leaf_rank",
    "best_match_is_best",
    "best_match_score",
    "best_match_level",
    "coverage",
    "covbases",
    "meanmapq",
    "numreads",
    "error_rate",
    "total_uniq_reads",
    "is_trash",
    "max_shared",
    "is_input_taxid",
    "survives_recall_filter",
    "survives_fixed_filter",
    "in_clade_pre",
    "in_clade_post",
    "in_clade_fixed",
)


def _read_tsv(path: str) -> pd.DataFrame | None:
    """Read a TSV, returning None when the file is absent or unreadable."""
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception as e:  # pragma: no cover - defensive
        logger.warning(f"Could not read {path}: {e}")
        return None


# Placeholders that stand for "this classifier did not fire" in the legacy
# one-column-per-classifier schema.
_CLASSIFIER_PLACEHOLDERS = {"", "-", "--", "0", "0.0", "false", "nan", "none", "null", "na"}


def _classifier_name_columns(df: pd.DataFrame) -> list[str]:
    """Return the columns that hold classifier names, tolerating both schemas.

    Two classification schemas exist in the wild:
      - ``taxid, description, uniq_reads, classifiers`` (current studies)
      - ``description, taxid, classifier_centrifuge, classifier_kraken2,
        classification`` (older raw virus data)

    In the legacy schema the classifiers are spread across one column each, and
    ``classification`` is a 0/1 flag rather than a classifier name — so only the
    ``classifier_*`` columns are used.
    """
    if "classifiers" in df.columns:
        return ["classifiers"]
    return [c for c in df.columns if c.startswith("classifier_")]


def _classifier_names_by_taxid(df: pd.DataFrame) -> pd.Series | None:
    """Map taxid -> ``"A;B"`` of the classifiers that reported it.

    Returns None when the source table carries no classifier columns at all.
    Values that are just placeholders are dropped, so a taxid only ever lists
    classifiers that actually fired.
    """
    columns = _classifier_name_columns(df)
    if not columns:
        return None

    frames = []
    for col in columns:
        names = df[col].fillna("").astype(str).str.strip()
        frames.append(df.assign(_n=names).query("_n.str.lower() not in @_CLASSIFIER_PLACEHOLDERS"))
    stacked = pd.concat(frames, ignore_index=True)
    if stacked.empty:
        return None

    names = stacked.groupby("taxid")["_n"].apply(lambda s: ";".join(sorted(set(s))))
    return names[names.str.len() > 0]


def _uniq_reads_series(df: pd.DataFrame) -> pd.Series | None:
    """Return a taxid -> uniq_reads mapping, or None when unavailable.

    Returns None (rather than a zero-filled series) when the source table has
    no read-count column, so callers can emit NaN instead of a misleading 0.
    """
    for candidate in ("uniq_reads", "total_uniq_reads"):
        if candidate in df.columns:
            values = pd.to_numeric(df[candidate], errors="coerce")
            return df.assign(_uniq=values).groupby("taxid")["_uniq"].max()
    return None


def load_reference_evidence(
    study_output_filepath: str | os.PathLike,
    data_set_name: str,
    m_stats: pd.DataFrame,
    input_taxids: list[int],
    preloaded: dict[str, pd.DataFrame | None] | None = None,
) -> pd.DataFrame:
    """Build the wide per-reference evidence table for one dataset.

    Args:
        study_output_filepath: Root of the study output directory. Accepts
            ``str`` or ``Path`` (``EvaluatorConfig.study_output_filepath`` is a
            ``Path``).
        data_set_name: Dataset folder name.
        m_stats: m_stats matrix for the dataset, indexed by ``assid``
            (``filter_no_leaf=False``) so that the ``leaf`` column survives.
        input_taxids: The dataset's input taxids; one output row per entry.
        preloaded: Optional mapping of ``"matched"`` / ``"classification"`` /
            ``"coverage"`` to already-read DataFrames, to avoid re-reading
            files that ``get_m_stats_matrix`` has already loaded.

    Returns:
        DataFrame with one row per input taxid and the columns listed in
        :data:`REFERENCE_DETAIL_COLUMNS`. ``loss_stage`` and the
        :data:`LOSS_FLAGS` are left unset — see :func:`assign_loss_stages`.
    """
    preloaded = preloaded or {}
    base = os.fspath(study_output_filepath).rstrip("/")
    paths = {
        "classification": f"{base}/{data_set_name}/classification/{data_set_name}_merged_classification.tsv",
        "coverage": f"{base}/{data_set_name}/output/merged_coverage_statistics.tsv",
    }
    classification = preloaded.get("classification")
    coverage = preloaded.get("coverage")
    if classification is None and "classification" not in preloaded:
        classification = _read_tsv(paths["classification"])
    if coverage is None and "coverage" not in preloaded:
        coverage = _read_tsv(paths["coverage"])

    detail = pd.DataFrame({"taxid": pd.Series(list(dict.fromkeys(input_taxids)), dtype="object")})
    detail["taxid"] = pd.to_numeric(detail["taxid"], errors="coerce")

    detail["in_classification"] = False
    detail["classifier_reads_available"] = False
    detail["classifier_uniq_reads"] = np.nan
    detail["classifier_names"] = None
    detail["has_coverage_row"] = False

    if classification is not None and "taxid" in classification.columns:
        classification = classification.copy()
        classification["taxid"] = pd.to_numeric(classification["taxid"], errors="coerce")
        classification = classification.dropna(subset=["taxid"])
        detail["in_classification"] = detail["taxid"].isin(set(classification["taxid"]))

        uniq = _uniq_reads_series(classification)
        if uniq is not None:
            # Only then does ``uniq_reads > 0`` mean "detected"; see
            # ``dataset_processor.load_classifier_detected_taxids``, which credits
            # every listed taxid when the column is absent.
            detail["classifier_reads_available"] = True
            detail["classifier_uniq_reads"] = detail["taxid"].map(uniq)

        names = _classifier_names_by_taxid(classification)
        if names is not None:
            detail["classifier_names"] = detail["taxid"].map(names)

    # Accessions that actually produced a coverage record, i.e. were mapped.
    covered_accs: set[str] = set()
    if coverage is not None:
        rname = "#rname" if "#rname" in coverage.columns else ("index" if "index" in coverage.columns else None)
        if rname is not None:
            covered_accs = set(coverage[rname].dropna().astype(str))

    _attach_mapping_columns(detail, m_stats, covered_accs)

    detail["has_coverage_row"] = detail["n_covered_leaves"].fillna(0).astype(int) > 0

    for col in MAPPING_COLUMNS:
        if col not in detail.columns:
            detail[col] = np.nan
    for col in ("n_leaves", "n_clean_matches", "n_covered_leaves"):
        detail[col] = detail[col].fillna(0).astype(int)
    for flag in ("is_trash_best",):
        detail[flag] = detail[flag].fillna(False).astype(bool)

    return detail


def _attach_mapping_columns(detail: pd.DataFrame, m_stats: pd.DataFrame, covered_accs: set[str]) -> None:
    """Attach best-match leaf statistics and per-taxid leaf aggregates.

    Everything here keys on ``m_stats["best_match_taxid"]`` — the reference that
    a leaf was matched to — because that is what ``metrics.compute_recall`` and
    ``compute_clade_recall`` use. ``matched_assemblies.tsv``'s own ``taxid``
    column is the *matched assembly's* taxid, so joining the input taxids on it
    would report that no reference was ever matched.

    ``n_clean_matches`` mirrors the pipeline's recall criterion exactly: a
    reference counts as assembly-recalled when it has at least one leaf that is
    neither trash nor non-best (``is_trash == False and best_match_is_best == True``).
    """
    defaults: dict[str, Any] = {
        "best_accession": pd.Series([None] * len(detail), dtype="object"),
        "best_leaf": pd.Series([None] * len(detail), dtype="object"),
        "best_match_level": pd.Series([None] * len(detail), dtype="object"),
        "leaf_rank": np.nan,
        "numreads": np.nan,
        "total_uniq_reads": np.nan,
        "covbases": np.nan,
        "meanmapq": np.nan,
        "error_rate": np.nan,
        "best_match_score": np.nan,
        "max_coverage": np.nan,
        "max_shared": np.nan,
        "n_leaves": 0,
        "n_clean_matches": 0,
        "n_covered_leaves": 0,
        "is_trash_best": False,
    }
    assert set(defaults) == set(MAPPING_COLUMNS), "MAPPING_COLUMNS and defaults are out of sync"
    for col, default in defaults.items():
        detail[col] = default

    if m_stats is None or m_stats.empty or "best_match_taxid" not in m_stats.columns:
        return

    frame = m_stats.copy()
    frame["best_match_taxid"] = pd.to_numeric(frame["best_match_taxid"], errors="coerce")
    # m_stats is sorted by total_uniq_reads descending (node_stats._finalize_m_stats),
    # so the row position is the rank the recall filter truncates on.
    frame["leaf_rank"] = np.arange(len(frame))
    if "leaf" not in frame.columns:
        frame["leaf"] = frame.index.astype(str)
    frame["_accession"] = frame.index.astype(str)

    detail_taxids = set(detail["taxid"])

    # ``is_trash``/``best_match_is_best`` default to True when absent so that an
    # unexpected m_stats shape is treated as "not a clean match" rather than
    # silently counting as recalled.
    is_trash = frame.get("is_trash", pd.Series(True, index=frame.index)).fillna(True).astype(bool)
    is_best = frame.get("best_match_is_best", pd.Series(False, index=frame.index)).fillna(False).astype(bool)

    frame = frame.assign(_clean=(~is_trash) & is_best, _covered=frame["_accession"].isin(covered_accs))

    for taxid, group in frame.groupby("best_match_taxid", dropna=True):
        if taxid not in detail_taxids:
            continue
        mask = detail["taxid"] == taxid
        detail.loc[mask, "n_leaves"] = len(group)
        detail.loc[mask, "n_clean_matches"] = int(group["_clean"].sum())
        detail.loc[mask, "n_covered_leaves"] = int(group["_covered"].sum())
        if "coverage" in group.columns:
            detail.loc[mask, "max_coverage"] = float(pd.to_numeric(group["coverage"], errors="coerce").max())

        # The "best" leaf is the highest-ranked clean one, which is the row the
        # pipeline's recall and the leaf table both treat as the match.
        best = group[group["_clean"]]
        if best.empty:
            best = group
        if best.empty:
            continue
        row = best.iloc[0]
        detail.loc[mask, "best_accession"] = row["_accession"]
        detail.loc[mask, "best_leaf"] = row["leaf"]
        detail.loc[mask, "leaf_rank"] = int(row["leaf_rank"])
        for col in (
            "numreads",
            "total_uniq_reads",
            "covbases",
            "meanmapq",
            "error_rate",
            "best_match_score",
            "best_match_level",
            "max_shared",
        ):
            if col in group.columns:
                detail.loc[mask, col] = row[col]
        if "is_trash" in group.columns:
            detail.loc[mask, "is_trash_best"] = bool(row["is_trash"]) if pd.notna(row["is_trash"]) else False


def build_reference_leaf_table(
    data_set_name: str,
    m_stats: pd.DataFrame,
    input_taxids: list[int],
) -> pd.DataFrame:
    """Build the long per-leaf table for one dataset.

    One row per m_stats leaf whose ``best_match_taxid`` is an input taxid, plus
    leaves carrying a null ``best_match_taxid`` (unassigned/trash) so the
    denominator of the clustering output is visible.
    """
    columns = ["data_set"] + list(LEAF_COLUMNS)
    if m_stats is None or m_stats.empty:
        return pd.DataFrame(columns=columns)

    frame = m_stats.copy()
    frame["best_match_taxid"] = pd.to_numeric(frame["best_match_taxid"], errors="coerce")
    frame["leaf_rank"] = np.arange(len(frame))
    frame["assembly_accession"] = frame.index.astype(str)
    if "leaf" not in frame.columns:
        frame["leaf"] = frame.index.astype(str)

    input_set = {int(t) for t in input_taxids if pd.notna(t)}
    is_input = frame["best_match_taxid"].isin(input_set)
    keep = is_input | frame["best_match_taxid"].isna()
    frame = frame[keep]

    # m_stats already carries a `taxid` column holding the *assembly's* taxid;
    # rename it so the leaf table's `taxid` unambiguously means the input taxid.
    if "taxid" in frame.columns:
        frame = frame.rename(columns={"taxid": "assembly_taxid"})
    frame["taxid"] = frame["best_match_taxid"]

    out = frame.reindex(columns=[c for c in LEAF_COLUMNS if c in frame.columns])
    out.insert(0, "data_set", data_set_name)
    out["is_input_taxid"] = is_input[keep].values

    for col in ("survives_recall_filter", "survives_fixed_filter", "in_clade_pre", "in_clade_post", "in_clade_fixed"):
        out[col] = False

    for col in LEAF_COLUMNS:
        if col not in out.columns:
            out[col] = np.nan
    return out[columns]


def _clade_maps(results_df: pd.DataFrame | None, m_stats: pd.DataFrame) -> tuple[dict, dict]:
    """Map leaf -> (node, node_precision) and taxid -> set of nodes.

    A reference with several leaves can land in more than one predicted clade,
    so the taxid-level mapping is a set and the count is reported alongside.
    """
    leaf_to_clade: dict = {}
    taxid_to_nodes: dict = {}
    if results_df is None or results_df.empty:
        return leaf_to_clade, taxid_to_nodes

    best_taxid_by_leaf: dict = {}
    if m_stats is not None and not m_stats.empty and "leaf" in m_stats.columns:
        for leaf, bmt in zip(m_stats["leaf"], m_stats["best_match_taxid"]):
            if pd.notna(bmt):
                best_taxid_by_leaf[str(leaf)] = bmt

    for _, row in results_df.iterrows():
        node = row.get("node")
        precision = row.get("node_precision")
        leaves = row.get("leaves") or []
        if not isinstance(leaves, (list, tuple, set, np.ndarray)):
            leaves = [leaves]
        for leaf in leaves:
            leaf_to_clade.setdefault(str(leaf), (node, precision))
            taxid = best_taxid_by_leaf.get(str(leaf))
            if taxid is not None and pd.notna(taxid):
                taxid_to_nodes.setdefault(taxid, set()).add(node)

    return leaf_to_clade, taxid_to_nodes


def apply_clade_membership(
    detail: pd.DataFrame,
    results_df: pd.DataFrame | None,
    m_stats: pd.DataFrame,
    prefix: str,
) -> None:
    """Fill ``found_in_clade_*`` / ``clade_*_node`` / precision columns.

    Safe to call on a frame that has no ``best_leaf`` column (e.g. a bare
    ``taxid`` frame): membership is then derived from the taxid->nodes map only,
    with no per-row node or precision to report.
    """
    leaf_to_clade, taxid_to_nodes = _clade_maps(results_df, m_stats)
    found = f"found_in_clade_{prefix}"
    detail[found] = False
    detail[f"clade_{prefix}_node"] = np.nan
    detail[f"clade_{prefix}_precision"] = np.nan
    detail[f"clade_{prefix}_n_nodes"] = 0

    if not taxid_to_nodes:
        return

    n_nodes = detail["taxid"].map(lambda t: len(taxid_to_nodes.get(t, ())))
    detail[found] = n_nodes > 0
    detail[f"clade_{prefix}_n_nodes"] = n_nodes.astype(int)

    if "best_leaf" in detail.columns:
        best_leaf = detail["best_leaf"].astype("object")
        detail[f"clade_{prefix}_node"] = best_leaf.map(
            lambda leaf: leaf_to_clade.get(str(leaf), (None, np.nan))[0] if pd.notna(leaf) else None
        )
        detail[f"clade_{prefix}_precision"] = best_leaf.map(
            lambda leaf: leaf_to_clade.get(str(leaf), (None, np.nan))[1] if pd.notna(leaf) else np.nan
        )


def mark_clade_leaves(
    leaves: pd.DataFrame, results_df: pd.DataFrame | None, m_stats: pd.DataFrame, prefix: str
) -> None:
    """Set ``in_clade_*`` on the long leaf table."""
    col = f"in_clade_{prefix}"
    leaves[col] = False
    if results_df is None or results_df.empty or "leaves" not in results_df.columns:
        return
    if "leaf" not in leaves.columns:
        return
    members: set = set()
    for leaves_value in results_df["leaves"]:
        if isinstance(leaves_value, (list, tuple, set, np.ndarray)):
            members.update(str(x) for x in leaves_value)
        elif pd.notna(leaves_value):
            members.add(str(leaves_value))
    leaves[col] = leaves["leaf"].astype(str).isin(members)


def assign_loss_stages(
    detail: pd.DataFrame,
    keep_index: int | None = None,
    recall_filter_applied: bool = True,
) -> pd.DataFrame:
    """Assign the ordered ``loss_stage`` plus independent loss flags.

    ``loss_stage`` is a strict first-match cascade in pipeline order: the
    earliest stage that claims a reference wins, so the stage counts form a
    proper funnel. The boolean flags are deliberately *independent* — a
    reference can fail several conditions at once — so the flags can be used to
    verify the cascade rather than being a second copy of it.

    The first four stage names match
    ``analysis_scripts/diagnose_recall_gap.py`` exactly.
    """
    in_classification = detail["in_classification"] == True  # noqa: E712
    no_coverage_row = detail["has_coverage_row"] == False  # noqa: E712
    zero_reads = detail["numreads"].fillna(0).astype(float) == 0

    has_leaves = detail["n_leaves"].fillna(0).astype(int) > 0
    has_best = has_leaves & detail["best_accession"].notna()
    best_trash = has_best & (detail["is_trash_best"] == True)  # noqa: E712

    found_post = detail["found_in_clade_post"] == True  # noqa: E712

    # Classification-aware recall, matching the headline `recall_baseline`:
    # a reference counts as recalled with a clean assembly match **or** a
    # classifier call. `load_classifier_detected_taxids` credits every listed
    # taxid when the classification table has no read-count column, so the
    # `uniq_reads > 0` test only applies when a count is actually available.
    recalled_assembly = detail["n_clean_matches"].fillna(0).astype(int) > 0
    reads_available = detail["classifier_reads_available"] == True  # noqa: E712
    reads = pd.to_numeric(detail["classifier_uniq_reads"], errors="coerce")
    recalled_classifier = in_classification & (~reads_available | (reads > 0))
    recalled_classifier_only = recalled_classifier & ~recalled_assembly
    recalled = recalled_assembly | recalled_classifier

    truncated = pd.Series(False, index=detail.index)
    if recall_filter_applied and keep_index is not None:
        rank = pd.to_numeric(detail["leaf_rank"], errors="coerce")
        truncated = has_best & rank.notna() & (rank >= float(keep_index))

    # First match in pipeline order wins. The two `~recalled` gates make the
    # cascade agree with the classification-aware recall definition: a reference
    # the classifier found is never labelled a recall loss, it is credited.
    conditions = [
        ("absent_classification", ~recalled & ~in_classification),
        ("survived_classifier_only", recalled_classifier_only),
        ("classified_no_assembly", ~recalled & in_classification),
        ("assembly_no_coverage", has_best & no_coverage_row),
        ("assembly_zero_reads", has_best & ~no_coverage_row & zero_reads),
        ("lost_recall_truncation", truncated),
        ("lost_clade_prediction", ~found_post),
    ]

    stage = pd.Series(pd.NA, index=detail.index, dtype="object")
    for name, applies in conditions:
        # Fill only where still unclaimed, so earlier stages take precedence.
        stage = stage.mask(stage.isna() & applies, name)
    stage = stage.fillna("survived_post_cleanup")

    detail["loss_stage"] = pd.Categorical(stage, categories=list(LOSS_STAGES), ordered=True)

    flags = {
        "recalled_classification_aware": recalled,
        "recalled_assembly": recalled_assembly,
        "recalled_classifier_only": recalled_classifier_only,
        "lost_before_classification": ~in_classification,
        "lost_assembly_lookup": ~recalled_assembly,
        "lost_no_coverage": has_best & no_coverage_row,
        "lost_zero_reads": has_best & ~no_coverage_row & zero_reads,
        "lost_no_best_match": has_leaves & ~has_best,
        "lost_best_match_trash": best_trash,
        "lost_recall_truncation": truncated,
        "lost_clade_prediction": ~found_post,
        "survived_post_cleanup": found_post,
    }
    assert set(flags) == set(LOSS_FLAGS), "LOSS_FLAGS and flag assignments are out of sync"
    for flag, values in flags.items():
        detail[flag] = values.fillna(False).astype(bool)

    return detail


def mark_filter_survival(
    leaves: pd.DataFrame,
    keep_index: int | None,
    fixed_keep_index: int | None,
) -> pd.DataFrame:
    """Flag which leaves survive the recall and fixed filters.

    Both filters truncate on the same ``total_uniq_reads``-descending rank and
    additionally require ``coverage > 0`` (``dataset_processor._apply_*_filter``),
    so both conditions are applied here.

    A ``None`` keep index means that filter did not run (``--no-apply-recall-filter``),
    so only the coverage condition applies. It does *not* mean "nothing survived" —
    that would report a disabled filter as a total loss.
    """
    if leaves.empty:
        return leaves

    covered = pd.to_numeric(leaves.get("coverage"), errors="coerce") > 0
    rank = pd.to_numeric(leaves["leaf_rank"], errors="coerce")

    def survives(keep: int | None) -> pd.Series:
        if keep is None:
            return covered
        return covered & rank.lt(float(keep))

    leaves["survives_recall_filter"] = survives(keep_index)
    leaves["survives_fixed_filter"] = survives(fixed_keep_index)
    return leaves


def order_reference_detail(detail: pd.DataFrame) -> pd.DataFrame:
    """Return the wide table with :data:`REFERENCE_DETAIL_COLUMNS` in a fixed order.

    Any extra columns are appended after the known ones rather than dropped, so
    new fields can be added without silently disappearing from the TSV.
    """
    ordered = [c for c in REFERENCE_DETAIL_COLUMNS if c in detail.columns]
    remaining = [c for c in detail.columns if c not in ordered]
    return detail[ordered + remaining]
