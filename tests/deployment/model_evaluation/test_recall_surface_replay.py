"""Tests: recall probability surfaces can be replayed from the saved
explanatory/recall_model directory without refitting (analysis_data_extractor.py).

Requires statsmodels + patsy, which only exist in the main .venv (the extractor
runs there). Skipped wherever statsmodels is unavailable.
"""

import os

import numpy as np
import pandas as pd
import pytest

from deployment.model_evaluation.analysis_data_extractor import (
    RECALL_FORMULAS,
    _compute_surface_predictions,
    _recall_surface_grid,
    _render_surface,
    _replay_recall_surface,
    fit_recall_gee,
    plot_recall_surface,
    plot_recall_surface_from_saved,
    save_recall_summary,
)

statsmodels = pytest.importorskip("statsmodels")  # noqa: F841


@pytest.fixture(scope="module")
def recall_df():
    rng = np.random.RandomState(42)
    orders = ["Alphacoronavirus", "Betacoronavirus", "Gammacoronavirus", "Deltacoronavirus"]
    rows = []
    for ds in range(6):
        for ti in range(8):
            order = orders[(ds + ti) % len(orders)]
            reads = int(rng.randint(20, 20000))
            mutation_rate = float(rng.uniform(0.0, 0.5))
            recalled = int((rng.rand() + 0.3 * (mutation_rate + reads / 40000.0)) > 0.4)
            rows.append(
                {
                    "data_set": f"ds{ds}",
                    "taxid": 10_000 + ti,
                    "reads_simulated": reads,
                    "mutation_rate": mutation_rate,
                    "order": order,
                    "recalled": recalled,
                }
            )
    return pd.DataFrame(rows)


@pytest.fixture(scope="module")
def fitted_results(recall_df):
    results = {}
    for label, formula in RECALL_FORMULAS.items():
        result = fit_recall_gee(recall_df, formula)
        assert result is not None, f"GEE variant {label} failed to converge"
        results[label] = result
    return results


@pytest.fixture()
def saved_recall_dir(recall_df, fitted_results, tmp_path_factory):
    """Serialize recall_data.tsv + variant summaries exactly like the pipeline."""
    recall_dir = tmp_path_factory.mktemp("recall_model")
    recall_df.to_csv(recall_dir / "recall_data.tsv", sep="\t", index=False)
    for label, result in fitted_results.items():
        save_recall_summary(result, label, str(recall_dir))
    return recall_dir


def _grid(recall_df):
    recall_df = recall_df.copy()
    recall_df["log_reads"] = np.log1p(recall_df["reads_simulated"])
    return _recall_surface_grid(recall_df)


def test_replay_matches_fitted_surfaces(fitted_results, recall_df, saved_recall_dir, tmp_path):
    grid_df, mr_grid, lr_grid, order_levels = _grid(recall_df)
    recall_df_logged = recall_df.copy()
    recall_df_logged["log_reads"] = np.log1p(recall_df_logged["reads_simulated"])

    replay_out = tmp_path / "replay"
    plot_recall_surface(fitted_results, str(tmp_path / "live"), order_mode="reference")
    plot_recall_surface_from_saved(str(saved_recall_dir), output_dir=str(replay_out), order_mode="reference")

    for label, result in fitted_results.items():
        oracle, _ = _compute_surface_predictions(result.predict, grid_df, mr_grid, lr_grid, order_levels, "reference")
        replay, _ = _replay_recall_surface(
            recall_df_logged, str(saved_recall_dir), label, "reference", grid_df, mr_grid, lr_grid, order_levels
        )
        assert np.allclose(replay, oracle, atol=1e-9), f"variant {label} replay != fitted surface"

    assert (replay_out / "recall_probability_surface.png").exists()


def test_replay_average_matches_fitted_average(fitted_results, recall_df, saved_recall_dir):
    grid_df, mr_grid, lr_grid, order_levels = _grid(recall_df)
    recall_df_logged = recall_df.copy()
    recall_df_logged["log_reads"] = np.log1p(recall_df_logged["reads_simulated"])

    label = "A"
    oracle, _ = _compute_surface_predictions(
        fitted_results[label].predict, grid_df, mr_grid, lr_grid, order_levels, "average"
    )
    replay, _ = _replay_recall_surface(
        recall_df_logged, str(saved_recall_dir), label, "average", grid_df, mr_grid, lr_grid, order_levels
    )
    assert np.allclose(replay, oracle, atol=1e-9)


def test_variants_produce_distinct_surfaces(fitted_results, recall_df):
    grid_df, mr_grid, lr_grid, order_levels = _grid(recall_df)
    probs = {}
    for label in ["A", "B", "C"]:
        probs[label], _ = _compute_surface_predictions(
            fitted_results[label].predict, grid_df, mr_grid, lr_grid, order_levels, "reference"
        )
    assert not np.allclose(probs["A"], probs["B"], atol=1e-9)
    assert not np.allclose(probs["A"], probs["C"], atol=1e-9)


def test_average_order_mode_spans_orders_and_differs(fitted_results, recall_df):
    grid_df, mr_grid, lr_grid, order_levels = _grid(recall_df)
    assert len(order_levels) > 1

    label = "A"
    ref, _ = _compute_surface_predictions(
        fitted_results[label].predict, grid_df, mr_grid, lr_grid, order_levels, "reference"
    )
    avg, _ = _compute_surface_predictions(
        fitted_results[label].predict, grid_df, mr_grid, lr_grid, order_levels, "average"
    )
    assert not np.allclose(ref, avg, atol=1e-9)

    per_order = np.stack(
        [
            np.reshape(np.asarray(fitted_results[label].predict(grid_df.assign(order=o))), (len(mr_grid), len(mr_grid)))
            for o in order_levels
        ]
    )
    assert np.allclose(avg, per_order.mean(axis=0), atol=1e-9)


def test_render_is_deterministic(tmp_path):
    from deployment.model_evaluation.analysis_data_extractor import GRID_SIZE

    rng = np.random.RandomState(0)
    grid = np.sort(1 / (1 + np.exp(-rng.randn(GRID_SIZE, GRID_SIZE))))
    prob_by_label = {"A": (grid, "order=X"), "B": (grid[::-1], "order=Y")}
    out1 = tmp_path / "one"
    out2 = tmp_path / "two"
    out1.mkdir()
    out2.mkdir()
    _render_surface(prob_by_label, np.linspace(0, 10, GRID_SIZE), np.linspace(0, 1, GRID_SIZE), str(out1))
    _render_surface(prob_by_label, np.linspace(0, 10, GRID_SIZE), np.linspace(0, 1, GRID_SIZE), str(out2))
    assert (out1 / "recall_probability_surface.png").read_bytes() == (out2 / "recall_probability_surface.png").read_bytes()


def test_replay_missing_recall_data_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        plot_recall_surface_from_saved(str(tmp_path))


def test_replay_skips_missing_variant(recall_df, fitted_results, saved_recall_dir, tmp_path):
    os.remove(saved_recall_dir / "recall_variant_B_summary.tsv")
    out = tmp_path / "out"
    out.mkdir()
    plot_recall_surface_from_saved(str(saved_recall_dir), output_dir=str(out))
    assert (out / "recall_probability_surface.png").exists()
