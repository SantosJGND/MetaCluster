"""
Tests for feature-importance extraction/persistence in recall and composition modellers.
"""

import numpy as np
import pandas as pd
import pytest

from metagenomics_utils.overlap_manager.om_models import (
    BaseCompositionModeller,
    CutoffRecallModeller,
    DirectXGBRecallModeller,
    GBCompositionModeller,
    GPCLFRecallModeller,
    LRCompositionModeller,
    RecallModeller,
    RFCompositionModeller,
    XGBCompositionModeller,
)

STATS_COLS = ["n_leaves", "tax_diversity", "Min_Dist", "Min_Shared"]
TAXA_COLS = ["taxon_B", "taxon_A"]
COMPOSITION_COLS = STATS_COLS + TAXA_COLS


@pytest.fixture
def synthetic_composition_data():
    rng = np.random.RandomState(42)
    X = pd.DataFrame(rng.rand(80, len(COMPOSITION_COLS)), columns=COMPOSITION_COLS)
    y = pd.Series((X["Min_Dist"] + 0.2 * X["taxon_B"] > 0.5).astype(int))
    return X, y


class _DummyCompositionModeller(BaseCompositionModeller):
    """Composition modeller with a classifier lacking any importance attribute."""

    def _build_pipeline(self, X_train, y_train):
        from sklearn.dummy import DummyClassifier
        from sklearn.pipeline import Pipeline

        return Pipeline([("classifier", DummyClassifier(strategy="most_frequent"))]).fit(X_train, y_train)


class TestCompositionFeatureImportances:
    def test_tree_modellers_return_sorted_series(self, synthetic_composition_data):
        X, y = synthetic_composition_data
        for modeller in [
            XGBCompositionModeller(n_estimators=20, max_depth=3),
            RFCompositionModeller(n_estimators=20, max_depth=5),
            GBCompositionModeller(n_estimators=20, max_depth=3),
        ]:
            fitted = modeller.fit(X, y)
            importances = fitted.feature_importance_dataframe()

            assert isinstance(importances, pd.Series)
            assert len(importances) == len(COMPOSITION_COLS)
            assert set(importances.index) == set(COMPOSITION_COLS)
            assert importances.is_monotonic_decreasing
            assert importances.notna().all()

    def test_lr_modeller_uses_coefficients(self, synthetic_composition_data):
        X, y = synthetic_composition_data
        modeller = LRCompositionModeller(max_iter=2000).fit(X, y)

        importances = modeller.feature_importance_dataframe()

        assert isinstance(importances, pd.Series)
        assert set(importances.index) == set(STATS_COLS)
        assert (importances >= 0).all()

    def test_modeller_without_importance_returns_none(self, synthetic_composition_data):
        X, y = synthetic_composition_data
        modeller = _DummyCompositionModeller().fit(X, y)

        assert modeller.feature_importance_dataframe() is None

    def test_save_writes_tsv(self, synthetic_composition_data, tmp_path):
        X, y = synthetic_composition_data
        modeller = XGBCompositionModeller(n_estimators=20, max_depth=3).fit(X, y)

        modeller.save_feature_importances(str(tmp_path))

        out = tmp_path / "composition_feature_importances.tsv"
        assert out.exists()
        content = out.read_text()
        for feature in COMPOSITION_COLS:
            assert feature in content


class TestRecallFeatureImportances:
    def _stub_modeller(self, modeller, est_with_importances):
        class _FakeEst:
            feature_importances_ = np.array([0.9, 0.1])

        class _FakeMulti:
            estimators_ = [_FakeEst(), _FakeEst()]

        modeller.RecP_feature_cols = ["f1", "f2"]
        modeller.model = _FakeMulti() if est_with_importances else type("M", (), {})()
        return modeller

    def test_multiple_output_aggregated_series(self):
        modeller = self._stub_modeller(RecallModeller(data_set_divide=2), est_with_importances=True)

        importances = modeller.feature_importance_dataframe()

        assert isinstance(importances, pd.Series)
        assert list(importances.index) == ["f1", "f2"]
        assert importances["f1"] == pytest.approx(0.9)
        assert importances.is_monotonic_decreasing

    def test_direct_xgb_uses_feature_importances(self):
        modeller = DirectXGBRecallModeller(data_set_divide=2)
        X = pd.DataFrame(np.random.RandomState(0).rand(40, 2), columns=["f1", "f2"])
        y = pd.Series(X["f1"] * 0.5 + np.random.RandomState(1).rand(40) * 0.01)
        from sklearn.ensemble import RandomForestRegressor

        modeller.model = RandomForestRegressor(n_estimators=10, random_state=42).fit(X, y)
        modeller.RecP_feature_cols = ["f1", "f2"]

        importances = modeller.feature_importance_dataframe()

        assert isinstance(importances, pd.Series)
        assert set(importances.index) == {"f1", "f2"}
        assert importances.is_monotonic_decreasing

    def test_cutoff_uses_feature_importances(self):
        modeller = CutoffRecallModeller(data_set_divide=4, target_recall=1.0)
        X = pd.DataFrame(np.random.RandomState(0).rand(40, 3), columns=["f1", "f2", "f3"])
        y = pd.Series([i % 3 for i in range(40)])
        from sklearn.ensemble import RandomForestClassifier

        modeller.model = RandomForestClassifier(n_estimators=10, random_state=42).fit(X, y)
        modeller.RecP_feature_cols = ["f1", "f2", "f3"]

        importances = modeller.feature_importance_dataframe()

        assert isinstance(importances, pd.Series)
        assert set(importances.index) == {"f1", "f2", "f3"}
        assert importances.is_monotonic_decreasing

    def test_unsupported_returns_none(self):
        modeller = self._stub_modeller(RecallModeller(data_set_divide=2), est_with_importances=False)
        assert modeller.feature_importance_dataframe() is None

        gp_modeller = GPCLFRecallModeller(data_set_divide=2)
        gp_modeller.model = object()
        gp_modeller.RecP_feature_cols = ["f1", "f2"]
        assert gp_modeller.feature_importance_dataframe() is None

    def test_save_writes_tsv(self, tmp_path):
        modeller = self._stub_modeller(RecallModeller(data_set_divide=2), est_with_importances=True)

        modeller.save_feature_importances(str(tmp_path))

        out = tmp_path / "recall_feature_importances.tsv"
        assert out.exists()
        content = out.read_text()
        assert "f1" in content
        assert "f2" in content

    def test_save_skips_when_unsupported(self, tmp_path):
        modeller = self._stub_modeller(RecallModeller(data_set_divide=2), est_with_importances=False)

        modeller.save_feature_importances(str(tmp_path))

        assert not (tmp_path / "recall_feature_importances.tsv").exists()
