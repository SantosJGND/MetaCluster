from typing import Any

import requests


class MLAPIError(Exception):
    """Error raised when the ML API returns a non-2xx response."""

    def __init__(self, status_code: int, detail: str, response_body: Any = None):
        self.status_code = status_code
        self.detail = detail
        self.response_body = response_body
        super().__init__(f"HTTP {status_code}: {detail}")


class MLAPIClient:
    """Python client for the INSaFLU ML API.

    Models are referenced by composite keys ``{tax_level}_{category}_{variant}``,
    e.g. ``"order_recall_gp_clf"`` or ``"genus_composition_xgb"``.

    The API auto-discovers models by scanning the ``models/`` directory for pickle
    bundles containing ``model_category``, ``tax_level``, and ``model_type`` fields.
    """

    def __init__(self, base_url: str = "http://localhost:8000", timeout: int = 30):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _raise_for_response(self, r: requests.Response):
        try:
            err_body = r.json()
            detail = err_body.get("detail", r.text)
        except Exception:
            detail = r.text
        raise MLAPIError(r.status_code, detail, err_body if isinstance(detail, (dict, list)) else None)

    def _get(self, path: str, params: dict | None = None) -> dict[str, Any]:
        r = requests.get(f"{self.base_url}{path}", params=params, timeout=self.timeout)
        if not r.ok:
            self._raise_for_response(r)
        return r.json()

    def _post(self, path: str, body: dict | None = None) -> dict[str, Any]:
        r = requests.post(f"{self.base_url}{path}", json=body, timeout=self.timeout)
        if not r.ok:
            self._raise_for_response(r)
        return r.json()

    def health(self) -> dict[str, Any]:
        """GET /health — liveness check."""
        return self._get("/health")

    def models(self) -> dict[str, Any]:
        """GET /models — list cached models with version and stage info."""
        return self._get("/models")

    def reload(self, model_type: str | None = None) -> dict[str, Any]:
        """POST /reload or /reload/{model_type} — reload model cache.

        Args:
            model_type: Composite key (e.g. ``"order_recall_gp_clf"``).
                        Omit to reload all models.
        """
        if model_type:
            return self._post(f"/reload/{model_type}")
        return self._post("/reload")

    def composition_model_tax_level(self, model: str | None = None) -> dict[str, Any]:
        """GET /composition_model_tax_level — return the tax_level of a composition model.

        Args:
            model: Full composite key (e.g. ``"order_composition_rf"``) or short
                   variant (e.g. ``"rf"``). Omit to list all composition models.
        """
        return self._get("/composition_model_tax_level", params={"model": model} if model else None)

    def predict_recall_cutoff(
        self,
        rows: list[dict[str, Any]],
        model: str = "order_recall_gp_clf",
        target_recall: float | None = None,
        confidence: float | None = None,
    ) -> dict[str, Any]:
        """POST /predict_recall_cutoff_from_table — predict recall cutoff from raw table rows.

        Args:
            rows: List of dicts with keys ``taxid``, ``total_uniq_reads``, ``order``, ``family``.
            model: Full composite key from ``GET /models`` (e.g. ``"order_recall_gp_clf"``).
            target_recall: Optional target recall threshold.
            confidence: Optional confidence level for probability-guided cutoff.
        """
        body: dict[str, Any] = {
            "model": model,
            "rows": rows,
        }
        if target_recall is not None:
            body["target_recall"] = target_recall
        if confidence is not None:
            body["confidence"] = confidence
        return self._post("/predict_recall_cutoff_from_table", body)

    def predict_clustering_threshold(self, features: dict[str, Any]) -> dict[str, Any]:
        """POST /predict_televir_clustering_threshold — **deprecated**, returns 501."""
        return self._post("/predict_televir_clustering_threshold", features)

    def predict_composition_stop_traversal(self, features: dict[str, float], model: str = "xgb") -> dict[str, Any]:
        """POST /predict_composition_stop_traversal — predict stop_traversal from node features.

        Args:
            features: Dict of node feature values matching training column names.
            model: Composition variant (e.g. ``"xgb"``, ``"rf"``, ``"lr"``).
        """
        return self._post("/predict_composition_stop_traversal", {"model": model, "features": features})
