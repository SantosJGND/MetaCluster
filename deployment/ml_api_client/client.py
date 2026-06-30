import requests
from typing import Any, Dict, List, Optional


class MLAPIClient:
    def __init__(self, base_url: str = "http://localhost:8000", timeout: int = 30):
        self.base_url = base_url.rstrip("/")
        self.timeout = timeout

    def _get(self, path: str, params: dict | None = None) -> Dict[str, Any]:
        r = requests.get(f"{self.base_url}{path}", params=params, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def _post(self, path: str, body: dict | None = None) -> Dict[str, Any]:
        r = requests.post(f"{self.base_url}{path}", json=body, timeout=self.timeout)
        r.raise_for_status()
        return r.json()

    def health(self) -> Dict[str, Any]:
        return self._get("/health")

    def models(self) -> Dict[str, Any]:
        return self._get("/models")

    def reload(self, model_type: Optional[str] = None) -> Dict[str, Any]:
        if model_type:
            return self._post(f"/reload/{model_type}")
        return self._post("/reload")

    def predict_recall_cutoff(
        self,
        rows: List[Dict[str, Any]],
        model: str = "gp_clf",
        tax_level: str = "order",
        target_recall: Optional[float] = None,
        confidence: Optional[float] = None,
    ) -> Dict[str, Any]:
        body: Dict[str, Any] = {
            "model": model,
            "rows": rows,
            "tax_level": tax_level,
        }
        if target_recall is not None:
            body["target_recall"] = target_recall
        if confidence is not None:
            body["confidence"] = confidence
        return self._post("/predict_recall_cutoff_from_table", body)

    def predict_clustering_threshold(self, features: Dict[str, Any]) -> Dict[str, Any]:
        return self._post("/predict_televir_clustering_threshold", features)

    def predict_composition_stop_traversal(self, features: Dict[str, float]) -> Dict[str, Any]:
        return self._post("/predict_composition_stop_traversal", {"features": features})
