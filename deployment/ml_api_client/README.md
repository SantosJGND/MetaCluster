# INSaFLU ML API Client

Python client and CLI tool for the INSaFLU ML prediction API.

## Installation

```bash
pip install -r requirements.txt
```

## Python usage

```python
from ml_api_client import MLAPIClient

c = MLAPIClient(base_url="http://localhost:8000")

# Liveness
c.health()           # -> {"status": "ok"}

# Model cache
c.models()           # -> {"models": [...]}
c.reload()           # reload all
c.reload("order_recall_gp_clf")  # reload one model (composite key)

# Recall cutoff prediction from raw table rows
rows = [
    {"taxid": 2697049, "total_uniq_reads": 2348731, "best_match_is_best": False},
]
c.predict_recall_cutoff(rows, model="gp_clf", tax_level="order", target_recall=0.95)

# Composition stop-traversal prediction from node features
c.predict_composition_stop_traversal({
    "n_leaves": 12,
    "tax_diversity": 1.2,
    "Min_Dist": 0.15,
    "Min_Shared": 0.3,
    "Coronaviridae": 0.021,
}, model="xgb")
```

> **Note:** The `/predict_televir_clustering_threshold` endpoint is **deprecated** and returns HTTP 501.
> The corresponding method `c.predict_clustering_threshold()` is kept for backward compatibility
> but will always raise an error.

## CLI

```bash
python -m ml_api_client.cli --help
```

Run with `PYTHONPATH=deployment` from the repo root, or install the package.

### Commands

| Command | Description |
|---|---|
| `mlapi health` | Liveness check |
| `mlapi models` | List cached models with version/stage |
| `mlapi reload [model_type]` | Reload model cache (all or one — uses composite keys like `order_recall_gp_clf`) |
| `mlapi predict-recall --rows file.csv [options]` | Predict recall cutoff from raw table |
| `mlapi predict-composition-stop --features file.json` | Predict stop_traversal from node features |

**Deprecated:**

| Command | Description |
|---|---|
| `mlapi predict-clustering --features file.json` | **Deprecated** — endpoint returns 501 |

### predict-recall options

| Flag | Default | Description |
|---|---|---|
| `--rows` | (required) | CSV with `taxid,total_uniq_reads[,best_match_is_best]` |
| `--model` | `gp_clf` | Recall model variant (`gp_clf`, `xgb_direct`, `xgb_multi`) |
| `--tax-level` | `order` | Taxonomic level the model was trained on |
| `--target-recall` | — | Target recall threshold (0–1) |
| `--confidence` | — | Confidence threshold for probability-guided cutoff (0–1) |

### CSV format (predict-recall)

```csv
taxid,total_uniq_reads,best_match_is_best
2697049,2348731,false
1122334,50000,true
```

### JSON format (predict-composition-stop)

Flat dict of node features matching training column names:

```json
{
  "n_leaves": 12,
  "tax_diversity": 1.2,
  "Min_Dist": 0.15,
  "Min_Shared": 0.3,
  "Coronaviridae": 0.021,
  "Peduoviridae": 0.0,
  "Straboviridae": 0.0
}
```

### Example workflow

```bash
# 1. Check the API is up
mlapi health

# 2. Load models into cache
mlapi reload

# 3. Predict recall cutoff from raw table rows
mlapi predict-recall --rows hits.csv --model gp_clf --tax-level order

# 4. Predict stop_traversal from node features
mlapi predict-composition-stop --features node_features.json
```

## Error handling

All errors are printed as JSON to stderr with exit code 1:

```bash
$ mlapi health
{"error": "HTTPConnectionPool(host='localhost', port=8000): ..."}
```
