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
c.reload("recall_gp_clf")  # reload one model variant

# Recall cutoff prediction from raw table rows
rows = [
    {"taxid": 2697049, "total_uniq_reads": 2348731, "best_match_is_best": False},
]
c.predict_recall_cutoff(rows, model="gp_clf", target_recall=0.95)

# Clustering threshold prediction
c.predict_clustering_threshold({
    "classifiers": ["xgb"],
    "taxa": ["Coronaviridae"],
    "taxon_hits": [0, 0, 0.021, 0, 0, 0, 0, 0, 0, 0, 0, 0.979],
    "taxonomic_diversity": 0.7,
    "n_leaves": 5,
    "proportion_shared_hits": 0.3,
    "proportion_unique_hits": 0.7,
})

# Composition stop-traversal prediction from node features
c.predict_composition_stop_traversal({
    "n_leaves": 12,
    "tax_diversity": 1.2,
    "Min_Dist": 0.15,
    "Min_Shared": 0.3,
    "Coronaviridae": 0.021,
})
```

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
| `mlapi reload [model_type]` | Reload model cache (all or one) |
| `mlapi predict-recall --rows file.csv` | Predict recall cutoff from raw table |
| `mlapi predict-clustering --features file.json` | Predict TELEVIR clustering threshold |
| `mlapi predict-composition-stop --features file.json` | Predict stop_traversal from node features |

### CSV format (predict-recall)

```csv
taxid,total_uniq_reads,best_match_is_best
2697049,2348731,false
1122334,50000,true
```

### JSON format (predict-clustering)

```json
{
  "classifiers": ["xgb"],
  "taxa": ["Coronaviridae"],
  "taxon_hits": [0, 0, 0.021, 0, 0, 0, 0, 0, 0, 0, 0, 0.979],
  "taxonomic_diversity": 0.7,
  "n_leaves": 5,
  "proportion_shared_hits": 0.3,
  "proportion_unique_hits": 0.7
}
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
mlapi predict-recall --rows hits.csv --model gp_clf

# 4. Predict clustering threshold
mlapi predict-clustering --features sample_cluster.json

# 5. Predict stop_traversal from node features
mlapi predict-composition-stop --features node_features.json
```

## Error handling

All errors are printed as JSON to stderr with exit code 1:

```bash
$ mlapi health
{"error": "HTTPConnectionPool(host='localhost', port=8000): ..."}
```
