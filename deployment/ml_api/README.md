# INSaFLU ML API

FastAPI service serving ML models registered in MLflow for the INSaFLU platform.

## Endpoints

| Method | Path | Description |
|---|---|---|
| `GET`  | `/health`                              | Liveness check |
| `GET`  | `/models`                              | List cached models with version/stage |
| `POST` | `/reload`                              | Reload all models from MLflow registry (or local files if MLflow unreachable) |
| `POST` | `/reload/{model_type}`                 | Reload a single model |
| `POST` | `/predict_recall_cutoff_from_table`    | Recall cutoff prediction from raw table rows |
| `POST` | `/predict_composition_stop_traversal`  | Stop-traversal prediction from node features |
| `POST` | `/predict_televir_clustering_threshold` | TELEVIR clustering threshold prediction |

## Model types

| Cache key | Model | Class |
|---|---|---|
| `recall_gp_clf` | GP+CLF recall | `GPCLFRecallModeller` |
| `recall_xgb_direct` | Direct XGBoost recall | `DirectXGBRecallModeller` |
| `recall_xgb_multi` | Multi-output XGBoost recall | `RecallModeller` |
| `composition` | Stop-traversal classifier | `XGBCompositionModeller` (or `rf`, `gb`, `lr`, `xgb_optimized`) |
| `televir_clustering` | Clustering threshold | `ClusteringPipeline` |

## Model serving

Models are loaded from local pickle files in `models/` (or from MLflow if reachable) into an in-memory cache at startup.
To reload after replacing a pickle file:

    curl -X POST http://localhost:8000/reload/recall_gp_clf
    curl -X POST http://localhost:8000/reload/composition

## Environment

- `MLFLOW_TRACKING_URI` — MLflow server URL (default: `http://mlflow:5000`)
- `MLFLOW_EXPERIMENT` — MLflow experiment name (default: `INSaFLU_ML_API`)

## Build & run

    docker build -t ml_api .
    docker run -p 8000:8000 -e MLFLOW_TRACKING_URI=http://localhost:5000 ml_api

## Workflow

### 1. Train a model

```bash
# Recall — train a GP+CLF model variant
python train_recall.py --model gp_clf

# Composition — train via the model_evaluation pipeline
# (trains during evaluate.py, saves composition_xgb_bundle.pkl)
```

### 2. Copy the pickled bundle to ml_api/models/

```bash
# Recall variants are saved to models/ by train_recall.py
# Composition same: copy from training output
cp /path/to/trained/composition_xgb_bundle.pkl deployment/ml_api/models/
```

### 3. Test the endpoint

```bash
# Recall cutoff from raw table rows
curl -s -X POST 'http://localhost:8000/predict_recall_cutoff_from_table' \
  -H 'Content-Type: application/json' \
  -d '{
    "model": "gp_clf",
    "rows": [{"taxid": 2697049, "total_uniq_reads": 2348731, "best_match_is_best": false}],
    "target_recall": 0.95
  }'

# Composition stop-traversal from node features
curl -s -X POST 'http://localhost:8000/predict_composition_stop_traversal' \
  -H 'Content-Type: application/json' \
  -d '{
    "features": {
      "n_leaves": 12,
      "tax_diversity": 1.2,
      "Min_Dist": 0.15,
      "Min_Shared": 0.3
    }
  }'
```

For more examples, see `test_recall_payloads.json` and `deployment/ml_api_client/README.md`.
