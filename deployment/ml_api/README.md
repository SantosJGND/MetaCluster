# INSaFLU ML API

FastAPI service serving ML models for the INSaFLU platform.

## Endpoints

| Method | Path | Description |
|---|---|---|
| `GET`  | `/health`                              | Liveness check |
| `GET`  | `/models`                              | List cached models with version/stage |
| `GET`  | `/composition_model_tax_level`         | Taxonomic level of composition model(s). Optional `?model=` (composite key or short variant). |
| `POST` | `/reload`                              | Reload all models from MLflow registry (or local files if MLflow unreachable) |
| `POST` | `/reload/{model_type}`                 | Reload a single model (composite key, see below) |
| `POST` | `/predict_recall_cutoff_from_table`    | Recall cutoff prediction from raw table rows. Accepts `tax_level` in body. |
| `POST` | `/predict_composition_stop_traversal`  | Stop-traversal prediction from node features. Accepts `tax_level` in body. |

## Model discovery and cache keys

At startup, `discover_models()` scans the `models/` directory for `.pkl` files,
reads each bundle dict, and registers them under a composite key:

```
{tax_level}_{category}_{model_type}
```

Examples: `order_recall_gp_clf`, `genus_composition_xgb`, `family_recall_xgb_multi`.

### Required bundle metadata fields

Each pickle bundle must be a **dict** containing at least:

| Field | Type | Description |
|---|---|---|
| `model_category` | `str` | `"recall"`, `"composition"`, or `"crosshit"` |
| `tax_level` | `str` | Taxonomic level the model was trained on (`"order"`, `"family"`, `"genus"`, etc.) |
| `model_type` | `str` | Variant identifier (e.g. `"gp_clf"`, `"xgb"`, `"rf"`) |
| `description` | `str` | Human-readable description of the model |
| `date_trained` | `str` | ISO 8601 timestamp of training |
| `pipeline` / `model` | estimator | The fitted sklearn-compatible estimator |

### Currently registered model variants

| Composite key pattern | Category | Model types |
|---|---|---|
| `{tax_level}_recall_gp_clf` | recall | GP+CLF recall |
| `{tax_level}_recall_xgb_direct` | recall | Direct XGBoost recall |
| `{tax_level}_recall_xgb_multi` | recall | Multi-output XGBoost recall |
| `{tax_level}_composition_xgb` | composition | XGBoost stop-traversal |
| `{tax_level}_composition_xgb_optimized` | composition | XGBoost + Optuna stop-traversal |
| `{tax_level}_composition_rf` | composition | Random Forest stop-traversal |
| `{tax_level}_composition_gb` | composition | Gradient Boosting stop-traversal |
| `{tax_level}_composition_lr` | composition | Logistic Regression stop-traversal (stats-only) |

## Model serving

Models are loaded from local pickle files in `models/` (or from MLflow if reachable) into an in-memory cache at startup.
To reload after replacing a pickle file, use the composite key:

    curl -X POST http://localhost:8000/reload/order_recall_gp_clf
    curl -X POST http://localhost:8000/reload/order_composition_xgb

To look up the taxonomic level of a composition model (composite key or short variant):

    curl -s http://localhost:8000/composition_model_tax_level?model=rf
    # {"key": "order_composition_rf", "model_type": "rf", "tax_level": "order", ...}

    curl -s http://localhost:8000/composition_model_tax_level
    # {"models": [{"key": "order_composition_rf", "model_type": "rf", "tax_level": "order", ...}, ...]}

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
python train_recall.py --model gp_clf --tax-level order

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
# Recall cutoff from raw table rows (tax_level defaults to "order")
curl -s -X POST 'http://localhost:8000/predict_recall_cutoff_from_table' \
  -H 'Content-Type: application/json' \
  -d '{
    "model": "gp_clf",
    "tax_level": "order",
    "rows": [{"taxid": 2697049, "total_uniq_reads": 2348731, "best_match_is_best": false}],
    "target_recall": 0.95
  }'

# Composition stop-traversal from node features
curl -s -X POST 'http://localhost:8000/predict_composition_stop_traversal' \
  -H 'Content-Type: application/json' \
  -d '{
    "model": "xgb",
    "tax_level": "order",
    "features": {
      "n_leaves": 12,
      "tax_diversity": 1.2,
      "Min_Dist": 0.15,
      "Min_Shared": 0.3
    }
  }'
```

For more examples see `deployment/ml_api_client/README.md`.
