import socket
import sys
import threading
import time
from pathlib import Path

import joblib

_THIS_DIR = Path(__file__).parent
sys.path.append(str(_THIS_DIR))
sys.path.append(str(_THIS_DIR.parent.parent))

from config import (
    MLFLOW_EXPERIMENT,
    MLFLOW_TRACKING_URI,
    MODEL_TYPE_MAP,
    get_logger,
    registry_name,
)

logger = get_logger(__name__)

MODELS_DIR = Path(__file__).parent / "models"

_model_cache: dict[str, dict] = {}
_cache_lock = threading.Lock()

_ALL_MODELS: dict[str, dict] = {}

def discover_models():
    """Scan models/ directory, read bundle metadata, build composite key map.

    Each key is ``{tax_level}_{category}_{variant}`` (e.g. ``order_recall_xgb``).
    """
    discovered = {}
    if not MODELS_DIR.is_dir():
        logger.warning("Models directory does not exist: %s", MODELS_DIR)
        return discovered

    GENUS_ALIASES = {"genus", "g", "genera"}

    for p in sorted(MODELS_DIR.iterdir()):
        if p.suffix not in (".pkl",):
            continue
        try:
            bundle = joblib.load(str(p))
        except Exception as exc:
            logger.warning("Could not load %s: %s", p.name, exc)
            continue
        if not isinstance(bundle, dict):
            logger.warning("Skipping %s: not a dict bundle", p.name)
            continue

        category = bundle.get("model_category", "unknown")
        tax_level = bundle.get("tax_level", "unknown")
        model_type = bundle.get("model_type", p.stem)
        description = bundle.get("description", "")
        date_trained = bundle.get("date_trained", "")

        if tax_level in GENUS_ALIASES:
            tax_level = "genus"

        key = f"{tax_level}_{category}_{model_type}"
        discovered[key] = {
            "path": str(p),
            "model_category": category,
            "tax_level": tax_level,
            "model_type": model_type,
            "description": description,
            "date_trained": date_trained,
            "bundle": bundle,
        }

    _ALL_MODELS.clear()
    _ALL_MODELS.update(discovered)
    logger.info("Discovered %d models", len(discovered))
    for k, v in discovered.items():
        logger.debug("  %s — %s", k, v["description"])
    return discovered


def all_model_keys() -> list:
    keys = list(_ALL_MODELS.keys())
    if not keys:
        discover_models()
        keys = list(_ALL_MODELS.keys())
    return keys


def _mlflow_reachable(timeout=3):
    host = MLFLOW_TRACKING_URI.replace("http://", "").replace("https://", "").split(":")[0]
    port = int(MLFLOW_TRACKING_URI.split(":")[-1])
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(timeout)
        result = s.connect_ex((host, port))
        s.close()
        return result == 0
    except Exception:
        return False


def _load_local_model(model_type):
    if model_type in _ALL_MODELS:
        return _ALL_MODELS[model_type]["bundle"]
    raise RuntimeError(f"No local model found for '{model_type}' — run discover_models() first")


def load_production_model(model_type):
    if _mlflow_reachable():
        try:
            import mlflow

            mlflow.set_tracking_uri(MLFLOW_TRACKING_URI)
            mlflow.set_experiment(MLFLOW_EXPERIMENT)
            name = registry_name(model_type)
            stages = ["Production", "Staging"]

            for stage in stages:
                uri = f"models:/{name}/{stage}"
                try:
                    model = mlflow.sklearn.load_model(uri)
                    logger.info(f"Loaded {name} from {stage}")
                    return model
                except Exception as e:
                    logger.warning(f"Failed to load {name} from {stage}: {e}")
        except Exception as e:
            logger.warning(f"MLflow error: {e}")
    else:
        logger.info(f"MLflow not reachable at {MLFLOW_TRACKING_URI}, using local file")

    return _load_local_model(model_type)


def get_model_version(model_type):
    if _mlflow_reachable():
        try:
            import mlflow
            from mlflow.tracking import MlflowClient

            mlflow.set_tracking_uri(MLFLOW_TRACKING_URI)
            mlflow.set_experiment(MLFLOW_EXPERIMENT)
            name = registry_name(model_type)
            client = MlflowClient()
            for stage in ["Production", "Staging"]:
                versions = client.search_model_versions(f"name='{name}'")
                if versions:
                    v = versions[0]
                    return {"version": v.version, "stage": v.current_stage, "run_id": v.run_id}
        except Exception:
            pass
    return {"version": "local", "stage": "local", "run_id": None}


def load_and_cache(model_type: str) -> dict:
    model = load_production_model(model_type)
    version_info = get_model_version(model_type)
    meta = _ALL_MODELS.get(model_type, {})
    entry = {
        "model": model,
        "type": MODEL_TYPE_MAP.get(model_type, "unknown"),
        "version": version_info.get("version"),
        "stage": version_info.get("stage"),
        "run_id": version_info.get("run_id"),
        "description": meta.get("description", ""),
        "date_trained": meta.get("date_trained", ""),
        "tax_level": meta.get("tax_level", "unknown"),
        "loaded_at": time.time(),
    }
    with _cache_lock:
        _model_cache[model_type] = entry
    registered = registry_name(model_type)
    logger.info(f"Cached {registered} v{entry['version']} ({entry['stage']})")
    return entry


def get_cached_model(model_type: str) -> tuple:
    with _cache_lock:
        entry = _model_cache.get(model_type)
    if entry is None:
        raise RuntimeError(f"Model '{model_type}' not in cache. Use POST /reload/{model_type} to load it.")
    return entry["model"], {
        "version": entry["version"],
        "stage": entry["stage"],
        "run_id": entry["run_id"],
        "tax_level": entry["tax_level"],
    }


def invalidate_cache(model_type: str = None):
    with _cache_lock:
        if model_type is None:
            _model_cache.clear()
            logger.info("Cleared entire model cache")
        else:
            _model_cache.pop(model_type, None)
            logger.info(f"Cleared cache for {model_type}")


def cache_status() -> list[dict]:
    with _cache_lock:
        items = list(_model_cache.items())
    return [
        {
            "model_type": mt,
            "type": e["type"],
            "version": e["version"],
            "stage": e["stage"],
            "description": e.get("description", ""),
            "date_trained": e.get("date_trained", ""),
            "tax_level": e.get("tax_level", "unknown"),
            "loaded_seconds_ago": round(time.time() - e["loaded_at"], 1),
        }
        for mt, e in sorted(items)
    ]


def composition_models_tax_level(model: str | None = None) -> list[dict]:
    """Return tax_level metadata for composition models.

    Args:
        model: Optional full composite key (e.g. ``order_composition_rf``) or
               short model_type variant (e.g. ``rf``). Omit for all composition models.

    Returns:
        List of dicts with ``key``, ``model_type``, ``tax_level``, ``category``,
        ``description``, and ``date_trained``.
    """
    if not _ALL_MODELS:
        discover_models()
    comp = [
        (k, v)
        for k, v in sorted(_ALL_MODELS.items())
        if v.get("model_category") == "composition"
    ]
    if model:
        comp = [(k, v) for k, v in comp if k == model or v.get("model_type") == model]
    return [
        {
            "key": k,
            "model_type": v.get("model_type"),
            "tax_level": v.get("tax_level", "unknown"),
            "category": v.get("model_category"),
            "description": v.get("description", ""),
            "date_trained": v.get("date_trained", ""),
        }
        for k, v in comp
    ]
