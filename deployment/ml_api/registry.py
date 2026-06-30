import sys
import time
import socket
import threading
from pathlib import Path
import joblib

_THIS_DIR = Path(__file__).parent
sys.path.append(str(_THIS_DIR))
sys.path.append(str(_THIS_DIR.parent.parent))

from config import get_logger, MLFLOW_TRACKING_URI, registry_name, MLFLOW_EXPERIMENT, ModelFile, RECALL_MODEL_VARIANTS, PROJECT_MODELS
from metagenomics_utils.overlap_manager.om_models import (
    ClusteringPipeline,
    GPCLFThreshold,
    RecallModeller,
    CutoffRecallModeller,
    DirectXGBRecallModeller,
    GPCLFRecallModeller,
    XGBCompositionModeller,
    OptunaXGBCompositionModeller,
    RFCompositionModeller,
    GBCompositionModeller,
    LRCompositionModeller,
)

logger = get_logger(__name__)

MODELS_DIR = Path(__file__).parent / "models"

_model_cache: dict[str, dict] = {}
_cache_lock = threading.Lock()


def all_model_keys() -> list:
    """Return all concrete cache keys (televir + recall variants)."""
    keys = list(PROJECT_MODELS)
    for variant in RECALL_MODEL_VARIANTS:
        keys.append(f"recall_{variant}")
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
    filename = ModelFile.project_files.get(model_type)
    if not filename:
        raise RuntimeError(f"No local filename configured for model type '{model_type}'")
    model_path = MODELS_DIR / filename
    if not model_path.exists():
        raise RuntimeError(f"Local model file not found: {model_path}")
    model = joblib.load(str(model_path))
    logger.info(f"Loaded local model '{model_type}' from {model_path}")
    return model


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
                versions = client.search_model_versions(
                    f"name='{name}'"
                )
                if versions:
                    v = versions[0]
                    return {"version": v.version, "stage": v.current_stage, "run_id": v.run_id}
        except Exception:
            pass
    return {"version": "local", "stage": "local", "run_id": None}


def load_and_cache(model_type: str) -> dict:
    model = load_production_model(model_type)
    version_info = get_model_version(model_type)
    entry = {
        "model": model,
        "version": version_info.get("version"),
        "stage": version_info.get("stage"),
        "run_id": version_info.get("run_id"),
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
        raise RuntimeError(
            f"Model '{model_type}' not in cache. Use POST /reload/{model_type} to load it."
        )
    return entry["model"], {
        "version": entry["version"],
        "stage": entry["stage"],
        "run_id": entry["run_id"],
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
            "version": e["version"],
            "stage": e["stage"],
            "loaded_seconds_ago": round(time.time() - e["loaded_at"], 1),
        }
        for mt, e in sorted(items)
    ]
