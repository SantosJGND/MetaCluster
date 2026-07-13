import logging
import os

PROJECT_NAME = "INSaFLU_ML_API"
RANDOM_SEED = 42
DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

MLFLOW_TRACKING_URI = os.getenv("MLFLOW_TRACKING_URI", "http://mlflow:5000")
MLFLOW_EXPERIMENT = os.getenv("MLFLOW_EXPERIMENT", PROJECT_NAME)
MLRUNS_DIR = os.getenv("MLRUNS_DIR", "./mlruns")
FEATURES_DIR = os.path.join(DATA_DIR, "features")

RECALL_MODEL_VARIANTS = {
    "xgb_direct": {"file": "direct_xgb_bundle.pkl", "cls": "DirectXGBRecallModeller"},
    "xgb_multi": {"file": "recall_xgb_bundle.pkl", "cls": "RecallModeller"},
    "gp_clf": {"file": "recall_gp_clf_pipeline.pkl", "cls": "GPCLFRecallModeller"},
}

COMPOSITION_MODEL_VARIANTS = {
    "xgb": {"file": "composition_xgb_bundle.pkl", "cls": "XGBCompositionModeller"},
    "xgb_optimized": {"file": "composition_optuna_bundle.pkl", "cls": "OptunaXGBCompositionModeller"},
    "rf": {"file": "composition_rf_bundle.pkl", "cls": "RFCompositionModeller"},
    "gb": {"file": "composition_gb_bundle.pkl", "cls": "GBCompositionModeller"},
    "lr": {"file": "composition_lr_bundle.pkl", "cls": "LRCompositionModeller"},
}

MODEL_TYPE_MAP = {
    **{f"recall_{k}": "recall" for k in RECALL_MODEL_VARIANTS},
    **{f"composition_{k}": "composition" for k in COMPOSITION_MODEL_VARIANTS},
}


def registry_name(model_type):
    return f"{PROJECT_NAME}_{model_type}"


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")


def get_logger(name):
    return logging.getLogger(name)
