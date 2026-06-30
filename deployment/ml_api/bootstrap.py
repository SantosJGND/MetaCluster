import sys
import joblib
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))
import mlflow

from config import get_logger, MLFLOW_TRACKING_URI, MLFLOW_EXPERIMENT, registry_name, ModelFile
from registry import all_model_keys

logger = get_logger(__name__)

MODELS_DIR = Path(__file__).parent / "models"

def _get_client() -> "MlflowClient":
    from mlflow.tracking import MlflowClient
    mlflow.set_tracking_uri(MLFLOW_TRACKING_URI)
    return MlflowClient()


def transition_to_production(registered_name, version):
    client = _get_client()
    client.transition_model_version_stage(
        name=registered_name,
        version=version,
        stage="Production"
    )
    logger.info(f"{registered_name} v{version} → Production")
    for v in client.search_model_versions(f"name='{registered_name}'"):
        if v.version != version and v.current_stage == "Production":
            client.transition_model_version_stage(
                name=registered_name,
                version=v.version,
                stage="Archived"
            )
            logger.info(f"{registered_name} v{v.version} → Archived")


def _unwrap_model(bundle):
    """Extract the inner sklearn-compatible model from a dict bundle or pass through."""
    if isinstance(bundle, dict):
        return bundle.get('pipeline') or bundle.get('model') or bundle
    return bundle


def bootstrap_model(model_type):
    import mlflow
    import mlflow.sklearn

    bundle_name = ModelFile.project_files.get(model_type, f"{model_type}.joblib")
    model_path = MODELS_DIR / bundle_name
    if not model_path.exists():
        logger.warning(f"No model file found for {model_type} at {model_path}")
        return False

    bundle = joblib.load(str(model_path))
    model = _unwrap_model(bundle)
    registered_name = registry_name(model_type)

    mlflow.set_tracking_uri(MLFLOW_TRACKING_URI)
    mlflow.set_experiment(MLFLOW_EXPERIMENT)

    with mlflow.start_run(run_name=f"bootstrap_{model_type}") as run:
        run_id = run.info.run_id
        mlflow.log_param("model_type", model_type)
        mlflow.log_param("source", "bootstrap")
        mlflow.sklearn.log_model(sk_model=model, artifact_path="model")

        result = mlflow.register_model(f"runs:/{run_id}/model", registered_name)
        logger.info(f"{model_type} registered as {registered_name} v{result.version}")

        transition_to_production(registered_name, result.version)
        logger.info(f"{model_type} promoted to Production as v{result.version}")

    return True


def main():
    logger.info("Bootstrapping existing models into MLflow Registry...")
    success = 0
    for mt in all_model_keys():
        if bootstrap_model(mt):
            success += 1
            
    logger.info(f"Bootstrapped {success}/{len(all_model_keys())} models")
    if success == 0:
        logger.warning("No models were bootstrapped. Train models first!")
        sys.exit(1)


if __name__ == "__main__":
    main()
