from fastapi import FastAPI, HTTPException, Response
import numpy as np
import pandas as pd
from validation.schemas import (
    RecallCutoffFromTableRequest,
    RecallCutoffFromTableRow,
    TelevirClusteringThresholdRequest,
    ClusteringThresholdResult,
    RecallCutoffResult,
    CompositionStopTraversalRequest,
    CompositionStopTraversalResult,
)
import time
from monitoring import log_prediction, get_metrics
from registry import load_and_cache, get_cached_model, invalidate_cache, cache_status, all_model_keys
from config import registry_name, get_logger, RECALL_MODEL_VARIANTS

logger = get_logger(__name__)

app = FastAPI(title="INSaFLU ML API", description="API for INSaFLU machine learning models", version="1.0.0")


@app.on_event("startup")
def warmup_cache():
    for mt in all_model_keys():
        try:
            load_and_cache(mt)
        except Exception as e:
            logger.warning(f"Failed to preload '{mt}': {e}")
    logger.info(f"Cache warmup complete ({cache_status()})")


@app.get("/health")
def health():
    return {"status": "ok"}


@app.get("/models")
def list_models():
    return {"models": cache_status()}


@app.post("/reload")
def reload_all():
    invalidate_cache()
    count = 0
    errors = []
    for mt in all_model_keys():
        try:
            load_and_cache(mt)
            count += 1
        except Exception as e:
            errors.append({"model_type": mt, "error": str(e)})
    return {"reloaded": count, "errors": errors}


@app.post("/reload/{model_type}")
def reload_one(model_type: str):
    known = all_model_keys()
    if model_type not in known:
        raise HTTPException(
            status_code=404,
            detail=f"Unknown model type: '{model_type}'. Available: {known}",
        )
    invalidate_cache(model_type)
    try:
        load_and_cache(model_type)
    except RuntimeError as e:
        raise HTTPException(status_code=503, detail=str(e))
    return {"reloaded": model_type}


def _interpolate_recall_at_cutoff(cutoff, raw_recalls, n_divisions):
    """Interpolate recall at a given cutoff fraction using per-division predictions."""
    fractions = np.arange(1, n_divisions + 1) / n_divisions
    for i in range(n_divisions):
        if cutoff <= fractions[i]:
            if i == 0:
                return float(raw_recalls[i] * (cutoff / fractions[i]))
            f_low, f_high = fractions[i - 1], fractions[i]
            r_low, r_high = raw_recalls[i - 1], raw_recalls[i]
            if f_high == f_low:
                return float(r_high)
            return float(r_low + (cutoff - f_low) * (r_high - r_low) / (f_high - f_low))
    return float(raw_recalls[-1])


_RECALL_VARIANT_KEYS = {v: f"recall_{v}" for v in RECALL_MODEL_VARIANTS}


@app.post("/predict_recall_cutoff_from_table", response_model=RecallCutoffResult)
def predict_recall_cutoff_from_table(request: RecallCutoffFromTableRequest, response: Response):
    start = time.time()
    variant_key = _RECALL_VARIANT_KEYS.get(request.model)
    if variant_key is None:
        raise HTTPException(422, detail=f"Unknown recall model '{request.model}'. Choose from: {list(_RECALL_VARIANT_KEYS)}")
    try:
        bundle, version_info = get_cached_model(variant_key)
    except RuntimeError as e:
        raise HTTPException(status_code=503, detail=str(e))

    transformer = bundle.get('transformer')
    if transformer is None:
        raise HTTPException(503, detail=f"Model '{request.model}' has no transformer — not trained via the new pipeline")

    pipeline = bundle.get('pipeline') or bundle.get('model')
    if pipeline is None:
        raise HTTPException(503, detail=f"Model '{request.model}' has no pipeline")

    df = pd.DataFrame([r.model_dump() for r in request.rows])
    df[request.tax_level] = request.tax_level
    features = transformer.transform(df)
    feat_cols = bundle.get('feature_names', transformer.get_feature_names_out())
    X = features[feat_cols]

    tau = request.target_recall
    xc = request.confidence

    if hasattr(pipeline, 'predict_raw'):
        predicted_cutoff = pipeline.predict(X.values, tau=tau, x_thresh=xc)[0]
    else:
        predicted_cutoff = pipeline.predict(X.values)[0]

    predicted_recall_at_cutoff = 0.0
    if hasattr(pipeline, 'predict_raw') and hasattr(pipeline, 'n_divisions'):
        raw = pipeline.predict_raw(X.values)
        recall_indices = getattr(pipeline, 'recall_indices_', None)
        if recall_indices is not None:
            raw_recalls = raw[0, recall_indices]
        else:
            raw_recalls = raw[0]
        predicted_recall_at_cutoff = _interpolate_recall_at_cutoff(
            predicted_cutoff, raw_recalls, pipeline.n_divisions
        )

    latency_ms = (time.time() - start) * 1000
    log_prediction(
        model_version=f"{registry_name(variant_key)} v{version_info.get('version', '?')}",
        prediction=float(predicted_cutoff),
        latency_ms=latency_ms,
    )

    response.headers["X-Model-Version"] = str(version_info.get("version", ""))
    response.headers["X-Model-Stage"] = str(version_info.get("stage", ""))
    return RecallCutoffResult(
        predicted_cutoff=float(predicted_cutoff),
        predicted_recall_at_cutoff=predicted_recall_at_cutoff,
        target_recall=tau,
        confidence_score=0.9,
        model_version=str(version_info.get("version")),
        model_stage=version_info.get("stage"),
    )


@app.post("/predict_composition_stop_traversal", response_model=CompositionStopTraversalResult)
def predict_composition_stop_traversal(request: CompositionStopTraversalRequest, response: Response):
    start = time.time()
    try:
        bundle, version_info = get_cached_model("composition")
    except RuntimeError as e:
        raise HTTPException(status_code=503, detail=str(e))

    pipeline = bundle.get("pipeline") or bundle.get("model")
    if pipeline is None:
        raise HTTPException(503, detail="Composition model bundle missing 'pipeline' or 'model'")
    fn = bundle.get("feature_names")
    if fn is None and hasattr(pipeline, 'feature_names_in_'):
        fn = pipeline.feature_names_in_
    if fn is None and hasattr(pipeline, 'feature_names_'):
        fn = pipeline.feature_names_
    if isinstance(fn, np.ndarray):
        fn = fn.tolist()

    if fn:
        filled = {**{name: 0.0 for name in fn}, **request.features}
        X = pd.DataFrame([filled])
    else:
        X = pd.DataFrame([request.features])

    pred = pipeline.predict(X)[0]
    proba = pipeline.predict_proba(X)[0]
    confidence_score = float(max(proba))

    latency_ms = (time.time() - start) * 1000
    log_prediction(
        model_version=f"{registry_name('composition')} v{version_info.get('version', '?')}",
        prediction=float(pred),
        latency_ms=latency_ms,
    )

    response.headers["X-Model-Version"] = str(version_info.get("version", ""))
    response.headers["X-Model-Stage"] = str(version_info.get("stage", ""))
    return CompositionStopTraversalResult(
        stop_traversal=bool(pred),
        probability=float(proba[1]),
        confidence_score=confidence_score,
        model_version=str(version_info.get("version")),
        model_stage=version_info.get("stage"),
    )


@app.post("/predict_televir_clustering_threshold", response_model=ClusteringThresholdResult)
def predict_televir_clustering_threshold(request: TelevirClusteringThresholdRequest, response: Response):
    start = time.time()
    model_type = "televir_clustering"
    try:
        model, version_info = get_cached_model(model_type)
    except RuntimeError as e:
        raise HTTPException(status_code=503, detail=str(e))

    X = request.to_array(feature_names=getattr(model, 'feature_names_', None))
    prediction = model.predict([X])[0]
    proba = model.predict_proba([X])[0]
    confidence_score = float(max(proba))
    latency_ms = (time.time() - start) * 1000
    log_prediction(
        model_version=f"{registry_name(model_type)} v{version_info.get('version', '?')}",
        prediction=float(prediction),
        latency_ms=latency_ms,
    )

    response.headers["X-Model-Version"] = str(version_info.get("version", ""))
    response.headers["X-Model-Stage"] = str(version_info.get("stage", ""))
    return ClusteringThresholdResult(
        is_cluster=bool(prediction),
        confidence_score=confidence_score,
        model_version=str(version_info.get("version")),
        model_stage=version_info.get("stage"),
    )
