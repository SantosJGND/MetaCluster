import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import joblib
import pandas as pd
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score
from sklearn.model_selection import train_test_split

from metagenomics_utils.overlap_manager.om_models import ClusteringPipeline

DATA_PATH = Path(__file__).parent / "training_results_cache.parquet"
MODELS_DIR = Path(__file__).parent / "models"
MODEL_SAVE_PATH = MODELS_DIR / "clustering_xgb_bundle.pkl"

TAXON_COLS = [
    "Baculoviridae",
    "Coronaviridae",
    "Flaviviridae",
    "Herelleviridae",
    "Orthoherpesviridae",
    "Papillomaviridae",
    "Peduoviridae",
    "Picornaviridae",
    "Rhabdoviridae",
    "Schitoviridae",
    "Straboviridae",
    "unclassified",
]

NUMERIC_COLS = [
    "taxonomic_diversity",
    "n_leaves",
    "proportion_shared_hits",
    "proportion_unique_hits",
]

FEATURE_NAMES = NUMERIC_COLS + TAXON_COLS


def train():
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(DATA_PATH)
    print(f"Loaded training data: {df.shape}")

    X = pd.DataFrame()
    X["taxonomic_diversity"] = df["tax_diversity"]
    X["n_leaves"] = df["n_leaves"].astype(float)
    X["proportion_shared_hits"] = df["Min_Shared"]
    X["proportion_unique_hits"] = 1.0 - df["Min_Shared"]
    for col in TAXON_COLS:
        X[col] = df[col].astype(float)

    y = df["precision_increased"].astype(int)

    print(f"Feature matrix: {X.shape}")
    print(f"Feature names: {FEATURE_NAMES}")
    print(f"Target distribution:\n{y.value_counts()}")

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        random_state=42,
        stratify=y,
    )

    pipeline = ClusteringPipeline(
        numeric_cols=NUMERIC_COLS,
        taxon_cols=TAXON_COLS,
        feature_names=FEATURE_NAMES,
        use_optuna=True,
    )
    pipeline.fit(X_train, y_train)

    y_pred = pipeline.predict(X_test.values)
    y_prob = pipeline.predict_proba(X_test.values)[:, 1]

    print(f"\nTest accuracy: {accuracy_score(y_test, y_pred):.4f}")
    print(f"Test ROC-AUC: {roc_auc_score(y_test, y_prob):.4f}")
    print("\nClassification report:")
    print(classification_report(y_test, y_pred, target_names=["not_cluster", "is_cluster"]))

    joblib.dump(pipeline, str(MODEL_SAVE_PATH))
    print(f"\nModel saved to {MODEL_SAVE_PATH}")


if __name__ == "__main__":
    train()
