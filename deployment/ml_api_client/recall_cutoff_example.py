"""
Client-side usage examples for the INSaFLU ML API.

Prerequisites:
    pip install -r deployment/ml_api_client/requirements.txt

Run:
    python deployment/ml_api_client/recall_cutoff_example.py [--url URL]
"""
import argparse

from ml_api_client import MLAPIClient


def main():
    parser = argparse.ArgumentParser(description="INSaFLU ML API client examples")
    parser.add_argument("--url", default="http://localhost:8000", help="API base URL")
    args = parser.parse_args()

    c = MLAPIClient(base_url=args.url)

    # --- 1. Health check ---
    print("=== health ===")
    print(c.health())

    # --- 2. List loaded models ---
    print("\n=== models ===")
    models_resp = c.models()
    for m in models_resp.get("models", []):
        print(f"  {m['model_type']:40s}  v{m['version']}  ({m['stage']})")

    # --- 3. Recall cutoff — minimal (default model) ---
    print("\n=== predict_recall_cutoff (default model) ===")
rows = [
    {"taxid": 2697049, "total_uniq_reads": 2348731,
        "order": "Mononegavirales", "family": "Filoviridae", "best_match_is_best": True},
    {"taxid": 1122334, "total_uniq_reads": 295652,
        "order": "Picornaviridae", "family": "Picornaviridae", "best_match_is_best": True},
]
    r = c.predict_recall_cutoff(rows, model="order_recall_direct_xgb")
    print(f"  predicted_cutoff:            {r['predicted_cutoff']:.4f}")
    print(f"  predicted_recall_at_cutoff:  {r['predicted_recall_at_cutoff']:.4f}")
    print(f"  model_version:               {r['model_version']}")

    # --- 4. Recall cutoff — with target_recall and confidence ---
    print("\n=== predict_recall_cutoff (tau=0.95, X=0.8) ===")
    r = c.predict_recall_cutoff(rows, target_recall=0.95, confidence=0.8)
    print(f"  predicted_cutoff:  {r['predicted_cutoff']:.4f}")
    print(f"  target_recall:     {r['target_recall']}")

    # --- 5. Recall cutoff — different variant ---
    print("\n=== predict_recall_cutoff (xgb_direct) ===")
    r = c.predict_recall_cutoff(rows, model="order_recall_xgb_direct")
    print(f"  predicted_cutoff:  {r['predicted_cutoff']:.4f}")

    # --- 6. Composition stop-traversal ---
    print("\n=== predict_composition_stop_traversal ===")
    node_features = {
        "n_leaves": 12,
        "tax_diversity": 1.2,
        "Min_Dist": 0.15,
        "Min_Shared": 0.3,
        "Coronaviridae": 0.021,
        "Peduoviridae": 0.0,
        "Straboviridae": 0.0,
    }
    r = c.predict_composition_stop_traversal(node_features, model="xgb")
    print(f"  stop_traversal:  {r['stop_traversal']}")
    print(f"  probability:     {r['probability']:.4f}")

    # --- 7. Reload a single model ---
    print("\n=== reload order_recall_gp_clf ===")
    print(c.reload("order_recall_gp_clf"))

    print("\nDone.")


if __name__ == "__main__":
    main()
