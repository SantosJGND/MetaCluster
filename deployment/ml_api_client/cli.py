import argparse
import csv
import json
import sys

from .client import MLAPIError, MLAPIClient


def _read_rows(csv_path: str) -> list[dict]:
    rows = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                {
                    "taxid": int(row["taxid"]),
                    "total_uniq_reads": float(row["total_uniq_reads"]),
                    "order": row["order"],
                    "family": row["family"],
                }
            )
    return rows


def _read_json(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser(prog="mlapi", description="INSaFLU ML API client")
    parser.add_argument("--base-url", default="http://localhost:8000", help="API base URL")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("health", help="Liveness check")

    sub.add_parser("models", help="List cached models")

    reload_p = sub.add_parser("reload", help="Reload model cache")
    reload_p.add_argument("model_type", nargs="?", help="Specific model to reload (omit for all)")

    comp_tax_p = sub.add_parser("composition-tax-level", help="Get tax_level of composition model(s)")
    comp_tax_p.add_argument("--model", help="Full composite key or short variant (e.g. order_composition_rf or rf)")

    recall_p = sub.add_parser("predict-recall", help="Predict recall cutoff from raw table")
    recall_p.add_argument("--rows", required=True, help="CSV with taxid,total_uniq_reads,order,family")
    recall_p.add_argument("--model", default="gp_clf", choices=["gp_clf", "xgb_direct", "xgb_multi"],
                          help="Recall model variant (default: gp_clf)")
    recall_p.add_argument("--target-recall", type=float,
                          help="Target recall threshold (0-1)")
    recall_p.add_argument("--confidence", type=float,
                          help="Confidence threshold for probability-guided cutoff (0-1)")

    clust_p = sub.add_parser("predict-clustering", help="Predict TELEVIR clustering threshold")
    clust_p.add_argument("--features", required=True, help="JSON file with clustering features")

    comp_p = sub.add_parser("predict-composition-stop", help="Predict stop_traversal from node features")
    comp_p.add_argument("--features", required=True, help="JSON file with node feature dict")
    comp_p.add_argument("--model", default="order_composition_xgb",
                        help="Composition model composite key (default: order_composition_xgb)")

    args = parser.parse_args()
    client = MLAPIClient(base_url=args.base_url)

    try:
        if args.command == "health":
            result = client.health()
        elif args.command == "models":
            result = client.models()
        elif args.command == "reload":
            result = client.reload(args.model_type)
        elif args.command == "composition-tax-level":
            result = client.composition_model_tax_level(args.model)
        elif args.command == "predict-recall":
            rows = _read_rows(args.rows)
            result = client.predict_recall_cutoff(
                rows,
                model=args.model,
                target_recall=args.target_recall,
                confidence=args.confidence,
            )
        elif args.command == "predict-clustering":
            features = _read_json(args.features)
            result = client.predict_clustering_threshold(features)
        elif args.command == "predict-composition-stop":
            features = _read_json(args.features)
            result = client.predict_composition_stop_traversal(features, model=args.model)
        else:
            parser.print_help()
            sys.exit(1)
    except MLAPIError as e:
        print(json.dumps({"error": e.detail, "status_code": e.status_code}, default=str), file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(json.dumps({"error": str(e)}), file=sys.stderr)
        sys.exit(1)

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
