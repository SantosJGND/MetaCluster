import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import pandas as pd

from metagenomics_utils.overlap_manager import OverlapManager
from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix
from metagenomics_utils.overlap_manager.om_models import (
    DirectXGBRecallModeller,
    GPCLFRecallModeller,
    RecallModeller,
)

MODEL_CLASSES = {
    "gp_clf": GPCLFRecallModeller,
    "xgb_direct": DirectXGBRecallModeller,
    "xgb_multi": RecallModeller,
}


def train():
    parser = argparse.ArgumentParser(description="Train a recall model variant")
    parser.add_argument("--study-dir", required=True, help="Study output directory")
    parser.add_argument("--folders", nargs="+", required=True, help="Training folder names")
    parser.add_argument("--taxids-file", required=True, help="Path to taxids_to_use.parquet")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: <study-dir>/models)")
    parser.add_argument("--data-set-divide", type=int, default=20, help="Number of recall divisions")
    parser.add_argument("--model", default="gp_clf", choices=list(MODEL_CLASSES), help="Model variant to train")
    parser.add_argument("--tax-level", default="order", help="Taxonomic level used for training")
    parser.add_argument("--description", default=None, help="Optional description for the trained model")

    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else Path(args.study_dir) / "models"
    output_dir.mkdir(parents=True, exist_ok=True)

    taxids_to_use = pd.read_parquet(args.taxids_file)

    recall_matrices = []
    for folder in args.folders:
        clustering_dir = os.path.join(args.study_dir, folder, "clustering")
        if not os.path.isdir(clustering_dir):
            print(f"Skipping {folder}: clustering dir not found")
            continue
        om = OverlapManager(clustering_dir)
        if om.m_stats_matrix.empty:
            print(f"Skipping {folder}: no mapped reads")
            continue
        matrix = get_m_stats_matrix(folder, args.study_dir, None, om)
        if not matrix.empty:
            recall_matrices.append(matrix)
            print(f"  {folder}: {matrix.shape[0]} rows")

    if not recall_matrices:
        print("ERROR: No training data collected.")
        sys.exit(1)

    print(f"Collected {len(recall_matrices)} raw recall matrices")

    cls = MODEL_CLASSES[args.model]
    modeller = cls(
        data_set_divide=args.data_set_divide,
        tax_level=args.tax_level,
        description=args.description,
    )

    print(f"Training {args.model} recall model ...")
    modeller.fit(recall_matrices, taxids_to_use)
    modeller.save_model(str(output_dir))

    bundle_path = output_dir / modeller.model_save_filename
    print(f"Model saved to {bundle_path}")
    print(f"  model_type: {args.model}")
    print(f"  feature names ({len(modeller.RecP_feature_cols)}):")
    for fn in modeller.RecP_feature_cols:
        print(f"    - {fn}")

    return modeller


if __name__ == "__main__":
    train()
