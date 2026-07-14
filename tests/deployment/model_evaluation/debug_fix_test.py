"""
Debug: verify clade prediction fix on single-accession (cov_filtered=1) datasets.
Requires study output and trained models from a prior run.
"""

import os
import sys

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

from pathlib import Path

import pandas as pd

from deployment.model_evaluation.config import EvaluatorConfig, TrainedModels
from deployment.model_evaluation.data_loader import DataLoader
from deployment.model_evaluation.dataset_processor import DatasetProcessor
from metagenomics_utils.overlap_manager.om_models import (
    CrossHitModeller,
    DirectXGBRecallModeller,
    GBCompositionModeller,
    GPCLFRecallModeller,
    LRCompositionModeller,
    OptunaXGBCompositionModeller,
    RFCompositionModeller,
    RecallModeller,
    XGBCompositionModeller,
)

# --- CONFIG ---
STUDY_OUTPUT = Path("/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/virus/output_study")
TAXID_PLAN = Path("/home/bioinf/Desktop/INSA/Manuscripts/Clustering_Clinical_Metagenomics/data/Panel/viral_assess.tsv")
MODELS_DIR = Path(
    "/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/virus/filter_direct_xgb_order/models"
)

config = EvaluatorConfig(
    study_output_filepath=STUDY_OUTPUT,
    taxid_plan_filepath=TAXID_PLAN,
    analysis_output_filepath=MODELS_DIR.parent,
    enable_cross_hit=False,
)

# --- Load Data ---
loader = DataLoader(config).initialize()
ncbi_wrapper = loader.get_ncbi_wrapper()
input_tax_df = loader.get_input_tax_df()
taxids_to_use = loader.get_taxids_to_use()

# --- Load Models ---
recall_bundle_map = {
    "direct_xgb_bundle.pkl": DirectXGBRecallModeller,
    "recall_xgb_bundle.pkl": RecallModeller,
    "recall_gp_clf_pipeline.pkl": GPCLFRecallModeller,
}
recall_modeller = None
for fname, cls in recall_bundle_map.items():
    if (MODELS_DIR / fname).exists():
        recall_modeller = cls(data_set_divide=config.data_set_divide)
        recall_modeller.load_model(str(MODELS_DIR))
        break

composition_bundle_map = {
    "composition_rf_bundle.pkl": RFCompositionModeller,
    "composition_xgb_bundle.pkl": XGBCompositionModeller,
    "composition_gb_bundle.pkl": GBCompositionModeller,
    "composition_lr_bundle.pkl": LRCompositionModeller,
    "composition_optuna_bundle.pkl": OptunaXGBCompositionModeller,
}
composition_modeller = None
for fname, cls in composition_bundle_map.items():
    if (MODELS_DIR / fname).exists():
        composition_modeller = cls()
        composition_modeller.load_model(str(MODELS_DIR))
        break

crosshit_modeller = CrossHitModeller(
    prediction_trainning_results_df=pd.DataFrame(columns=["taxid", "leaf", "is_trash"])
)
crosshit_modeller.load_model(str(MODELS_DIR))

models = TrainedModels(
    recall_modeller=recall_modeller,
    composition_modeller=composition_modeller,
    crosshit_modeller=crosshit_modeller,
    taxids_to_use=taxids_to_use,
)

processor = DatasetProcessor(config, models, ncbi_wrapper, input_tax_df, taxids_to_use)

# --- Test datasets with output_cov_filtered==1 from evaluation_results.json ---
DATASETS = [
    "dataset_0505_plan",
    "dataset_0690_plan",
    "dataset_1355_plan",
    "dataset_1475_plan",
    "dataset_1983_plan",
    "dataset_0370_plan",
]

print(f"Testing {len(DATASETS)} single-accession datasets...")
all_ok = True
for ds in DATASETS:
    result = processor.process(ds)
    if result is None:
        print(f"  {ds}: SKIPPED")
        all_ok = False
        continue

    fixed = result.predicted_clades_fixed
    post = result.predicted_clades_post
    prec_fixed = result.precision.clade_precision_fixed
    prec_post = result.precision.clade_precision_post

    ok = fixed >= 1 and post >= 1
    status = "OK" if ok else "FAIL"
    print(f"  {ds}: fixed={fixed} ({prec_fixed:.2f}) post={post} ({prec_post:.2f}) [{status}]")
    if not ok:
        all_ok = False

print()
if all_ok:
    print("All datasets produce at least 1 clade with this fix.")
else:
    print("Some datasets still have 0 predicted clades.")
