"""
Metrics Debug: Recall & Precision Analysis
Extracted from debug_metrics_test2.ipynb for offline debugging.
"""
import sys, os
project_root = os.path.abspath(os.path.join(os.getcwd(), '..', '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import json
from pathlib import Path

from deployment.model_evaluation.config import EvaluatorConfig, TrainedModels
from deployment.model_evaluation.data_loader import DataLoader
from deployment.model_evaluation.dataset_processor import DatasetProcessor
from deployment.model_evaluation.metrics import (
    compute_purity,
    compute_mstats_precision,
    compute_recall,
    compute_clade_recall,
)

print('Imports OK')

# --- CONFIG ---
STUDY_OUTPUT = Path('/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/virus/output_study')
TAXID_PLAN   = Path('/home/bioinf/Desktop/INSA/Manuscripts/Clustering_Clinical_Metagenomics/data/Panel/viral_assess.tsv')
ANALYSIS_OUT = Path("/home/bioinf/Desktop/INSA/Projectos/CLUSTER_EVAL/study/analysis/virus/filter_xgb_family/")

config = EvaluatorConfig(
    study_output_filepath=STUDY_OUTPUT,
    taxid_plan_filepath=TAXID_PLAN,
    analysis_output_filepath=ANALYSIS_OUT,
    enable_cross_hit=False,
)
print(f'Config OK → tax_level={config.tax_level}')

# --- Load Data ---
loader = DataLoader(config).initialize()
ncbi_wrapper = loader.get_ncbi_wrapper()
input_tax_df = loader.get_input_tax_df()
taxids_to_use = loader.get_taxids_to_use()

# --- Load Models ---
models_dir = config.models_dir
model = 'xgb'

from metagenomics_utils.overlap_manager.om_models import (
    RecallModeller, GPCLFRecallModeller, DirectXGBRecallModeller,
    XGBCompositionModeller, RFCompositionModeller,
    GBCompositionModeller, LRCompositionModeller,
    OptunaXGBCompositionModeller, CrossHitModeller,
)

recall_modeller = None
if models_dir.exists():
    recall_bundle_map = {
        "recall_gp_clf_pipeline.pkl": GPCLFRecallModeller,
        "recall_xgb_bundle.pkl":      RecallModeller,
        "direct_xgb_bundle.pkl":      DirectXGBRecallModeller,
    }
    for fname, cls in recall_bundle_map.items():
        if (models_dir / fname).exists():
            recall_modeller = cls(data_set_divide=config.data_set_divide)
            recall_modeller.load_model(str(models_dir))
            break

    composition_bundle_map = {
        "composition_xgb_bundle.pkl":     XGBCompositionModeller,
        "composition_rf_bundle.pkl":      RFCompositionModeller,
        "composition_gb_bundle.pkl":      GBCompositionModeller,
        "composition_lr_bundle.pkl":      LRCompositionModeller,
        "composition_optuna_bundle.pkl":  OptunaXGBCompositionModeller,
    }
    composition_modeller = None
    for fname, cls in composition_bundle_map.items():
        if (models_dir / fname).exists():
            composition_modeller = cls()
            composition_modeller.load_model(str(models_dir))
            break
    if composition_modeller is None:
        print(f'WARNING: no composition model bundle found in {models_dir}')

    crosshit_modeller = CrossHitModeller(prediction_trainning_results_df=pd.DataFrame(columns=['taxid', 'leaf', 'is_trash']))
    crosshit_modeller.load_model(str(models_dir))

    models = TrainedModels(
        recall_modeller=recall_modeller,
        composition_modeller=composition_modeller,
        crosshit_modeller=crosshit_modeller,
        taxids_to_use=taxids_to_use,
    )
    print(f'Models loaded from {models_dir}')
else:
    print(f'WARNING: models dir not found at {models_dir}')
    models = None

processor = DatasetProcessor(config, models, ncbi_wrapper, input_tax_df, taxids_to_use)
print(f'config.models_dir = {config.models_dir}')

# --- Section 3: Compute m_stats Matrix ---
from metagenomics_utils.overlap_manager.node_stats import get_m_stats_matrix

DATASET_NAME = 'dataset_0005_plan'

overlap_manager, input_summary = processor._load_data(DATASET_NAME)
print(f'OverlapManager leaves: {len(overlap_manager.leaves)}')
print(f'Dataset: {DATASET_NAME}')
print(f'Input taxids: {input_summary["taxid"].nunique()}')

m_stats = get_m_stats_matrix(
    DATASET_NAME,
    str(config.study_output_filepath),
    ncbi_wrapper,
    overlap_manager,
    filter_no_leaf=False,
)

print(f'm_stats shape: {m_stats.shape}')
print(f'Columns: {list(m_stats.columns)}')
print()
print('is_trash distribution:')
print(m_stats['is_trash'].value_counts())
print()
print('best_match_is_best distribution:')
print(m_stats['best_match_is_best'].value_counts())
print()
print('coverage > 0 count:', (m_stats['coverage'] > 0).sum())
print('best_match_taxid non-null:', m_stats['best_match_taxid'].notna().sum())

# --- Section 4: Old vs New Metrics ---
input_taxids = set(input_summary['taxid'].unique())

# Old (buggy) computations — using noutput (count of rows in clean m_stats) instead of intersection
def old_fuzzy_precision(m_stats):
    total = len(m_stats)
    if total == 0: return 0.0, 0.0
    raw = len(m_stats[m_stats['is_trash'] == False]) / total
    cov = len(m_stats[(m_stats['is_trash'] == False) & (m_stats['coverage'] > 0)]) / total
    return raw, cov

def old_mstats_precision(m_stats):
    if len(m_stats) == 0: return 0.0
    u = m_stats.dropna(subset=['best_match_taxid'])
    u = u[(u['is_trash'] == False) & (u['best_match_is_best'] == True)]
    if len(u) == 0: return 0.0
    return u['best_match_taxid'].nunique() / len(u)

def old_recall(m_stats, input_summary):
    unique_taxids = input_summary['taxid'].nunique()
    if unique_taxids == 0: return 0.0, 0.0, [], []
    u = m_stats.dropna(subset=['best_match_taxid'])
    u = u[(u['is_trash'] == False) & (u['best_match_is_best'] == True)]
    recall_raw = u.shape[0] / unique_taxids
    cov = u[u['coverage'] > 0].shape[0] / unique_taxids
    return recall_raw, cov, u['best_match_taxid'].dropna().unique().tolist(), []

old_raw, old_cov = old_fuzzy_precision(m_stats)
old_prec = old_mstats_precision(m_stats)
old_rec, old_rec_cov, _, _ = old_recall(m_stats, input_summary)

print('\n=== OLD (buggy) metrics ===')
print(f'  Purity (fuzzy):         raw={old_raw:.4f}  cov_filtered={old_cov:.4f}')
print(f'  Overall precision:      {old_prec:.4f}')
print(f'  Recall:                 raw={old_rec:.4f}  cov_filtered={old_rec_cov:.4f}')

# New (fixed) computations
new_purity_raw, purity_best_match, new_purity_cov = compute_purity(m_stats)
new_prec = compute_mstats_precision(m_stats, input_summary)
new_rec, new_rec_cov, found_taxids, found_taxids_cov = compute_recall(m_stats, input_summary)

print()
print('=== NEW (fixed) metrics ===')
print(f'  Purity:                 raw={new_purity_raw:.4f}  cov_filtered={new_purity_cov:.4f}')
print(f'  Overall precision:      {new_prec:.4f}')
print(f'  Recall:                 raw={new_rec:.4f}  cov_filtered={new_rec_cov:.4f}')
print()
print('=== Summary ===')
print(f'  Input taxids: {len(input_taxids)}')
print(f'  Precision Δ: {new_prec - old_prec:+.4f}')
print(f'  Recall Δ:    {new_rec - old_rec:+.4f}')

# --- Section 5: Set Intersection Analysis ---
best_matches = m_stats.dropna(subset=['best_match_taxid'])
best_matches = best_matches[(best_matches['is_trash'] == False) & (best_matches['best_match_is_best'] == True)]
output_taxids = set(best_matches['best_match_taxid'].unique())

tp = input_taxids & output_taxids
fp = output_taxids - input_taxids
fn = input_taxids - output_taxids

print(f'\nTrue Positives  (in input & output): {len(tp)} taxids')
print(f'False Positives (in output only):    {len(fp)} taxids')
print(f'False Negatives (in input only):     {len(fn)} taxids')
if fn:
    print('Missed taxids:      ', sorted(fn)[:20])

expected_recall = len(tp) / len(input_taxids) if len(input_taxids) > 0 else 0
expected_prec  = len(tp) / len(output_taxids) if len(output_taxids) > 0 else 0
print(f'  Expected recall (TP/|input|):      {expected_recall:.4f}  (matches new: {abs(expected_recall - new_rec) < 1e-10})')
print(f'  Expected precision (TP/|output|):  {expected_prec:.4f}   (matches new: {abs(expected_prec - new_prec) < 1e-10})')

# --- Section 6: Filter Effects ---
def metrics_row(label, m_stats, input_summary):
    p_raw, r_raw_best_match, p_cov = compute_purity(m_stats)
    prec = compute_mstats_precision(m_stats, input_summary)
    rec, rec_cov, _, _ = compute_recall(m_stats, input_summary)
    return {
        'step': label,
        'n_leaves': len(m_stats),
        'purity_raw': round(p_raw, 4),
        'purity_cov': round(p_cov, 4),
        'precision': round(prec, 4),
        'recall_raw': round(rec, 4),
        'recall_cov': round(rec_cov, 4),
    }

rows = [metrics_row('baseline', m_stats, input_summary)]

if models is not None:
    from metagenomics_utils.overlap_manager.om_models import cut_off_recall_prediction
    print(f'\nRecall modeller: {models.recall_modeller}')
    try:
        filtered_om, metrics_dict = cut_off_recall_prediction(
            str(config.study_output_filepath),
            DATASET_NAME, models.recall_modeller,
            config.data_set_divide, m_stats,
            taxids_to_use, tax_level=config.tax_level,
            target_recall=config.target_recall,
            confidence=config.cutoff_confidence,
        )
        new_m_stats = get_m_stats_matrix(DATASET_NAME, str(config.study_output_filepath), ncbi_wrapper, filtered_om, filter_no_leaf=False)
        new_m_stats = new_m_stats.head(metrics_dict.get('keep_index', len(new_m_stats)))
        new_m_stats = new_m_stats[new_m_stats['coverage'] > 0]
        rows.append(metrics_row('after_recall_filter', new_m_stats, input_summary))
    except Exception as e:
        print(f'  Recall filter failed: {e}')

print('\nFilter effects:')
print(pd.DataFrame(rows).to_string(index=False))

# --- Section 7: Clade Precision ---
from metagenomics_utils.overlap_manager.om_models import (
    predict_data_set_clades_composition,
    predict_data_set_clades_fixed,
    calculate_clade_precision,
)

m_stats_leaf = get_m_stats_matrix(DATASET_NAME, str(config.study_output_filepath), ncbi_wrapper, overlap_manager, filter_no_leaf=True)
print(f'\nm_stats_leaf shape: {m_stats_leaf.shape}')

print('\nInput summary:')
print(input_summary.to_string())

# Selected method (composition model)
results_df = None
if models is not None:
    try:
        results_df = predict_data_set_clades_composition(
            DATASET_NAME, m_stats_leaf, overlap_manager,
            models.composition_modeller, taxids_to_use, tax_level=config.tax_level
        )
        clade_prec = calculate_clade_precision(results_df, input_summary)
        clade_rec = compute_clade_recall(results_df, input_summary)
        print(f'\nSelected method:  {len(results_df)} clades, precision={clade_prec:.4f}, recall={clade_rec:.4f}')
    except Exception as e:
        print(f'\nSelected method FAILED: {e}')
        import traceback
        traceback.print_exc()
else:
    print('\nNo models loaded — skipping selected method')

# Fixed method
fixed_df = predict_data_set_clades_fixed(DATASET_NAME, m_stats_leaf, overlap_manager, min_dist_threshold=0.6)
clade_prec_fixed = calculate_clade_precision(fixed_df, input_summary)
print(f'Fixed method:     {len(fixed_df)} clades, precision={clade_prec_fixed:.4f}')

# Comparison
input_t = set(input_summary['taxid'].unique())
if results_df is not None:
    predicted_s = set(results_df['best_taxid_match'].dropna().unique())
else:
    predicted_s = set()
predicted_f = set(fixed_df['best_taxid_match'].dropna().unique())

for label, predicted in [('Selected', predicted_s), ('Fixed', predicted_f)]:
    tp = len(input_t & predicted)
    fp = len(predicted - input_t)
    fn = len(input_t - predicted)
    print(f'{label}: TP={tp}  FP={fp}  FN={fn}')

# Save comparison figure
fig, ax = plt.subplots(figsize=(8, 4))
metrics_names = ['Selected\nPrecision', 'Fixed\nPrecision', 'Clade\nRecall']
vals = [clade_prec if results_df is not None else 0.0,
        clade_prec_fixed,
        clade_rec if results_df is not None else 0.0]
bars = ax.bar(range(len(metrics_names)), vals, color=['steelblue', 'seagreen', 'coral'], alpha=0.8)
ax.set_ylabel('Score')
ax.set_title(f'Clade Precision & Recall — {DATASET_NAME}')
ax.set_xticks(range(len(metrics_names)))
ax.set_xticklabels(metrics_names)
ax.set_ylim(0, 1.1)
ax.grid(axis='y', alpha=0.3)
for bar in bars:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=10)
plt.tight_layout()
plt.savefig('clade_precision_comparison.png', dpi=150)
print('\nSaved clade_precision_comparison.png')

# --- Section 8: Visual Comparison ---
metrics_names2 = ['Purity', 'Precision', 'Recall']
old_vals = [old_raw, old_prec, old_rec]
new_vals = [new_purity_raw, new_prec, new_rec]

x = np.arange(len(metrics_names2))
width = 0.35

fig, ax = plt.subplots(figsize=(10, 5))
bars1 = ax.bar(x - width/2, old_vals, width, label='Old (buggy)', alpha=0.8)
bars2 = ax.bar(x + width/2, new_vals, width, label='New (fixed)', alpha=0.8)
ax.set_ylabel('Score')
ax.set_title(f'Old vs New Metrics — {DATASET_NAME}')
ax.set_xticks(x)
ax.set_xticklabels(metrics_names2)
ax.legend()
ax.grid(axis='y', alpha=0.3)
for bar in bars1:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=9)
for bar in bars2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
            f'{bar.get_height():.3f}', ha='center', va='bottom', fontsize=9)
plt.tight_layout()
plt.savefig('old_vs_new_metrics.png', dpi=150)
print('Saved old_vs_new_metrics.png')

print('\n=== DONE ===')
print('Key findings:')
print(f'  1. Old precision stuck at {old_prec:.4f} (always 1.0 due to best_match_is_best grouping)')
print(f'  2. New precision = {new_prec:.4f}  (meaningful: fraction of output taxids in input)')
print(f'  3. Old recall = {old_rec:.4f} (may overcount by including hallucinated taxids)')
print(f'  4. New recall = {new_rec:.4f}  (correct: fraction of input taxids recovered)')

# --- Section 9: Batch Clade Precision Diagnosis ---
print('\n=== Section 9: Batch Clade Precision Diagnosis ===')
from metagenomics_utils.overlap_manager.om_models import (
    predict_data_set_clades_composition,
    predict_data_set_clades_fixed,
    calculate_clade_precision,
)
from deployment.model_evaluation.metrics import compute_clade_recall
from tqdm import tqdm
import traceback
from deployment.model_evaluation.result_models import DatasetResult

all_datasets = sorted([
    d for d in os.listdir(config.study_output_filepath)
    if os.path.isdir(os.path.join(config.study_output_filepath, d))
])
print(f'Found {len(all_datasets)} datasets')

records = []
MAX_DS = min(200, len(all_datasets))

#for ds in tqdm(all_datasets[:MAX_DS], desc='Batch clade precision'):
for ds in ['dataset_0071_plan']:
    try:
        om, inp = processor._load_data(ds)
        if om is None: 
            continue
        m_stats = get_m_stats_matrix(
            ds, str(config.study_output_filepath),
            ncbi_wrapper, om, filter_no_leaf=True
        )
        n_leaves = len(om.leaves)
        n_input = inp['taxid'].nunique()

        # FULL (pre-cleanup, composition model)
        try:
            full_df = predict_data_set_clades_composition(
                ds, m_stats, om,
                models.composition_modeller, taxids_to_use,
                tax_level=config.tax_level
            )

            print("pre-cleanup clades predicted: ", len(full_df))

            prec_full = calculate_clade_precision(full_df, inp)
            if prec_full == 0.0:

                print(om.distance_mat)
                print(om.root_nodes)
                print(om.leaves)
                print("##################3")
                m_stats = get_m_stats_matrix(
                    ds, str(config.study_output_filepath),
                    ncbi_wrapper, om, filter_no_leaf=False, log = True
                )
                full_df = predict_data_set_clades_composition(
                    ds, m_stats, om,
                    models.composition_modeller, taxids_to_use,
                    tax_level=config.tax_level
                )
                print("pre-cleanup clades predicted: ", len(full_df))

                prec_full = calculate_clade_precision(full_df, inp)
                print(om.m_stats_matrix.head())
                print(m_stats.head())
                print(f'  Dataset {ds}: {len(full_df)} clades predicted (pre-cleanup)')
                print(full_df.head())
                print(f'  Dataset {ds}: precision={prec_full:.4f} (pre-cleanup)')
            rec_full = compute_clade_recall(full_df, inp)
            n_full = len(full_df.drop_duplicates(subset=['node']))
            n_pred_full = full_df['best_taxid_match'].dropna().nunique()
            failure_full = 'A' if n_pred_full == 0 else ('B' if prec_full == 0 else 'OK')
        except Exception:
            prec_full, rec_full, n_full, n_pred_full = 0.0, 0.0, 0, 0
            failure_full = 'C'

        # POST (after cross-hit cleanup)
        try:
            result2 = processor._apply_crosshit_cleanup(ds, om, DatasetResult(ds))
            om2 = result2.overlap_manager if hasattr(result2, 'overlap_manager') else om
            
            post_df = predict_data_set_clades_composition(
                ds, m_stats, om2,
                models.composition_modeller, taxids_to_use,
                tax_level=config.tax_level
            )
            prec_post = calculate_clade_precision(post_df, inp)
            n_post = len(post_df.drop_duplicates(subset=['node']))
            n_pred_post = post_df['best_taxid_match'].dropna().nunique()
            failure_post = 'A' if n_pred_post == 0 else ('B' if prec_post == 0 else 'OK')
        except Exception:
            prec_post, n_post, n_pred_post = 0.0, 0, 0
            failure_post = 'C'

        # FIXED (post-cleanup, distance threshold)
        try:
            result2 = processor._apply_crosshit_cleanup(ds, om, DatasetResult(ds))
            om2 = result2.overlap_manager if hasattr(result2, 'overlap_manager') else om
            

            fixed_df = predict_data_set_clades_fixed(
                ds, m_stats, om2,
                min_dist_threshold=0.6
            )
            prec_fixed = calculate_clade_precision(fixed_df, inp)
            if prec_fixed == 0.0:
                print("##################3")
                print(f'  Dataset {ds}: {len(fixed_df)} clades predicted (post-cleanup, fixed)')
                print(fixed_df.head())
                print(f'  Dataset {ds}: precision={prec_fixed:.4f} (post-cleanup, fixed)')
            n_fixed = len(fixed_df.drop_duplicates(subset=['node']))
            n_pred_fixed = fixed_df['best_taxid_match'].dropna().nunique()
            failure_fixed = 'A' if n_pred_fixed == 0 else ('B' if prec_fixed == 0 else 'OK')
        except Exception as e:
            print(f'  Dataset {ds}: fixed method failed with exception')
            print(e)
            prec_fixed, n_fixed, n_pred_fixed = 0.0, 0, 0
            failure_fixed = 'C'

        records.append({
            'dataset': ds, 'n_leaves': n_leaves, 'n_input': n_input,
            'prec_full': prec_full, 'rec_full': rec_full,
            'n_nodes_full': n_full, 'n_pred_full': n_pred_full,
            'failure_full': failure_full,
            'prec_post': prec_post, 'n_nodes_post': n_post,
            'n_pred_post': n_pred_post, 'failure_post': failure_post,
            'prec_fixed': prec_fixed, 'n_nodes_fixed': n_fixed,
            'n_pred_fixed': n_pred_fixed, 'failure_fixed': failure_fixed,
        })
    except Exception as e:
        traceback.print_exc()
        continue

diag = pd.DataFrame(records)
print(f'Processed {len(diag)} datasets')

print('\n=== Clade Precision: Zero proportion by method ===')
for method, col in [('full', 'prec_full'), ('post', 'prec_post'), ('fixed', 'prec_fixed')]:
    zero_count = (diag[col] == 0.0).sum()
    total = len(diag)
    print(f'  {method:6s}: {zero_count}/{total} = {zero_count/total*100:.1f}% precision=0')

print('\n=== Failure mode breakdown by method ===')
for method, fcol in [('full', 'failure_full'), ('post', 'failure_post'), ('fixed', 'failure_fixed')]:
    print(f'\n{method.upper()}:')
    counts = diag[fcol].value_counts()
    for k in ['OK', 'A', 'B', 'C']:
        v = counts.get(k, 0)
        lbl = {'OK': '>0 precision', 'A': 'empty predicted set',
               'B': 'no intersection', 'C': 'exception'}[k]
        print(f'    {lbl:20s}: {v}')

print('\n=== Transition: full→post precision delta ===')
dropped = (diag['prec_full'] > 0) & (diag['prec_post'] == 0)
print(f'  Datasets where full>0 but post=0: {dropped.sum()}')
improved = (diag['prec_full'] == 0) & (diag['prec_post'] > 0)
print(f'  Datasets where full=0 but post>0: {improved.sum()}')

# Diagnostic plots
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle('Batch Clade Precision Diagnosis', fontsize=16)

for idx, (method, col) in enumerate([('full', 'prec_full'), ('post', 'prec_post'), ('fixed', 'prec_fixed')]):
    ax = axes[0, idx]
    vals = diag[col].values
    ax.hist(vals[vals > 0], bins=20, color='steelblue', alpha=0.7, edgecolor='white')
    ax.axvline(0, color='red', linestyle='--', linewidth=1)
    n_zero = (vals == 0).sum()
    ax.set_title(f'{method} (n_zero={n_zero})')
    ax.set_xlabel('Precision')
    ax.set_ylabel('Datasets')
    ax.set_yscale('log')

for idx, (method, fcol) in enumerate([('full', 'failure_full'), ('post', 'failure_post'), ('fixed', 'failure_fixed')]):
    ax = axes[1, idx]
    counts = diag[fcol].value_counts()
    colors = {'OK': 'steelblue', 'A': 'coral', 'B': 'orange', 'C': 'grey'}
    bar_items = ax.bar(
        [c for c in ['OK', 'A', 'B', 'C'] if c in counts.index],
        [counts.get(c, 0) for c in ['OK', 'A', 'B', 'C'] if c in counts.index],
        color=[colors[c] for c in ['OK', 'A', 'B', 'C'] if c in counts.index],
        edgecolor='white'
    )
    ax.set_title(f'{method} failures')
    ax.set_ylabel('Datasets')
    for b in bar_items:
        ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.5,
                str(int(b.get_height())), ha='center', fontsize=10)

plt.tight_layout()
plt.savefig('batch_clade_precision_diagnosis.png', dpi=150)
print('Saved batch_clade_precision_diagnosis.png')

# Save diagnostics to CSV
diag.to_csv('clade_precision_diagnostics.csv', index=False)
print('Saved clade_precision_diagnostics.csv')

print('\n=== ALL DONE ===')
