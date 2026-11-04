#!/usr/bin/env python3
import argparse
import csv
import json
import logging
import os
import re
import statistics
import sys
from collections import defaultdict

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except Exception:
    HAS_MPL = False

HERE = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.dirname(HERE)
TAXA_DB = os.path.join(PIPELINE_DIR, 'taxa.db')
LOG_NAME = 'esviritu_sim_evaluation.log'
IDENTIFIER_CACHE_NAME = 'esviritu_identifier_taxid_cache.json'

logger = logging.getLogger('evaluate_esviritu_sim')


def norm_taxid(v):
    if v is None: return None
    s = str(v).strip()
    if not s or s.lower() in {'nan','none','null','na'}: return None
    m = re.findall(r'\d+', s)
    return m[0] if m else s


def norm_species(v):
    """Strip EsViritu rank prefix (s__) and normalise case/whitespace."""
    if v is None: return None
    s = re.sub(r'^[a-z]__', '', str(v).strip()).strip().lower()
    if s in {'', 'nan', 'none', 'null', 'na'}: return None
    return s


GENERIC_TAXONOMY_NAMES = {
    'viruses', 'unclassified viruses', 'root', 'cellular organisms', 'unclassified'
}


def is_generic_taxonomy_name(v):
    """True for placeholder labels that cannot support taxonomic matching."""
    if v is None: return True
    n = re.sub(r'\s+', ' ', str(v).strip().replace('_', ' ')).strip().lower()
    return n in GENERIC_TAXONOMY_NAMES or n.startswith('unclassified')


def norm_rank(v):
    """Strip rank prefix, normalise separators, and reject generic labels."""
    if v is None: return None
    s = re.sub(r'^[a-z]__', '', str(v).strip()).strip()
    s = re.sub(r'\s+', ' ', s.replace('_', ' ')).strip()
    if s.lower() in {'', 'nan', 'none', 'null', 'na'}: return None
    if is_generic_taxonomy_name(s): return None
    return s


def to_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def read_table(path):
    with open(path, newline='', encoding='utf-8-sig') as fh:
        first = fh.readline()
        delim = '\t' if first.count('\t')>=first.count(',') else ','
        fh.seek(0)
        return list(csv.DictReader(fh, delimiter=delim))

def load_gt_taxids(sim_root, dataset):
    gt = set()
    gt_info = {}
    input_path = os.path.join(sim_root, dataset, 'input', f'{dataset}.tsv')
    if not os.path.exists(input_path):
        logger.warning(f'{dataset}: ground-truth file not found: {input_path}')
        return gt, gt_info
    for row in read_table(input_path):
        taxid = norm_taxid(row.get('taxid') or row.get('TaxID') or row.get('taxID') or row.get('taxon'))
        if not taxid:
            # try parse from path
            p = row.get('path') or row.get('fasta') or row.get('assembly_file') or ''
            m = re.search(r'(\d{4,})_', os.path.basename(str(p)))
            if m: taxid = m.group(1)
        if taxid:
            gt.add(taxid)
            gt_info[taxid] = {'taxid': taxid}
    return gt, gt_info


def load_esviritu_rows(esviritu_dir, dataset):
    """One dict per EsViritu assembly row, keeping identifiers for NCBI lookup."""
    path = os.path.join(esviritu_dir, f'{dataset}.detected_virus.info.tsv')
    if not os.path.exists(path):
        logger.warning(f'{dataset}: EsViritu table not found: {path}')
        return []
    rows = []
    for row in read_table(path):
        rows.append({
            'assembly': str(row.get('Assembly') or '').strip(),
            'accession': str(row.get('Accession') or '').strip(),
            'read_count': to_float(row.get('read_count')),
            'covered_bases': to_float(row.get('covered_bases')),
            'mean_coverage': to_float(row.get('mean_coverage')),
            'esv_species': norm_species(row.get('species') or row.get('Species')),
            'esv_order': norm_rank(row.get('order') or row.get('Order')),
        })
    return rows


def load_identifier_cache(path):
    """Load the persistent accession→taxid cache, ignoring malformed files."""
    if not path or not os.path.exists(path):
        return {}
    try:
        with open(path, encoding='utf-8') as fh:
            data = json.load(fh)
        if not isinstance(data, dict):
            raise ValueError('cache root is not an object')
        return {str(k): int(v) for k, v in data.items() if v is not None}
    except Exception as e:
        logger.warning(f'ignoring identifier cache {path}: {e}')
        return {}


def save_identifier_cache(path, cache):
    """Atomically write successful accession→taxid lookups to disk."""
    if not path:
        return
    try:
        parent = os.path.dirname(path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        tmp = f'{path}.tmp'
        with open(tmp, 'w', encoding='utf-8') as fh:
            json.dump(cache, fh, indent=2, sort_keys=True)
        os.replace(tmp, path)
    except Exception as e:
        logger.warning(f'could not write identifier cache {path}: {e}')


def resolve_prediction_taxids(rows, cache_path=None):
    """Attach an NCBI taxid to each row via Accession, then Assembly.

    Species and order are derived from the resolved identifier through NCBI
    taxonomy. The EsViritu ``species``/``order`` columns are retained as
    fallbacks for rows whose identifier cannot be resolved. Unresolved
    identifiers are logged and leave ``pred_taxid`` as None. Successful
    identifier→taxid lookups are persisted to ``cache_path`` so repeated runs
    skip Entrez.
    """
    if not rows:
        return rows
    if PIPELINE_DIR not in sys.path:
        sys.path.insert(0, PIPELINE_DIR)
    try:
        from metagenomics_utils.ncbi_tools import normalize_accession, retrieve_taxid_from_accession
    except Exception as e:
        logger.error(f'NCBI tools unavailable ({e}); no prediction identifier can be resolved')
        for row in rows:
            row['pred_taxid'] = None
        return rows

    cache = load_identifier_cache(cache_path)
    misses = set()
    unresolved = set()
    cache_hits = 0
    dirty = False
    for row in rows:
        row['pred_taxid'] = None
        candidates = []
        for value in (row.get('accession'), row.get('assembly')):
            if not value:
                continue
            accession, db = normalize_accession(value)
            if accession and (accession, db) not in candidates:
                candidates.append((accession, db))
        for accession, db in candidates:
            key = f'{db}:{accession}'
            if key in cache:
                cache_hits += 1
                taxid = cache[key]
            elif key in misses:
                taxid = None
            else:
                try:
                    taxid = retrieve_taxid_from_accession(accession, db)
                except Exception as e:
                    logger.error(f'identifier lookup failed for {accession} ({db}): {e}')
                    taxid = None
                if taxid is None:
                    misses.add(key)
                else:
                    cache[key] = int(taxid)
                    dirty = True
            if taxid is not None:
                row['pred_taxid'] = str(taxid)
                break
        if row['pred_taxid'] is None:
            label = row.get('accession') or row.get('assembly') or '<no identifier>'
            if label not in unresolved:
                unresolved.add(label)
                logger.warning(f'unresolved identifier {label}: no NCBI taxid; using EsViritu species/order fallback if available')
    if dirty:
        save_identifier_cache(cache_path, cache)
    resolved = sum(1 for row in rows if row['pred_taxid'])
    logger.info(f'prediction identifiers resolved: {resolved}/{len(rows)} (cache hits: {cache_hits})')
    return rows


def resolve_taxonomy(taxids):
    """Map taxids to species-level match keys and orders via the NCBI lineage cache.

    Species-rank names are preferred. When a taxid has no species-rank lineage
    entry, the deepest available taxonomic name is used as the match key. This
    keeps ground-truth and predicted taxids symmetric for taxa absent from the
    species-rank cache. Failures without a usable fallback remain unmatchable.
    """
    info = {t: {'species': None, 'order': None} for t in taxids}
    ids = sorted({int(t) for t in taxids if str(t).isdigit()})
    if not ids:
        return info
    if PIPELINE_DIR not in sys.path:
        sys.path.insert(0, PIPELINE_DIR)
    try:
        from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
    except Exception as e:
        logger.error(f'taxonomy lookup unavailable ({e}); no taxid can be matched to species')
        return info
    try:
        wrapper = NCBITaxonomistWrapper(db=TAXA_DB if os.path.exists(TAXA_DB) else None)
        wrapper.resolve_lineages(ids)
    except Exception as e:
        logger.error(f'lineage resolution failed ({e}); no taxid can be matched to species')
        return info
    fallback_count = 0
    for t in taxids:
        if not str(t).isdigit():
            logger.error(f'taxid {t!r} is not numeric; cannot resolve species')
            continue
        tid = int(t)
        if tid not in wrapper.lineages:
            logger.error(f'taxid {t}: no lineage found; cannot resolve species')
            continue
        species = norm_species(wrapper.get_level(tid, NCBITaxonomistWrapper.TAX_SPECIES))
        info[t]['order'] = norm_rank(wrapper.get_level(tid, NCBITaxonomistWrapper.TAX_ORDER))
        if species is None:
            fallback = norm_species(wrapper.get_name(tid))
            if fallback and not is_generic_taxonomy_name(fallback):
                species = fallback
                fallback_count += 1
                logger.warning(
                    f'taxid {t}: no species-rank lineage entry; '
                    f'using deepest available name {fallback!r} as match key'
                )
            else:
                logger.error(f'taxid {t}: no usable species match key')
        info[t]['species'] = species
    logger.info(f'taxonomy species-key fallback used for {fallback_count} taxids')
    return info


def classify_dataset(ds, gt, tax_info, rows, min_reads):
    """Per-dataset records using the node_stats rule: TP only for best_match_is_best rows.

    - best row per GT species key = max(covered>0, read_count, mean_coverage); TP for GT taxa with a best row, else FN
    - other rows of a matched species key = CROSSHIT (not counted in precision/recall)
    - rows with no species match to any GT taxon = FP, retaining their predicted species and order
    - rows below min_reads are ignored
    """
    gt_by_species = defaultdict(list)
    for t in sorted(gt):
        sp = tax_info.get(t, {}).get('species')
        if sp:
            gt_by_species[sp].append(t)

    by_species = defaultdict(list)
    unmatched = []
    for r in rows:
        if r['read_count'] < min_reads:
            continue
        if r['species'] and r['species'] in gt_by_species:
            by_species[r['species']].append(r)
        else:
            unmatched.append(r)

    records = []
    best_of = {}
    for sp, rs in by_species.items():
        ranked = sorted(
            rs, key=lambda r: (r['covered_bases'] > 0, r['read_count'], r['mean_coverage']), reverse=True
        )
        best_of[sp] = ranked[0]
        first_gt = gt_by_species[sp][0]
        for r in ranked[1:]:
            records.append({
                'dataset': ds, 'taxid': ';'.join(gt_by_species[sp]), 'order': tax_info[first_gt]['order'],
                'species': sp, 'assembly': r['assembly'], 'read_count': r['read_count'], 'status': 'CROSSHIT',
            })

    for t in sorted(gt):
        sp = tax_info.get(t, {}).get('species')
        best = best_of.get(sp) if sp else None
        records.append({
            'dataset': ds, 'taxid': t, 'order': tax_info.get(t, {}).get('order'),
            'species': sp or '', 'assembly': best['assembly'] if best else '',
            'read_count': best['read_count'] if best else 0, 'status': 'TP' if best else 'FN',
        })

    for r in unmatched:
        records.append({
            'dataset': ds, 'taxid': '', 'order': r.get('order'), 'species': r.get('species') or '',
            'assembly': r['assembly'], 'read_count': r['read_count'], 'status': 'FP',
        })
    return records


def count_metrics(records):
    tp = sum(1 for r in records if r['status'] == 'TP')
    fp = sum(1 for r in records if r['status'] == 'FP')
    fn = sum(1 for r in records if r['status'] == 'FN')
    cross = sum(1 for r in records if r['status'] == 'CROSSHIT')
    prec = tp/(tp+fp) if (tp+fp)>0 else 0.0
    rec = tp/(tp+fn) if (tp+fn)>0 else 0.0
    f1 = 2*prec*rec/(prec+rec) if (prec+rec)>0 else 0.0
    return {'tp': tp, 'fp': fp, 'fn': fn, 'crosshit': cross, 'precision': prec, 'recall': rec, 'f1': f1}


def distribution_stats(values):
    """Mean/median/population-std/min/max for a list of numeric values."""
    vals = [float(v) for v in values]
    if not vals:
        return None
    return {
        'n_datasets': len(vals),
        'mean': statistics.mean(vals),
        'median': statistics.median(vals),
        'std': statistics.pstdev(vals),
        'min': min(vals),
        'max': max(vals),
    }


def setup_logging(output_dir, verbose):
    os.makedirs(output_dir, exist_ok=True)
    handlers = [
        logging.FileHandler(os.path.join(output_dir, LOG_NAME), mode='w'),
        logging.StreamHandler(sys.stderr),
    ]
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s %(levelname)s %(message)s',
        handlers=handlers,
        force=True,
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--esviritu-dir', required=True)
    p.add_argument('--sim-root', required=True)
    p.add_argument('--output-dir', required=True)
    p.add_argument('--datasets', help='comma separated')
    p.add_argument('--datasets-file')
    p.add_argument('--min-reads', type=float, default=1.0)
    p.add_argument('--plots', action='store_true')
    p.add_argument('--no-plots', action='store_true')
    p.add_argument('--verbose', action='store_true')
    args = p.parse_args()
    setup_logging(args.output_dir, args.verbose)
    plots = args.plots and not args.no_plots
    # discover datasets
    datasets = []
    if args.datasets:
        datasets = [x.strip() for x in args.datasets.split(',') if x.strip()]
    elif args.datasets_file and os.path.exists(args.datasets_file):
        with open(args.datasets_file) as fh:
            datasets = [x.strip() for x in fh if x.strip()]
    else:
        for d in os.listdir(args.sim_root):
            dp = os.path.join(args.sim_root, d)
            if os.path.isdir(dp) and os.path.exists(os.path.join(args.esviritu_dir, f'{d}.detected_virus.info.tsv')):
                datasets.append(d)
    datasets = sorted(set(datasets))
    logger.info(f'datasets: {len(datasets)}')

    gt_by_ds = {}
    rows_by_ds = {}
    all_taxids = set()
    for ds in datasets:
        gt, _ = load_gt_taxids(args.sim_root, ds)
        gt_by_ds[ds] = gt
        all_taxids.update(gt)
        rows = resolve_prediction_taxids(
            load_esviritu_rows(args.esviritu_dir, ds),
            cache_path=os.path.join(args.output_dir, IDENTIFIER_CACHE_NAME),
        )
        rows_by_ds[ds] = rows
        all_taxids.update(row['pred_taxid'] for row in rows if row['pred_taxid'])
    tax_info = resolve_taxonomy(sorted(all_taxids))
    for rows in rows_by_ds.values():
        for row in rows:
            info = tax_info.get(row['pred_taxid'], {}) if row.get('pred_taxid') else {}
            row['species'] = info.get('species') or row.get('esv_species')
            row['order'] = info.get('order') or row.get('esv_order')

    detailed = []
    per_ds = []
    for ds in datasets:
        gt = gt_by_ds[ds]
        rows = rows_by_ds[ds]
        records = classify_dataset(ds, gt, tax_info, rows, args.min_reads)
        detailed.extend(records)
        m = count_metrics(records)
        logger.info(f"{ds}: gt={len(gt)} tp={m['tp']} fp={m['fp']} fn={m['fn']} crosshit={m['crosshit']}")
        per_ds.append({'dataset': ds, 'n_gt': len(gt), 'n_pred_matched': m['tp'], **m})

    # overall (micro: pooled over all rows)
    m_all = count_metrics(detailed)
    overall = [{'metric': k, 'value': m_all[k]} for k in ('precision', 'recall', 'f1', 'tp', 'fp', 'fn', 'crosshit')]

    # macro: distribution of per-dataset precision/recall/f1
    aggregate = []
    for metric in ('precision', 'recall', 'f1'):
        stats = distribution_stats([row[metric] for row in per_ds])
        if stats:
            aggregate.append({'metric': metric, **stats})
    for row in aggregate:
        logger.info(
            f"macro {row['metric']}: n={row['n_datasets']} mean={row['mean']:.3f} "
            f"median={row['median']:.3f} std={row['std']:.3f} min={row['min']:.3f} max={row['max']:.3f}"
        )

    # order level (TP/FP/FN only, as before)
    order_stats = defaultdict(lambda: {'tp':0,'fp':0,'fn':0})
    for r in detailed:
        o = r.get('order') or 'unknown'
        if r['status']=='TP': order_stats[o]['tp']+=1
        elif r['status']=='FP': order_stats[o]['fp']+=1
        elif r['status']=='FN': order_stats[o]['fn']+=1
    order_out = []
    for o,s in order_stats.items():
        tp=s['tp']; fp=s['fp']; fn=s['fn']
        prec=tp/(tp+fp) if (tp+fp)>0 else 0.0
        rec=tp/(tp+fn) if (tp+fn)>0 else 0.0
        f1=2*prec*rec/(prec+rec) if (prec+rec)>0 else 0.0
        order_out.append({'order': o, 'tp': tp, 'fp': fp, 'fn': fn, 'precision': prec, 'recall': rec, 'f1': f1})
    order_out.sort(key=lambda x: x['order'])

    # write
    def wcsv(path, rows):
        if not rows: return
        with open(path,'w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter='\t')
            w.writeheader(); w.writerows(rows)
    wcsv(os.path.join(args.output_dir,'esviritu_sim_evaluation_detailed.tsv'), detailed)
    wcsv(os.path.join(args.output_dir,'esviritu_sim_evaluation_per_dataset.tsv'), per_ds)
    wcsv(os.path.join(args.output_dir,'esviritu_sim_evaluation_overall.tsv'), overall)
    wcsv(os.path.join(args.output_dir,'esviritu_sim_evaluation_aggregate.tsv'), aggregate)
    wcsv(os.path.join(args.output_dir,'esviritu_sim_evaluation_by_order.tsv'), order_out)
    meta = {
        'datasets': datasets,
        'min_reads': args.min_reads,
        'n_datasets': len(datasets),
        'taxonomy_species_fallback': 'deepest_available_name',
        'prediction_taxonomy_fallback': 'esviritu_species_and_order',
        'fp_orders_preserved': True,
    }
    with open(os.path.join(args.output_dir,'esviritu_sim_evaluation_metadata.json'),'w') as fh:
        json.dump(meta, fh, indent=2)
    if plots and HAS_MPL and per_ds:
        try:
            plt.figure(figsize=(max(6,len(per_ds)*0.6),4))
            ds=[r['dataset'] for r in per_ds]; rec=[r['recall'] for r in per_ds]; prec=[r['precision'] for r in per_ds]
            x=range(len(ds)); plt.bar([i-0.2 for i in x], rec, width=0.4, label='Recall')
            plt.bar([i+0.2 for i in x], prec, width=0.4, label='Precision')
            plt.xticks(x, ds, rotation=45, ha='right'); plt.ylim(0,1); plt.legend(); plt.tight_layout()
            plt.savefig(os.path.join(args.output_dir,'esviritu_sim_precision_recall_per_dataset.png'))
            plt.close()
        except Exception as e:
            logger.warning(f'per-dataset plot failed: {e}')
        try:
            if order_out:
                plt.figure(figsize=(max(6,len(order_out)*0.8),4))
                orders=[r['order'] for r in order_out]; rec=[r['recall'] for r in order_out]
                plt.bar(orders, rec); plt.xticks(rotation=45,ha='right'); plt.ylim(0,1); plt.title('Recall by Order'); plt.tight_layout()
                plt.savefig(os.path.join(args.output_dir,'esviritu_sim_recall_by_order.png'))
                plt.close()
        except Exception as e:
            logger.warning(f'order plot failed: {e}')
        try:
            if aggregate:
                labels = [r['metric'].capitalize() for r in aggregate]
                means = [r['mean'] for r in aggregate]
                errors = [r['std'] for r in aggregate]
                colors = ['#4c72b0', '#dd8452', '#55a868'][:len(labels)]
                plt.figure(figsize=(6,4))
                plt.bar(labels, means, yerr=errors, capsize=4, color=colors)
                plt.ylim(0,1)
                plt.ylabel('Score')
                plt.title(f'Aggregate metrics (mean +/- SD, n={len(per_ds)} datasets)')
                plt.tight_layout()
                plt.savefig(os.path.join(args.output_dir,'esviritu_sim_aggregate_metrics.png'))
                plt.close()
        except Exception as e:
            logger.warning(f'aggregate metrics plot failed: {e}')
        try:
            data = [
                [r['precision'] for r in per_ds],
                [r['recall'] for r in per_ds],
                [r['f1'] for r in per_ds],
            ]
            plt.figure(figsize=(6,4))
            plt.boxplot(data, showmeans=True)
            plt.xticks([1, 2, 3], ['Precision', 'Recall', 'F1'])
            plt.ylim(0,1)
            plt.ylabel('Score')
            plt.title(f'Per-dataset metric distribution (n={len(per_ds)})')
            plt.tight_layout()
            plt.savefig(os.path.join(args.output_dir,'esviritu_sim_metric_distribution.png'))
            plt.close()
        except Exception as e:
            logger.warning(f'metric distribution plot failed: {e}')
        try:
            plt.figure(figsize=(5,4))
            bottom = 0
            for value, label, color in (
                (m_all['tp'], 'TP', '#55a868'),
                (m_all['fp'], 'FP', '#c44e52'),
                (m_all['fn'], 'FN', '#8172b2'),
            ):
                plt.bar(['All datasets'], [value], bottom=bottom, label=label, color=color)
                bottom += value
            plt.ylabel('Count')
            plt.title('Aggregate match totals')
            plt.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(args.output_dir,'esviritu_sim_tp_fp_fn_totals.png'))
            plt.close()
        except Exception as e:
            logger.warning(f'match totals plot failed: {e}')
    logger.info(f'overall: {m_all}')
    print(f'Wrote evaluation to {args.output_dir}')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
