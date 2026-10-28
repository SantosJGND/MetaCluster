#!/usr/bin/env python3
import argparse
import csv
import json
import os
import re
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

def norm_taxid(v):
    if v is None: return None
    s = str(v).strip()
    if not s or s.lower() in {'nan','none','null','na'}: return None
    m = re.findall(r'\d+', s)
    return m[0] if m else s

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

def load_esviritu_hits(esviritu_dir, dataset, min_reads):
    path = os.path.join(esviritu_dir, f'{dataset}.detected_virus.info.tsv')
    if not os.path.exists(path):
        return [], []
    hits = []
    for row in read_table(path):
        name = str(row.get('species') or row.get('Name') or row.get('name') or '').strip()
        taxid = norm_taxid(row.get('taxid') or row.get('TaxID'))
        reads = 0.0
        try:
            reads = float(row.get('read_count') or row.get('reads') or 0)
        except Exception:
            reads = 0.0
        hits.append({'taxid': taxid, 'name': name, 'reads': reads, 'row': row})
    detected = [h for h in hits if h['reads'] >= min_reads]
    return hits, detected

def get_taxonomy_order(taxids):
    try:
        from metagenomics_utils.ncbi_tools import NCBITaxonomistWrapper
        wrapper = NCBITaxonomistWrapper(db=TAXA_DB if os.path.exists(TAXA_DB) else None)
        wrapper.resolve_lineages([int(t) for t in taxids if str(t).isdigit()])
        out = {}
        for t in taxids:
            try:
                tid = int(t)
                lin = wrapper.lineages.get(tid, {})
                order = None
                for k,v in lin.items():
                    if isinstance(v, dict) and v.get('level')==3:  # rough
                        pass
                # try common
                if 'order' in lin and isinstance(lin['order'], dict):
                    order = lin['order'].get('name')
                out[t] = order
            except Exception:
                out[t] = None
        return out
    except Exception as e:
        return {t: None for t in taxids}

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
    os.makedirs(args.output_dir, exist_ok=True)
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
    detailed = []
    per_ds = []
    all_gt = set()
    all_pred = set()
    # get orders
    all_taxids = []
    for ds in datasets:
        gt,_ = load_gt_taxids(args.sim_root, ds)
        all_taxids.extend(list(gt))
    order_map = get_taxonomy_order(sorted(set(all_taxids)))
    for ds in datasets:
        gt, gt_info = load_gt_taxids(args.sim_root, ds)
        all_gt.update(gt)
        hits, detected = load_esviritu_hits(args.esviritu_dir, ds, args.min_reads)
        # try to match detected to gt by taxid if present
        pred_taxids = set()
        for h in detected:
            if h['taxid'] and h['taxid'] in gt:
                pred_taxids.add(h['taxid'])
            elif h['taxid']:
                # if taxid present but not in gt, still count as pred taxid for off-panel? but compare to gt
                pred_taxids.add(h['taxid'])
        all_pred.update([p for p in pred_taxids if p in gt])
        tp = len(gt & pred_taxids)
        fp = len(pred_taxids - gt)
        fn = len(gt - pred_taxids)
        prec = tp/(tp+fp) if (tp+fp)>0 else 0.0
        rec = tp/(tp+fn) if (tp+fn)>0 else 0.0
        f1 = 2*prec*rec/(prec+rec) if (prec+rec)>0 else 0.0
        per_ds.append({'dataset': ds, 'n_gt': len(gt), 'n_pred_matched': len(pred_taxids & gt), 'tp': tp, 'fp': fp, 'fn': fn, 'precision': prec, 'recall': rec, 'f1': f1})
        for tid in gt:
            matched = tid in pred_taxids
            status = 'TP' if matched else 'FN'
            detailed.append({'dataset': ds, 'taxid': tid, 'order': order_map.get(tid), 'status': status})
        for pt in pred_taxids - gt:
            detailed.append({'dataset': ds, 'taxid': pt, 'order': order_map.get(pt), 'status': 'FP'})
    # overall
    tp_all = sum(1 for r in detailed if r['status']=='TP')
    fp_all = sum(1 for r in detailed if r['status']=='FP')
    fn_all = sum(1 for r in detailed if r['status']=='FN')
    prec_all = tp_all/(tp_all+fp_all) if (tp_all+fp_all)>0 else 0.0
    rec_all = tp_all/(tp_all+fn_all) if (tp_all+fn_all)>0 else 0.0
    f1_all = 2*prec_all*rec_all/(prec_all+rec_all) if (prec_all+rec_all)>0 else 0.0
    overall = [{'metric': 'precision', 'value': prec_all}, {'metric': 'recall', 'value': rec_all}, {'metric': 'f1', 'value': f1_all}, {'metric': 'tp', 'value': tp_all}, {'metric': 'fp', 'value': fp_all}, {'metric': 'fn', 'value': fn_all}]
    # order level
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
    wcsv(os.path.join(args.output_dir,'esviritu_sim_evaluation_by_order.tsv'), order_out)
    meta = {'datasets': datasets, 'min_reads': args.min_reads, 'n_datasets': len(datasets)}
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
        except Exception:
            pass
        try:
            if order_out:
                plt.figure(figsize=(max(6,len(order_out)*0.8),4))
                orders=[r['order'] for r in order_out]; rec=[r['recall'] for r in order_out]
                plt.bar(orders, rec); plt.xticks(rotation=45,ha='right'); plt.ylim(0,1); plt.title('Recall by Order'); plt.tight_layout()
                plt.savefig(os.path.join(args.output_dir,'esviritu_sim_recall_by_order.png'))
                plt.close()
        except Exception:
            pass
    print(f'Wrote evaluation to {args.output_dir}')
    return 0

if __name__=='__main__':
    raise SystemExit(main())
