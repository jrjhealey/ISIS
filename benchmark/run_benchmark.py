#!/usr/bin/env python3
"""
Full 50-antigen benchmark of ISIS against real IEDB-confirmed epitope positions.
"""
import sys
sys.path.insert(0, '/home/user/ISIS/src')

import json
import numpy as np
from collections import defaultdict
import os

os.makedirs('benchmark/output', exist_ok=True)

with open('benchmark/dataset/final_benchmark_set.json') as f:
    ANTIGENS = json.load(f)

print(f"Loaded {len(ANTIGENS)} benchmark antigens")

from isis import predict, available_methods
from isis.models.mhc_predictor import MHCPredictor

METHODS = available_methods()
predictor = MHCPredictor()
MHC_ALLELES = ['HLA-A*02:01', 'HLA-B*07:02', 'HLA-DRB1*01:01']


def calculate_metrics(pred_positions, gt_positions, seq_len):
    all_pos = set(range(seq_len))
    tp = len(pred_positions & gt_positions)
    fp = len(pred_positions - gt_positions)
    fn = len(gt_positions - pred_positions)
    tn = len(all_pos - pred_positions - gt_positions)

    sens = tp / (tp + fn) if (tp + fn) > 0 else 0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    f1 = 2 * ppv * sens / (ppv + sens) if (ppv + sens) > 0 else 0
    denom = np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    mcc = (tp*tn - fp*fn) / denom if denom > 0 else 0
    acc = (tp + tn) / (tp + fp + fn + tn)

    return {'sens': sens, 'spec': spec, 'ppv': ppv, 'npv': npv, 'f1': f1,
            'mcc': mcc, 'acc': acc, 'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn}


all_results = []
method_perf = defaultdict(list)

print("\nRunning predictions on all antigens...")
for idx, antigen in enumerate(ANTIGENS):
    seq = antigen['sequence']
    gt_positions = set(antigen['epitope_positions'])
    seq_len = antigen['seq_len']

    print(f"  [{idx+1}/{len(ANTIGENS)}] {antigen['accession']} ({antigen['organism'][:30]}, {seq_len}aa)...")

    # B-cell linear methods
    bcell_results = {}
    for method in METHODS:
        try:
            result = predict(seq, method=method)
            bcell_results[method] = [(ep.start, ep.end) for ep in result.epitopes]
        except Exception:
            bcell_results[method] = []

    # Per-method metrics
    for method, epitopes in bcell_results.items():
        pred_pos = set()
        for start, end in epitopes:
            for i in range(start-1, min(end, seq_len)):
                pred_pos.add(i)
        m = calculate_metrics(pred_pos, gt_positions, seq_len)
        method_perf[method].append(m)

    # Consensus votes
    consensus = np.zeros(seq_len)
    for method, epitopes in bcell_results.items():
        for start, end in epitopes:
            for i in range(start-1, min(end, seq_len)):
                consensus[i] += 1

    # MHC-I/II binder positions (as additional evidence layer)
    mhc_strong_positions = set()
    for allele in MHC_ALLELES:
        try:
            if 'DRB' in allele:
                pred = predictor.predict_mhc2(seq, allele=allele)
            else:
                pred = predictor.predict_mhc1(seq, allele=allele)
            for pep in pred['peptides']:
                if pep['binding'] in ('strong', 'weak'):
                    for i in range(pep['start']-1, min(pep['end'], seq_len)):
                        mhc_strong_positions.add(i)
        except Exception:
            pass

    threshold_results = {}
    for thresh in [1, 2, 3, 4, 5]:
        pred_pos = set(i for i, v in enumerate(consensus) if v >= thresh)
        m = calculate_metrics(pred_pos, gt_positions, seq_len)
        threshold_results[thresh] = m

    # Combined: consensus>=3 OR mhc binder overlap (tests improvement hypothesis)
    combined_pos = set(i for i, v in enumerate(consensus) if v >= 3) | mhc_strong_positions
    combined_metrics = calculate_metrics(combined_pos, gt_positions, seq_len)

    all_results.append({
        'accession': antigen['accession'],
        'organism': antigen['organism'],
        'seq_len': seq_len,
        'gt_len': len(gt_positions),
        'coverage': antigen['coverage'],
        'threshold_results': threshold_results,
        'combined_metrics': combined_metrics,
    })

print("\nAll predictions complete. Computing summary statistics...\n")

with open('benchmark/output/raw_results.json', 'w') as f:
    json.dump(all_results, f, default=lambda x: float(x) if isinstance(x, np.floating) else x)


def summarize(key_extractor, results):
    vals = [key_extractor(r) for r in results]
    return {'mean': np.mean(vals), 'std': np.std(vals), 'median': np.median(vals),
            'min': np.min(vals), 'max': np.max(vals)}

print("="*80)
print(f"BENCHMARK SUMMARY: {len(all_results)} antigens, "
      f"{sum(r['seq_len'] for r in all_results):,} total residues, "
      f"{sum(r['gt_len'] for r in all_results):,} ground-truth epitope positions")
print("="*80)

for thresh in [1, 2, 3, 4, 5]:
    sens = summarize(lambda r: r['threshold_results'][thresh]['sens'], all_results)
    spec = summarize(lambda r: r['threshold_results'][thresh]['spec'], all_results)
    ppv = summarize(lambda r: r['threshold_results'][thresh]['ppv'], all_results)
    f1 = summarize(lambda r: r['threshold_results'][thresh]['f1'], all_results)
    mcc = summarize(lambda r: r['threshold_results'][thresh]['mcc'], all_results)
    print(f"\nThreshold={thresh}:")
    print(f"  Sensitivity: {sens['mean']:.3f} +/- {sens['std']:.3f}")
    print(f"  Specificity: {spec['mean']:.3f} +/- {spec['std']:.3f}")
    print(f"  PPV:         {ppv['mean']:.3f} +/- {ppv['std']:.3f}")
    print(f"  F1:          {f1['mean']:.3f} +/- {f1['std']:.3f}")
    print(f"  MCC:         {mcc['mean']:.3f} +/- {mcc['std']:.3f}")

print("\n" + "="*80)
print("COMBINED STRATEGY (consensus>=3 OR MHC binder overlap):")
sens = summarize(lambda r: r['combined_metrics']['sens'], all_results)
spec = summarize(lambda r: r['combined_metrics']['spec'], all_results)
f1 = summarize(lambda r: r['combined_metrics']['f1'], all_results)
mcc = summarize(lambda r: r['combined_metrics']['mcc'], all_results)
print(f"  Sensitivity: {sens['mean']:.3f} +/- {sens['std']:.3f}")
print(f"  Specificity: {spec['mean']:.3f} +/- {spec['std']:.3f}")
print(f"  F1:          {f1['mean']:.3f} +/- {f1['std']:.3f}")
print(f"  MCC:         {mcc['mean']:.3f} +/- {mcc['std']:.3f}")

print("\n" + "="*80)
print("PER-METHOD PERFORMANCE:")
print(f"{'Method':<22} {'Sens':>8} {'Spec':>8} {'PPV':>8} {'F1':>8} {'MCC':>8}")
print("-"*70)
for method in METHODS:
    perf = method_perf[method]
    sens = np.mean([p['sens'] for p in perf])
    spec = np.mean([p['spec'] for p in perf])
    ppv = np.mean([p['ppv'] for p in perf])
    f1 = np.mean([p['f1'] for p in perf])
    mcc = np.mean([p['mcc'] for p in perf])
    print(f"{method:<22} {sens:>8.3f} {spec:>8.3f} {ppv:>8.3f} {f1:>8.3f} {mcc:>8.3f}")

# ============================================================================
# FIGURES - all styling comes from isis.plotting so these match the figures
# produced for individual sequences (same palette, axes, fonts, spacing).
# ============================================================================
from isis import plotting as P

SUB = (f"{len(all_results)} antigens · "
       f"{sum(r['seq_len'] for r in all_results):,} residues · "
       f"IEDB-confirmed epitope positions")

thresholds = [1, 2, 3, 4, 5]
sweep = {
    key: [np.mean([r['threshold_results'][t][key] for r in all_results])
          for t in thresholds]
    for key in ('sens', 'spec', 'ppv', 'f1', 'mcc')
}
P.save_figure(
    P.plot_threshold_sweep(thresholds, sweep, subtitle=SUB),
    'benchmark/output/threshold_sweep.png')
print("Saved: benchmark/output/threshold_sweep.png")

per_method = {
    m: {k: float(np.mean([p[k] for p in perf])) for k in ('sens', 'spec', 'ppv', 'f1')}
    for m, perf in method_perf.items()
}
P.save_figure(
    P.plot_metric_comparison(per_method, subtitle=SUB),
    'benchmark/output/method_comparison.png')
print("Saved: benchmark/output/method_comparison.png")

# MCC is the headline: it is zero for random guessing, so it is the metric that
# shows several of these methods carry no real signal.
mcc_by_method = {
    m: float(np.mean([p['mcc'] for p in perf])) for m, perf in method_perf.items()
}
P.save_figure(
    P.plot_mcc_ranking(mcc_by_method,
                       subtitle=SUB + " · 0 = indistinguishable from random"),
    'benchmark/output/mcc_ranking.png')
print("Saved: benchmark/output/mcc_ranking.png")

P.save_figure(
    P.plot_score_distribution(
        [r['threshold_results'][3]['f1'] for r in all_results],
        xlabel="F1 score", title="F1 distribution across antigens",
        subtitle=f"Consensus threshold = 3 · n = {len(all_results)} antigens"),
    'benchmark/output/f1_distribution.png')
print("Saved: benchmark/output/f1_distribution.png")

print("\n" + "="*80)
print("BENCHMARK COMPLETE")
print("="*80)
