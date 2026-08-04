#!/usr/bin/env python3
"""
Full 50-antigen benchmark of ISIS against real IEDB-confirmed epitope positions.
"""
import sys
sys.path.insert(0, '/home/user/ISIS/src')

import json
import numpy as np
import matplotlib.pyplot as plt
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
# FIGURE 1: Threshold sweep + method comparison
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(15, 11))

thresholds = [1, 2, 3, 4, 5]
avg_sens = [np.mean([r['threshold_results'][t]['sens'] for r in all_results]) for t in thresholds]
avg_spec = [np.mean([r['threshold_results'][t]['spec'] for r in all_results]) for t in thresholds]
avg_ppv = [np.mean([r['threshold_results'][t]['ppv'] for r in all_results]) for t in thresholds]
avg_f1 = [np.mean([r['threshold_results'][t]['f1'] for r in all_results]) for t in thresholds]
avg_mcc = [np.mean([r['threshold_results'][t]['mcc'] for r in all_results]) for t in thresholds]

ax = axes[0, 0]
ax.plot(thresholds, avg_sens, 'g-o', label='Sensitivity (TPR)', linewidth=2, markersize=8)
ax.plot(thresholds, avg_spec, 'b-s', label='Specificity', linewidth=2, markersize=8)
ax.plot(thresholds, avg_ppv, 'orange', marker='^', label='PPV', linewidth=2, markersize=8)
ax.plot(thresholds, avg_f1, 'r-d', label='F1', linewidth=2, markersize=8)
ax.plot(thresholds, avg_mcc, 'purple', marker='*', label='MCC', linewidth=2, markersize=10)
ax.set_xlabel('B-cell Consensus Threshold (methods agreeing)')
ax.set_ylabel('Score')
ax.set_title(f'Threshold Sweep — N={len(all_results)} antigens (real IEDB ground truth)')
ax.legend()
ax.set_xticks(thresholds)
ax.grid(True, alpha=0.3)
ax.set_ylim(0, 1)

ax = axes[0, 1]
method_names = list(METHODS)
m_sens = [np.mean([p['sens'] for p in method_perf[m]]) for m in method_names]
m_spec = [np.mean([p['spec'] for p in method_perf[m]]) for m in method_names]
m_f1 = [np.mean([p['f1'] for p in method_perf[m]]) for m in method_names]

x = np.arange(len(method_names))
width = 0.25
ax.bar(x - width, m_sens, width, label='Sensitivity', color='green')
ax.bar(x, m_spec, width, label='Specificity', color='blue')
ax.bar(x + width, m_f1, width, label='F1', color='red')
ax.set_xticks(x)
ax.set_xticklabels(method_names, rotation=30, ha='right')
ax.set_ylabel('Score')
ax.set_title('Individual B-cell Method Performance (N=49)')
ax.legend()
ax.set_ylim(0, 1)

ax = axes[1, 0]
sens3 = [r['threshold_results'][3]['sens'] for r in all_results]
spec3 = [r['threshold_results'][3]['spec'] for r in all_results]
ax.scatter(spec3, sens3, alpha=0.6, s=60, c='steelblue', edgecolor='black')
ax.plot([0, 1], [1, 0], 'k--', alpha=0.3, label='Random')
ax.set_xlabel('Specificity')
ax.set_ylabel('Sensitivity')
ax.set_title('Per-Antigen Sens/Spec Distribution (Threshold=3)')
ax.legend()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
f1_dist = [r['threshold_results'][3]['f1'] for r in all_results]
ax.hist(f1_dist, bins=15, color='coral', edgecolor='black', alpha=0.8)
ax.axvline(np.mean(f1_dist), color='red', linestyle='--', linewidth=2,
           label=f'Mean={np.mean(f1_dist):.3f}')
ax.axvline(np.median(f1_dist), color='blue', linestyle='--', linewidth=2,
           label=f'Median={np.median(f1_dist):.3f}')
ax.set_xlabel('F1 Score')
ax.set_ylabel('Number of Antigens')
ax.set_title('F1 Score Distribution Across 49 Antigens')
ax.legend()

plt.tight_layout()
plt.savefig('benchmark/output/summary_threshold_sweep.png', dpi=150, bbox_inches='tight')
print("\nSaved: benchmark/output/summary_threshold_sweep.png")

# ============================================================================
# FIGURE 2: Strategy comparison
# ============================================================================
fig, ax = plt.subplots(figsize=(11, 7))

strategies = {
    'Consensus>=4\n(baseline)': [r['threshold_results'][4] for r in all_results],
    'Consensus>=3': [r['threshold_results'][3] for r in all_results],
    'Consensus>=2\n(lowered)': [r['threshold_results'][2] for r in all_results],
    'Consensus>=3\n+ MHC overlap': [r['combined_metrics'] for r in all_results],
}

metric_names = ['sens', 'spec', 'ppv', 'f1', 'mcc']
metric_labels = ['Sensitivity', 'Specificity', 'PPV', 'F1', 'MCC']
x = np.arange(len(metric_names))
width = 0.2

for i, (strategy_name, results) in enumerate(strategies.items()):
    vals = [np.mean([r[m] for r in results]) for m in metric_names]
    ax.bar(x + i*width - 1.5*width, vals, width, label=strategy_name)

ax.set_xticks(x)
ax.set_xticklabels(metric_labels)
ax.set_ylabel('Score')
ax.set_title('Strategy Comparison: Improving Sensitivity / Reducing False Negatives')
ax.legend(loc='upper right')
ax.set_ylim(0, 1)
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('benchmark/output/strategy_comparison.png', dpi=150, bbox_inches='tight')
print("Saved: benchmark/output/strategy_comparison.png")

# ============================================================================
# FIGURE 3: Full results table
# ============================================================================
fig, ax = plt.subplots(figsize=(14, 13))
ax.axis('off')

sorted_results = sorted(all_results, key=lambda r: -r['threshold_results'][3]['f1'])

table_data = []
for r in sorted_results:
    m = r['threshold_results'][3]
    table_data.append([
        r['accession'], r['organism'][:28], str(r['seq_len']),
        f"{m['sens']:.2f}", f"{m['spec']:.2f}", f"{m['ppv']:.2f}", f"{m['f1']:.2f}"
    ])

col_labels = ['Accession', 'Organism', 'Len', 'Sens', 'Spec', 'PPV', 'F1']
table = ax.table(cellText=table_data, colLabels=col_labels, loc='center',
                  cellLoc='center', colWidths=[0.15, 0.32, 0.08, 0.1, 0.1, 0.1, 0.1])
table.auto_set_font_size(False)
table.set_fontsize(8)
table.scale(1, 1.15)

for j in range(len(col_labels)):
    table[(0, j)].set_facecolor('#4472C4')
    table[(0, j)].set_text_props(color='white', fontweight='bold')

for i, r in enumerate(sorted_results):
    f1 = r['threshold_results'][3]['f1']
    color = '#C6EFCE' if f1 > 0.5 else '#FFEB9C' if f1 > 0.3 else '#FFC7CE'
    for j in range(len(col_labels)):
        table[(i+1, j)].set_facecolor(color)

ax.set_title(f'Full Benchmark Results — All {len(all_results)} Antigens (sorted by F1, Threshold=3)',
             fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('benchmark/output/full_results_table.png', dpi=150, bbox_inches='tight')
print("Saved: benchmark/output/full_results_table.png")

print("\n" + "="*80)
print("BENCHMARK COMPLETE")
print("="*80)
