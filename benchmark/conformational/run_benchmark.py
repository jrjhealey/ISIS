#!/usr/bin/env python3
"""
Benchmark the structure-based epitope predictors against observed
antibody-contact epitopes.

Reports AUC as the headline metric because that is what the conformational
epitope literature reports, which makes the numbers directly comparable:
SEMA 0.76, DiscoTope-2.0 and ElliPro below that, random 0.50.

Three controls are included deliberately, because the interesting question is
not "does DiscoTope beat random" but "does DiscoTope beat simply asking which
residues are on the surface":

  sasa_only        - relative solvent accessibility as the score
  protrusion_only  - protrusion index as the score
  contacts_only    - inverted contact number (fewer contacts = more exposed)

If a composite method cannot beat its own cheapest input feature, its extra
machinery is not earning anything. The sequence-only B-cell scales are run over
the same antigens too, so structure-based and sequence-based methods are
compared on identical ground truth.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

import numpy as np
from sklearn.metrics import roc_auc_score, matthews_corrcoef

from isis import predict, available_methods
from isis.methods.bcell_conformational import (
    DiscoTopePredictor, ElliProPredictor, SEPPAPredictor,
    calculate_protrusion_index,
)
from isis.methods.structural import calculate_contact_number

DATA = "benchmark/conformational/dataset/complexes.json"
OUTDIR = "benchmark/conformational/output"


def metrics(y_true, scores, y_pred=None):
    """AUC over the continuous score, plus MCC/F1 at a hard threshold."""
    out = {}
    if len(set(y_true)) > 1 and np.any(np.isfinite(scores)):
        s = np.nan_to_num(np.asarray(scores, dtype=float),
                          nan=float(np.nanmin(scores)) if np.any(np.isfinite(scores)) else 0.0)
        out["auc"] = float(roc_auc_score(y_true, s))
    else:
        out["auc"] = float("nan")

    if y_pred is not None:
        tp = int(((y_pred == 1) & (y_true == 1)).sum())
        fp = int(((y_pred == 1) & (y_true == 0)).sum())
        fn = int(((y_pred == 0) & (y_true == 1)).sum())
        tn = int(((y_pred == 0) & (y_true == 0)).sum())
        sens = tp / (tp + fn) if tp + fn else 0.0
        spec = tn / (tn + fp) if tn + fp else 0.0
        ppv = tp / (tp + fp) if tp + fp else 0.0
        f1 = 2 * ppv * sens / (ppv + sens) if ppv + sens else 0.0
        mcc = (matthews_corrcoef(y_true, y_pred)
               if len(set(y_pred)) > 1 and len(set(y_true)) > 1 else 0.0)
        out.update(sens=sens, spec=spec, ppv=ppv, f1=f1, mcc=float(mcc))
    return out


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    with open(DATA) as f:
        complexes = json.load(f)
    print(f"Loaded {len(complexes)} antibody-antigen complexes\n")

    per_method = {}

    def record(name, m):
        per_method.setdefault(name, []).append(m)

    linear_methods = list(available_methods())

    for i, cx in enumerate(complexes, 1):
        seq = cx["sequence"]
        n = len(seq)
        y = np.zeros(n, dtype=int)
        y[cx["epitope_idx"]] = 1

        sasa = np.array(cx["sasa"], dtype=float)
        coords = np.array(cx["ca_coords"], dtype=float)
        contacts = calculate_contact_number(coords)
        protrusion = calculate_protrusion_index(coords)

        print(f"[{i}/{len(complexes)}] {cx['pdb_id']}/{cx['antigen_chain']} "
              f"{n} aa, {len(cx['epitope_idx'])} epitope res")

        # ---- structure-based predictors ----
        for name, predictor in (("discotope", DiscoTopePredictor()),
                                ("ellipro", ElliProPredictor()),
                                ("seppa", SEPPAPredictor())):
            try:
                res = predictor.predict_with_structure(seq, sasa, contacts, coords)
                scores = np.asarray(res["scores"], dtype=float)
                thr = res.get("threshold")
                pred = (np.nan_to_num(scores, nan=-np.inf) >= thr).astype(int) \
                    if thr is not None else None
                record(name, metrics(y, scores, pred))
            except Exception as e:
                print(f"    {name}: FAILED {e}")

        # ---- single-feature controls ----
        record("sasa_only", metrics(y, sasa,
                                    (sasa >= np.nanmedian(sasa)).astype(int)))
        record("protrusion_only", metrics(y, protrusion,
                                         (protrusion >= np.nanmedian(protrusion)).astype(int)))
        record("contacts_only", metrics(y, -contacts.astype(float),
                                        (contacts <= np.median(contacts)).astype(int)))

        # ---- sequence-only B-cell scales, same antigens, same ground truth ----
        for method in linear_methods:
            try:
                p = predict(seq, method=method)
                scores = np.array([p.score_at(k + 1) if p.score_at(k + 1) is not None
                                   else np.nan for k in range(n)], dtype=float)
                pred = np.zeros(n, dtype=int)
                for ep in p.epitopes:
                    pred[ep.start - 1:ep.end] = 1
                record(f"seq:{method}", metrics(y, scores, pred))
            except Exception as e:
                print(f"    {method}: FAILED {e}")

    # ------------------------------------------------------------------
    summary = {}
    for name, runs in per_method.items():
        aucs = [r["auc"] for r in runs if not np.isnan(r["auc"])]
        summary[name] = {
            "n": len(runs),
            "auc": float(np.mean(aucs)) if aucs else float("nan"),
            "auc_sd": float(np.std(aucs)) if aucs else float("nan"),
            "auc_median": float(np.median(aucs)) if aucs else float("nan"),
            "frac_above_random": (float(np.mean([a > 0.5 for a in aucs]))
                                  if aucs else float("nan")),
            "mcc": float(np.mean([r.get("mcc", 0.0) for r in runs])),
            "f1": float(np.mean([r.get("f1", 0.0) for r in runs])),
            "sens": float(np.mean([r.get("sens", 0.0) for r in runs])),
            "spec": float(np.mean([r.get("spec", 0.0) for r in runs])),
        }

    with open(os.path.join(OUTDIR, "results.json"), "w") as f:
        json.dump({"summary": summary, "per_method": per_method}, f)

    order = sorted(summary, key=lambda k: -(summary[k]["auc"]
                                            if not np.isnan(summary[k]["auc"]) else -1))
    print("\n" + "=" * 92)
    print(f"CONFORMATIONAL EPITOPE BENCHMARK - {len(complexes)} antibody-antigen "
          f"complexes, {sum(len(c['sequence']) for c in complexes):,} residues")
    print("Ground truth: antigen residues within 4 A of antibody (Ponomarenko & Bourne)")
    print("=" * 92)
    print(f"{'method':<26}{'AUC':>8}{'±sd':>7}{'med':>7}{'>0.5':>7}"
          f"{'MCC':>8}{'F1':>7}{'Sens':>7}{'Spec':>7}")
    print("-" * 92)
    for name in order:
        s = summary[name]
        print(f"{name:<26}{s['auc']:>8.3f}{s['auc_sd']:>7.3f}{s['auc_median']:>7.3f}"
              f"{s['frac_above_random']:>7.0%}{s['mcc']:>8.3f}{s['f1']:>7.3f}"
              f"{s['sens']:>7.3f}{s['spec']:>7.3f}")
    print("-" * 92)
    print("Reference points from the literature: SEMA 0.76, DiscoTope-2.0 < that, "
          "random = 0.50")

    # ---- figures, in the shared ISIS style ----
    try:
        from isis import plotting as P
        sub = (f"{len(complexes)} antibody-antigen complexes · "
               f"{sum(len(c['sequence']) for c in complexes):,} residues · "
               f"epitope = antigen residue within 4 Å of antibody")

        auc_by_method = {k: summary[k]["auc"] for k in summary
                         if not np.isnan(summary[k]["auc"])}
        P.save_figure(
            P.plot_metric_ranking(auc_by_method, reference=0.5,
                                  title="Conformational epitope prediction — AUC",
                                  subtitle=sub, xlabel="ROC AUC"),
            os.path.join(OUTDIR, "auc_ranking.png"))
        print(f"\nSaved: {OUTDIR}/auc_ranking.png")

        best = [k for k in order if not np.isnan(summary[k]["auc"])][:6]
        P.save_figure(
            P.plot_metric_comparison(
                {k: {m: summary[k][m] for m in ("sens", "spec", "ppv", "f1")
                     if m in summary[k]} for k in best},
                metrics=("sens", "spec", "f1"),
                title="Conformational predictors vs controls",
                subtitle=sub),
            os.path.join(OUTDIR, "metric_comparison.png"))
        print(f"Saved: {OUTDIR}/metric_comparison.png")
    except Exception as e:
        print(f"\n(figures skipped: {e})")


if __name__ == "__main__":
    main()
