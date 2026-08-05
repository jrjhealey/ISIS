#!/usr/bin/env python3
"""
Fit and honestly evaluate a combined structural-epitope predictor.

Answers "which methods would you combine" with a measurement rather than an
opinion, by fitting weights over the available structural features and scoring
with leave-one-antigen-out cross-validation. Grouping by antigen matters: a
random split over residues would put residues of the same protein in both
train and test, and neighbouring residues of one protein are highly correlated,
so a residue-level split would report a number that does not survive contact
with a new antigen.

Feature sets are nested so the marginal value of each addition is visible,
which is the only way to tell whether a composite is earning its complexity.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon

from isis.methods.bcell_conformational import (
    DiscoTopePredictor, ElliProPredictor, SEPPAPredictor,
    calculate_protrusion_index,
)
from isis.methods.structural import calculate_contact_number

DATA = "benchmark/conformational/dataset/complexes.json"
OUTDIR = "benchmark/conformational/output"

FEATURES = ["sasa_rel", "protrusion", "contacts", "discotope", "ellipro", "seppa"]

SETS = {
    "sasa only":                    ["sasa_rel"],
    "sasa + protrusion":            ["sasa_rel", "protrusion"],
    "sasa + protrusion + contacts": ["sasa_rel", "protrusion", "contacts"],
    "3 published methods":          ["discotope", "ellipro", "seppa"],
    "everything":                   FEATURES,
}


def build():
    with open(DATA) as f:
        complexes = json.load(f)

    X_all, y_all, groups = [], [], []
    for gi, cx in enumerate(complexes):
        seq = cx["sequence"]
        n = len(seq)
        sasa = np.array(cx["sasa"], dtype=float)
        coords = np.array(cx["ca_coords"], dtype=float)
        contacts = calculate_contact_number(coords).astype(float)
        protrusion = calculate_protrusion_index(coords)

        cols = {
            # Normalised per protein: absolute SASA is not comparable across
            # antigens of different size, and the model must generalise to a
            # protein whose scale it has never seen.
            "sasa_rel": _z(sasa),
            "protrusion": _z(protrusion),
            "contacts": _z(contacts),
        }
        for name, pred in (("discotope", DiscoTopePredictor()),
                           ("ellipro", ElliProPredictor()),
                           ("seppa", SEPPAPredictor())):
            r = pred.predict_with_structure(seq, sasa, contacts, coords)
            cols[name] = _z(np.asarray(r["scores"], dtype=float))

        y = np.zeros(n, dtype=int)
        y[cx["epitope_idx"]] = 1

        X_all.append(np.column_stack([cols[f] for f in FEATURES]))
        y_all.append(y)
        groups.append(np.full(n, gi))

    return (np.vstack(X_all), np.concatenate(y_all), np.concatenate(groups),
            complexes)


def _z(a):
    a = np.asarray(a, dtype=float)
    finite = a[np.isfinite(a)]
    if finite.size == 0:
        return np.zeros_like(a)
    mu, sd = finite.mean(), finite.std()
    out = (a - mu) / sd if sd > 1e-12 else np.zeros_like(a)
    return np.nan_to_num(out, nan=0.0)


def loao_auc(X, y, groups, cols, model="lr"):
    """Leave-one-antigen-out: per-antigen AUC from a model that never saw it."""
    idx = [FEATURES.index(c) for c in cols]
    aucs = []
    for g in np.unique(groups):
        te = groups == g
        tr = ~te
        if len(set(y[te])) < 2:
            continue
        sc = StandardScaler().fit(X[tr][:, idx])
        if model == "lr":
            clf = LogisticRegression(max_iter=2000, class_weight="balanced")
        else:
            clf = RandomForestClassifier(n_estimators=200, max_depth=6,
                                         min_samples_leaf=20,
                                         class_weight="balanced",
                                         n_jobs=-1, random_state=0)
        clf.fit(sc.transform(X[tr][:, idx]), y[tr])
        p = clf.predict_proba(sc.transform(X[te][:, idx]))[:, 1]
        aucs.append(roc_auc_score(y[te], p))
    return np.array(aucs)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    X, y, groups, complexes = build()
    print(f"{len(complexes)} antigens, {len(y):,} residues, "
          f"{y.sum():,} epitope ({y.mean():.1%})\n")

    results = {}
    print(f"{'feature set':<32}{'LOAO AUC':>10}{'±sd':>7}{'>0.5':>7}")
    print("-" * 56)
    for label, cols in SETS.items():
        a = loao_auc(X, y, groups, cols)
        results[label] = a
        print(f"{label:<32}{a.mean():>10.3f}{a.std():>7.3f}{(a > 0.5).mean():>7.0%}")

    rf = loao_auc(X, y, groups, FEATURES, model="rf")
    results["everything (random forest)"] = rf
    print(f"{'everything (random forest)':<32}{rf.mean():>10.3f}{rf.std():>7.3f}"
          f"{(rf > 0.5).mean():>7.0%}")

    print("\n=== does each addition earn its place? (paired Wilcoxon) ===")
    base = results["sasa only"]
    for label in ("sasa + protrusion", "sasa + protrusion + contacts",
                  "3 published methods", "everything",
                  "everything (random forest)"):
        a = results[label]
        _s, p = wilcoxon(a, base)
        tag = ("better" if p < 0.05 and a.mean() > base.mean()
               else "worse" if p < 0.05 else "no different")
        print(f"  {label:<32} {a.mean():.3f} vs {base.mean():.3f} "
              f"delta {a.mean()-base.mean():+.3f}  p={p:.4f}  -> {tag} than SASA alone")

    # Weights from a model fitted on everything, for interpretation only
    sc = StandardScaler().fit(X)
    clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit(sc.transform(X), y)
    print("\n=== fitted weights (all antigens; interpretation only) ===")
    for f, w in sorted(zip(FEATURES, clf.coef_[0]), key=lambda t: -abs(t[1])):
        print(f"  {f:<14}{w:+.3f}")

    with open(os.path.join(OUTDIR, "consensus_results.json"), "w") as f:
        json.dump({k: v.tolist() for k, v in results.items()}, f)

    try:
        from isis import plotting as P
        P.save_figure(
            P.plot_metric_ranking(
                {k: float(v.mean()) for k, v in results.items()},
                reference=0.5,
                title="Combined structural predictor — leave-one-antigen-out AUC",
                subtitle=(f"{len(complexes)} antibody-antigen complexes · "
                          f"grouped CV, no antigen in its own training set"),
                xlabel="ROC AUC"),
            os.path.join(OUTDIR, "consensus_auc.png"))
        print(f"\nSaved: {OUTDIR}/consensus_auc.png")
    except Exception as e:
        print(f"(figure skipped: {e})")


if __name__ == "__main__":
    main()
