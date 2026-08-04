#!/usr/bin/env python3
"""
Train MHC binding prediction models on real IEDB data.

This script parses the actual IEDB MHC ligand database and trains
Random Forest models for MHC-I and MHC-II binding prediction.

Usage:
    python train_mhc_iedb.py --data /path/to/mhc_ligand_full.csv --output models/
"""

import argparse
import csv
import os
import pickle
import sys
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np

try:
    from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import mean_squared_error, r2_score, roc_auc_score
    import joblib
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


# Column indices in IEDB mhc_ligand_full.csv (0-indexed)
COL_PEPTIDE = 11          # Epitope Name (peptide sequence)
COL_MHC_ALLELE = 107      # MHC Allele Name
COL_QUANT_MEASURE = 97    # Quantitative measurement (IC50)
COL_QUAL_MEASURE = 95     # Qualitative Measurement
COL_MEASURE_INEQ = 96     # Measurement inequality (<, >, =)
COL_ASSAY_TYPE = 91       # Method/Assay type


# Amino acid encoding
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# BLOSUM62 for encoding
BLOSUM62 = {
    'A': [4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0],
    'R': [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3],
    'N': [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3],
    'D': [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3],
    'C': [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    'Q': [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2],
    'E': [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2],
    'G': [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3],
    'H': [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3],
    'I': [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3],
    'L': [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1],
    'K': [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2],
    'M': [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1],
    'F': [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1],
    'P': [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2],
    'S': [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2],
    'T': [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0],
    'W': [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3],
    'Y': [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1],
    'V': [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4],
}

# Qualitative to numeric mapping
QUAL_TO_SCORE = {
    'positive-high': 50,      # Strong binder
    'positive-intermediate': 200,
    'positive-low': 1000,
    'positive': 500,
    'negative': 20000,
}


def encode_peptide(peptide: str, max_length: int = 15) -> np.ndarray:
    """Encode peptide using one-hot + BLOSUM62 + length features."""
    # Pad or truncate to max_length
    peptide = peptide[:max_length].ljust(max_length, 'X')

    features = []

    # One-hot encoding
    for aa in peptide:
        onehot = [0] * 20
        if aa in AA_TO_IDX:
            onehot[AA_TO_IDX[aa]] = 1
        features.extend(onehot)

    # BLOSUM62 encoding
    for aa in peptide:
        if aa in BLOSUM62:
            features.extend(BLOSUM62[aa])
        else:
            features.extend([0] * 20)

    # Length feature
    actual_len = len(peptide.rstrip('X'))
    features.append(actual_len / max_length)

    # Position-specific features for anchor residues (P2, P9 for MHC-I)
    if len(peptide) >= 2:
        p2 = peptide[1]
        features.append(1 if p2 in 'LMI' else 0)  # Hydrophobic anchor
        features.append(1 if p2 in 'TS' else 0)   # Polar anchor
    else:
        features.extend([0, 0])

    if len(peptide) >= 9:
        p9 = peptide[min(8, actual_len-1)]
        features.append(1 if p9 in 'VLI' else 0)  # C-terminal anchor
        features.append(1 if p9 in 'YF' else 0)   # Aromatic anchor
    else:
        features.extend([0, 0])

    return np.array(features, dtype=np.float32)


def parse_iedb_data(filepath: str, max_records: int = None) -> Dict[str, List[Tuple[str, float, str]]]:
    """
    Parse IEDB MHC ligand data.

    Returns:
        Dict[allele] -> List[(peptide, ic50_nm, assay_type)]
    """
    allele_data = defaultdict(list)

    print(f"Parsing {filepath}...")

    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)

        # Skip header rows
        next(reader)
        next(reader)

        count = 0
        skipped = 0

        for row in reader:
            if max_records and count >= max_records:
                break

            try:
                if len(row) < 108:
                    skipped += 1
                    continue

                peptide = row[COL_PEPTIDE].strip()
                allele = row[COL_MHC_ALLELE].strip()
                quant = row[COL_QUANT_MEASURE].strip()
                qual = row[COL_QUAL_MEASURE].strip().lower()
                assay = row[COL_ASSAY_TYPE].strip().lower() if len(row) > COL_ASSAY_TYPE else ''

                # Validate
                if not peptide or not allele:
                    skipped += 1
                    continue
                if not (8 <= len(peptide) <= 15):
                    skipped += 1
                    continue
                if not peptide.isalpha() or not peptide.isupper():
                    skipped += 1
                    continue
                if 'HLA' not in allele and 'H-2' not in allele:
                    skipped += 1
                    continue

                # Parse IC50 value
                ic50 = None
                if quant:
                    try:
                        ic50 = float(quant)
                    except:
                        pass

                # Fall back to qualitative
                if ic50 is None and qual:
                    for key, value in QUAL_TO_SCORE.items():
                        if key in qual:
                            ic50 = value
                            break

                if ic50 is None:
                    ic50 = 500  # Default moderate

                # Clamp to reasonable range
                ic50 = max(1, min(ic50, 100000))

                allele_data[allele].append((peptide, ic50, assay))
                count += 1

                if count % 100000 == 0:
                    print(f"  Parsed {count} records ({skipped} skipped)...")

            except Exception as e:
                skipped += 1
                continue

    print(f"Parsed {count} records, skipped {skipped}")
    print(f"Found {len(allele_data)} unique alleles")

    return dict(allele_data)


def train_mhc_model(allele: str, data: List[Tuple[str, float, str]]) -> Tuple[object, object, dict]:
    """Train RF model for one allele."""

    peptides = [d[0] for d in data]
    ic50s = np.array([d[1] for d in data])

    # Determine max peptide length for this allele
    max_len = max(len(p) for p in peptides)
    max_len = min(max_len, 12)  # Cap at 12

    # Encode
    X = np.array([encode_peptide(p, max_len) for p in peptides])

    # Log transform IC50 for better regression
    y = np.log10(ic50s)

    # Split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # Scale
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    # Train
    model = RandomForestRegressor(
        n_estimators=100,
        max_depth=20,
        min_samples_leaf=3,
        n_jobs=-1,
        random_state=42
    )
    model.fit(X_train_s, y_train)

    # Evaluate
    y_pred = model.predict(X_test_s)
    r2 = r2_score(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    # Also compute AUC for binder/non-binder classification
    # Threshold: IC50 < 500 nM = binder
    y_test_binary = (10**y_test < 500).astype(int)
    y_pred_binary = (10**y_pred < 500).astype(int)

    try:
        auc = roc_auc_score(y_test_binary, 1 - 10**y_pred)  # Higher score = lower IC50 = better binder
    except:
        auc = 0.5

    metrics = {
        'allele': allele,
        'n_samples': len(data),
        'r2': r2,
        'rmse': rmse,
        'auc': auc,
        'max_length': max_len,
    }

    print(f"  {allele}: n={len(data)}, R²={r2:.3f}, RMSE={rmse:.3f}, AUC={auc:.3f}")

    return model, scaler, metrics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', required=True, help='Path to mhc_ligand_full.csv')
    parser.add_argument('--output', default='weights/', help='Output directory')
    parser.add_argument('--max-records', type=int, default=1000000, help='Max records to parse')
    parser.add_argument('--min-samples', type=int, default=500, help='Min samples per allele')
    parser.add_argument('--top-alleles', type=int, default=30, help='Train top N alleles')
    args = parser.parse_args()

    if not SKLEARN_AVAILABLE:
        print("Error: scikit-learn required")
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)

    # Parse data
    allele_data = parse_iedb_data(args.data, args.max_records)

    # Filter and sort by sample count
    valid = {a: d for a, d in allele_data.items() if len(d) >= args.min_samples}
    sorted_alleles = sorted(valid.items(), key=lambda x: -len(x[1]))[:args.top_alleles]

    print(f"\nTraining models for {len(sorted_alleles)} alleles with >= {args.min_samples} samples")

    # Separate MHC-I and MHC-II
    mhc1_models = {}
    mhc2_models = {}

    for allele, data in sorted_alleles:
        print(f"\nTraining {allele}...")
        model, scaler, metrics = train_mhc_model(allele, data)

        model_data = {
            'model': model,
            'scaler': scaler,
            'metrics': metrics,
        }

        # Classify as MHC-I or MHC-II
        if 'DRB' in allele or 'DQB' in allele or 'DPB' in allele:
            mhc2_models[allele] = model_data
        else:
            mhc1_models[allele] = model_data

    # Save models
    mhc1_path = os.path.join(args.output, 'mhc1_iedb_models.pkl')
    mhc2_path = os.path.join(args.output, 'mhc2_iedb_models.pkl')

    with open(mhc1_path, 'wb') as f:
        pickle.dump(mhc1_models, f)
    print(f"\nSaved {len(mhc1_models)} MHC-I models to {mhc1_path}")

    with open(mhc2_path, 'wb') as f:
        pickle.dump(mhc2_models, f)
    print(f"Saved {len(mhc2_models)} MHC-II models to {mhc2_path}")

    # Summary
    print("\n" + "="*60)
    print("TRAINING SUMMARY")
    print("="*60)

    print("\nMHC Class I:")
    for allele, data in sorted(mhc1_models.items(), key=lambda x: -x[1]['metrics']['n_samples']):
        m = data['metrics']
        print(f"  {allele}: n={m['n_samples']}, R²={m['r2']:.3f}, AUC={m['auc']:.3f}")

    print("\nMHC Class II:")
    for allele, data in sorted(mhc2_models.items(), key=lambda x: -x[1]['metrics']['n_samples']):
        m = data['metrics']
        print(f"  {allele}: n={m['n_samples']}, R²={m['r2']:.3f}, AUC={m['auc']:.3f}")

    avg_r2_1 = np.mean([d['metrics']['r2'] for d in mhc1_models.values()]) if mhc1_models else 0
    avg_r2_2 = np.mean([d['metrics']['r2'] for d in mhc2_models.values()]) if mhc2_models else 0
    avg_auc_1 = np.mean([d['metrics']['auc'] for d in mhc1_models.values()]) if mhc1_models else 0
    avg_auc_2 = np.mean([d['metrics']['auc'] for d in mhc2_models.values()]) if mhc2_models else 0

    print(f"\nMHC-I average: R²={avg_r2_1:.3f}, AUC={avg_auc_1:.3f}")
    print(f"MHC-II average: R²={avg_r2_2:.3f}, AUC={avg_auc_2:.3f}")


if __name__ == '__main__':
    main()
