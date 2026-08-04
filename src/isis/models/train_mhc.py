#!/usr/bin/env python3
"""
Train MHC binding prediction models on real IEDB data.

Usage:
    python train_mhc.py --data /path/to/mhc_ligand_full.csv --output models/
"""

import argparse
import csv
import os
import pickle
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import numpy as np

try:
    from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
    from sklearn.neural_network import MLPRegressor
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.metrics import mean_squared_error, r2_score
    import joblib
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    print("Warning: scikit-learn not available. Install with: pip install scikit-learn")


# Amino acid encoding
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# BLOSUM62 matrix for encoding
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


def encode_peptide_onehot(peptide: str, max_length: int = 15) -> np.ndarray:
    """One-hot encode a peptide sequence."""
    encoding = np.zeros((max_length, 20))
    for i, aa in enumerate(peptide[:max_length]):
        if aa in AA_TO_IDX:
            encoding[i, AA_TO_IDX[aa]] = 1
    return encoding.flatten()


def encode_peptide_blosum(peptide: str, max_length: int = 15) -> np.ndarray:
    """BLOSUM62 encode a peptide sequence."""
    encoding = np.zeros((max_length, 20))
    for i, aa in enumerate(peptide[:max_length]):
        if aa in BLOSUM62:
            encoding[i] = BLOSUM62[aa]
    return encoding.flatten()


def encode_peptide_combined(peptide: str, max_length: int = 15) -> np.ndarray:
    """Combined one-hot + BLOSUM encoding."""
    onehot = encode_peptide_onehot(peptide, max_length)
    blosum = encode_peptide_blosum(peptide, max_length)
    length_feature = np.array([len(peptide) / max_length])
    return np.concatenate([onehot, blosum, length_feature])


def parse_iedb_mhc_data(filepath: str, max_records: int = None) -> Dict[str, List[Tuple[str, float]]]:
    """
    Parse IEDB MHC ligand data file.

    Returns:
        Dict mapping allele name -> list of (peptide, ic50_nm) tuples
    """
    allele_data = defaultdict(list)

    print(f"Parsing {filepath}...")

    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        reader = csv.reader(f)

        # Read header
        header = next(reader)

        # Find relevant column indices
        # We need: peptide sequence, MHC allele, quantitative measurement
        peptide_col = None
        allele_col = None
        measurement_col = None
        measurement_type_col = None

        for i, col in enumerate(header):
            col_lower = col.lower()
            if 'description' in col_lower and peptide_col is None:
                peptide_col = i
            elif 'name' in col_lower and 'epitope' in header[i-1].lower() if i > 0 else False:
                peptide_col = i
            elif 'allele' in col_lower and 'name' in col_lower:
                allele_col = i
            elif 'quantitative' in col_lower and 'measurement' in col_lower:
                measurement_col = i
            elif 'measurement' in col_lower and 'inequality' not in col_lower:
                if measurement_col is None:
                    measurement_col = i

        # Fallback: search by position for known IEDB format
        # IEDB format has epitope name around column 10-12, allele around 45-50
        if peptide_col is None:
            for i, col in enumerate(header):
                if col == 'Name' or col == 'Description':
                    peptide_col = i
                    break

        print(f"Found columns - Peptide: {peptide_col}, Allele: {allele_col}, Measurement: {measurement_col}")

        # If we can't find columns, try a different approach
        if peptide_col is None or allele_col is None:
            print("Using positional parsing for IEDB format...")
            # Reset and use known positions
            f.seek(0)
            next(f)  # skip header

            records = 0
            for line in f:
                if max_records and records >= max_records:
                    break

                try:
                    parts = line.split(',')
                    if len(parts) < 50:
                        continue

                    # IEDB format: epitope name is typically around position 10-12
                    # Look for a peptide-like string (uppercase letters, 8-15 chars)
                    peptide = None
                    allele = None
                    measurement = None

                    for part in parts[8:20]:
                        part = part.strip().strip('"')
                        if part and 8 <= len(part) <= 15 and part.isalpha() and part.isupper():
                            peptide = part
                            break

                    # Look for HLA allele
                    for part in parts[40:60]:
                        part = part.strip().strip('"')
                        if 'HLA-' in part or 'H-2' in part:
                            allele = part
                            break

                    # Look for IC50 measurement (numeric value)
                    for part in parts[60:80]:
                        part = part.strip().strip('"')
                        try:
                            val = float(part)
                            if 0.1 <= val <= 100000:  # Reasonable IC50 range
                                measurement = val
                                break
                        except:
                            continue

                    if peptide and allele and measurement:
                        allele_data[allele].append((peptide, measurement))
                        records += 1

                        if records % 100000 == 0:
                            print(f"  Parsed {records} records...")

                except Exception as e:
                    continue

            print(f"Parsed {records} total records")
            return dict(allele_data)

        # Standard CSV parsing
        records = 0
        for row in reader:
            if max_records and records >= max_records:
                break

            try:
                if len(row) <= max(peptide_col or 0, allele_col or 0, measurement_col or 0):
                    continue

                peptide = row[peptide_col].strip().strip('"') if peptide_col else None
                allele = row[allele_col].strip().strip('"') if allele_col else None

                # Get measurement
                measurement = None
                if measurement_col and row[measurement_col]:
                    try:
                        measurement = float(row[measurement_col].strip().strip('"'))
                    except:
                        pass

                # Validate
                if not peptide or not allele:
                    continue
                if not (8 <= len(peptide) <= 15):
                    continue
                if not peptide.isalpha():
                    continue
                if not ('HLA' in allele or 'H-2' in allele):
                    continue

                # Use default measurement if not available
                if measurement is None:
                    measurement = 500  # Default moderate binder

                allele_data[allele].append((peptide, measurement))
                records += 1

                if records % 100000 == 0:
                    print(f"  Parsed {records} records...")

            except Exception as e:
                continue

    print(f"Parsed {records} total records across {len(allele_data)} alleles")
    return dict(allele_data)


def train_mhc_model(
    allele: str,
    peptides: List[str],
    ic50_values: List[float],
    model_type: str = 'rf'
) -> Tuple[object, object, dict]:
    """
    Train an MHC binding prediction model for one allele.

    Args:
        allele: MHC allele name
        peptides: List of peptide sequences
        ic50_values: List of IC50 values in nM
        model_type: 'rf' for RandomForest, 'mlp' for neural network

    Returns:
        (model, scaler, metrics)
    """
    print(f"Training model for {allele} with {len(peptides)} peptides...")

    # Encode peptides
    max_len = max(len(p) for p in peptides)
    max_len = min(max_len, 15)  # Cap at 15

    X = np.array([encode_peptide_combined(p, max_len) for p in peptides])

    # Convert IC50 to log scale (more stable for training)
    y = np.log10(np.clip(ic50_values, 1, 100000))

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )

    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Train model
    if model_type == 'rf':
        model = RandomForestRegressor(
            n_estimators=100,
            max_depth=15,
            min_samples_leaf=5,
            n_jobs=-1,
            random_state=42
        )
    elif model_type == 'gb':
        model = GradientBoostingRegressor(
            n_estimators=100,
            max_depth=8,
            learning_rate=0.1,
            random_state=42
        )
    else:  # mlp
        model = MLPRegressor(
            hidden_layer_sizes=(256, 128, 64),
            activation='relu',
            solver='adam',
            max_iter=500,
            early_stopping=True,
            random_state=42
        )

    model.fit(X_train_scaled, y_train)

    # Evaluate
    y_pred = model.predict(X_test_scaled)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    metrics = {
        'allele': allele,
        'n_samples': len(peptides),
        'mse': mse,
        'rmse': np.sqrt(mse),
        'r2': r2,
        'max_peptide_length': max_len
    }

    print(f"  {allele}: R² = {r2:.3f}, RMSE = {np.sqrt(mse):.3f} (log10 IC50)")

    return model, scaler, metrics


def main():
    parser = argparse.ArgumentParser(description='Train MHC binding models on IEDB data')
    parser.add_argument('--data', required=True, help='Path to mhc_ligand_full.csv')
    parser.add_argument('--output', default='models/', help='Output directory for models')
    parser.add_argument('--max-records', type=int, default=500000, help='Max records to parse')
    parser.add_argument('--min-samples', type=int, default=1000, help='Min samples per allele')
    parser.add_argument('--model-type', default='rf', choices=['rf', 'gb', 'mlp'])
    args = parser.parse_args()

    if not SKLEARN_AVAILABLE:
        print("Error: scikit-learn is required. Install with: pip install scikit-learn")
        sys.exit(1)

    os.makedirs(args.output, exist_ok=True)

    # Parse data
    allele_data = parse_iedb_mhc_data(args.data, args.max_records)

    # Filter alleles with enough data
    valid_alleles = {
        allele: data for allele, data in allele_data.items()
        if len(data) >= args.min_samples
    }

    print(f"\nFound {len(valid_alleles)} alleles with >= {args.min_samples} samples")

    # Train models for top alleles
    all_metrics = []
    models = {}

    # Sort by sample count, train top alleles
    sorted_alleles = sorted(valid_alleles.items(), key=lambda x: -len(x[1]))

    for allele, data in sorted_alleles[:20]:  # Top 20 alleles
        peptides = [d[0] for d in data]
        ic50s = [d[1] for d in data]

        model, scaler, metrics = train_mhc_model(
            allele, peptides, ic50s, args.model_type
        )

        models[allele] = {
            'model': model,
            'scaler': scaler,
            'metrics': metrics
        }
        all_metrics.append(metrics)

    # Save models
    output_file = os.path.join(args.output, 'mhc1_models.pkl')
    with open(output_file, 'wb') as f:
        pickle.dump(models, f)

    print(f"\nSaved models to {output_file}")

    # Summary
    print("\nTraining Summary:")
    print("-" * 60)
    for m in all_metrics:
        print(f"  {m['allele']}: {m['n_samples']} samples, R²={m['r2']:.3f}")

    avg_r2 = np.mean([m['r2'] for m in all_metrics])
    print(f"\nAverage R²: {avg_r2:.3f}")


if __name__ == '__main__':
    main()
