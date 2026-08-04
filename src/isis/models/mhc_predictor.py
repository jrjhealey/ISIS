"""
MHC Binding Prediction using IEDB-trained models.

This module provides MHC-I and MHC-II binding prediction using
Random Forest models trained on real IEDB binding data.
"""

import os
import pickle
from typing import List, Dict, Optional, Tuple
import numpy as np

# Try to import sklearn
try:
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


# Amino acid encoding
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

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


def encode_peptide(peptide: str, max_length: int = 12) -> np.ndarray:
    """Encode peptide for model input."""
    peptide = peptide[:max_length].ljust(max_length, 'X')
    features = []

    # One-hot
    for aa in peptide:
        onehot = [0] * 20
        if aa in AA_TO_IDX:
            onehot[AA_TO_IDX[aa]] = 1
        features.extend(onehot)

    # BLOSUM62
    for aa in peptide:
        if aa in BLOSUM62:
            features.extend(BLOSUM62[aa])
        else:
            features.extend([0] * 20)

    # Length
    actual_len = len(peptide.rstrip('X'))
    features.append(actual_len / max_length)

    # Anchor features
    if len(peptide) >= 2:
        p2 = peptide[1]
        features.append(1 if p2 in 'LMI' else 0)
        features.append(1 if p2 in 'TS' else 0)
    else:
        features.extend([0, 0])

    if len(peptide) >= 9:
        p9 = peptide[min(8, actual_len-1)]
        features.append(1 if p9 in 'VLI' else 0)
        features.append(1 if p9 in 'YF' else 0)
    else:
        features.extend([0, 0])

    return np.array(features, dtype=np.float32)


class MHCPredictor:
    """
    MHC binding prediction using IEDB-trained models.

    Supports both MHC Class I and Class II alleles.
    Models are loaded lazily on first use.
    """

    def __init__(self):
        self.mhc1_models = None
        self.mhc2_models = None
        self._weights_dir = os.path.join(os.path.dirname(__file__), 'weights')

    def _load_models(self, mhc_class: int):
        """Load models for MHC class."""
        if mhc_class == 1 and self.mhc1_models is None:
            path = os.path.join(self._weights_dir, 'mhc1_iedb_models.pkl')
            if os.path.exists(path):
                with open(path, 'rb') as f:
                    self.mhc1_models = pickle.load(f)
            else:
                self.mhc1_models = {}

        elif mhc_class == 2 and self.mhc2_models is None:
            path = os.path.join(self._weights_dir, 'mhc2_iedb_models.pkl')
            if os.path.exists(path):
                with open(path, 'rb') as f:
                    self.mhc2_models = pickle.load(f)
            else:
                self.mhc2_models = {}

    def available_alleles(self, mhc_class: int = None) -> List[str]:
        """Return list of available alleles."""
        alleles = []

        if mhc_class is None or mhc_class == 1:
            self._load_models(1)
            alleles.extend(self.mhc1_models.keys() if self.mhc1_models else [])

        if mhc_class is None or mhc_class == 2:
            self._load_models(2)
            alleles.extend(self.mhc2_models.keys() if self.mhc2_models else [])

        return sorted(alleles)

    def predict_mhc1(
        self,
        sequence: str,
        allele: str = 'HLA-A*02:01',
        peptide_length: int = 9
    ) -> Dict:
        """
        Predict MHC Class I binding for peptides in a sequence.

        Args:
            sequence: Protein sequence
            allele: MHC-I allele (e.g., 'HLA-A*02:01')
            peptide_length: Length of peptides to generate (8-11)

        Returns:
            Dict with 'peptides' (list of predictions) and 'allele'
        """
        self._load_models(1)

        if not self.mhc1_models:
            raise ValueError("No MHC-I models available. Train models first.")

        if allele not in self.mhc1_models:
            available = list(self.mhc1_models.keys())[:5]
            raise ValueError(f"Allele {allele} not available. Try: {available}")

        model_data = self.mhc1_models[allele]
        model = model_data['model']
        scaler = model_data['scaler']
        max_len = model_data['metrics'].get('max_length', 12)

        # Generate peptides
        peptides = []
        for i in range(len(sequence) - peptide_length + 1):
            pep = sequence[i:i + peptide_length]
            if all(aa in AMINO_ACIDS for aa in pep):
                peptides.append({
                    'peptide': pep,
                    'start': i + 1,
                    'end': i + peptide_length
                })

        if not peptides:
            return {'allele': allele, 'peptides': []}

        # Encode and predict
        X = np.array([encode_peptide(p['peptide'], max_len) for p in peptides])
        X_scaled = scaler.transform(X)
        log_ic50 = model.predict(X_scaled)

        # Convert to IC50 and add to results
        for i, pep in enumerate(peptides):
            ic50 = 10 ** log_ic50[i]
            pep['ic50'] = float(ic50)
            pep['log_ic50'] = float(log_ic50[i])
            pep['score'] = float(max(0, 1 - log_ic50[i] / 5))  # Normalize to 0-1

            # Classification
            if ic50 < 50:
                pep['binding'] = 'strong'
            elif ic50 < 500:
                pep['binding'] = 'weak'
            else:
                pep['binding'] = 'none'

        # Sort by IC50
        peptides.sort(key=lambda x: x['ic50'])

        # Calculate percentile ranks
        ic50s = [p['ic50'] for p in peptides]
        for i, pep in enumerate(peptides):
            pep['percentile'] = 100.0 * (i + 1) / len(peptides)

        return {
            'allele': allele,
            'peptides': peptides,
            'model_metrics': model_data['metrics']
        }

    def predict_mhc2(
        self,
        sequence: str,
        allele: str = 'HLA-DRB1*01:01'
    ) -> Dict:
        """
        Predict MHC Class II binding.

        MHC-II binds 15-mer peptides with a 9-mer binding core.

        Args:
            sequence: Protein sequence
            allele: MHC-II allele (e.g., 'HLA-DRB1*01:01')

        Returns:
            Dict with peptide predictions including binding core
        """
        self._load_models(2)

        if not self.mhc2_models:
            raise ValueError("No MHC-II models available. Train models first.")

        if allele not in self.mhc2_models:
            available = list(self.mhc2_models.keys())[:5]
            raise ValueError(f"Allele {allele} not available. Try: {available}")

        model_data = self.mhc2_models[allele]
        model = model_data['model']
        scaler = model_data['scaler']
        max_len = model_data['metrics'].get('max_length', 15)

        # Generate 15-mer peptides
        peptide_length = 15
        peptides = []

        for i in range(len(sequence) - peptide_length + 1):
            pep = sequence[i:i + peptide_length]
            if all(aa in AMINO_ACIDS for aa in pep):
                peptides.append({
                    'peptide': pep,
                    'start': i + 1,
                    'end': i + peptide_length
                })

        if not peptides:
            return {'allele': allele, 'peptides': []}

        # Encode and predict
        X = np.array([encode_peptide(p['peptide'], max_len) for p in peptides])
        X_scaled = scaler.transform(X)
        log_ic50 = model.predict(X_scaled)

        # Add predictions
        for i, pep in enumerate(peptides):
            ic50 = 10 ** log_ic50[i]
            pep['ic50'] = float(ic50)
            pep['log_ic50'] = float(log_ic50[i])
            pep['score'] = float(max(0, 1 - log_ic50[i] / 5))

            # Find binding core (9-mer with best anchor fit)
            best_core = self._find_binding_core(pep['peptide'], allele)
            pep['core'] = best_core['core']
            pep['core_start'] = best_core['position']

            if ic50 < 50:
                pep['binding'] = 'strong'
            elif ic50 < 500:
                pep['binding'] = 'weak'
            else:
                pep['binding'] = 'none'

        peptides.sort(key=lambda x: x['ic50'])

        return {
            'allele': allele,
            'peptides': peptides,
            'model_metrics': model_data['metrics']
        }

    def _find_binding_core(self, peptide: str, allele: str) -> Dict:
        """Find 9-mer binding core within a 15-mer."""
        # Anchor preferences by allele
        if 'DRB1*01:01' in allele:
            p1_pref = set('YFWLIV')
            p4_pref = set('AGST')
            p9_pref = set('LMIVF')
        elif 'DRB1*04:01' in allele:
            p1_pref = set('YFW')
            p4_pref = set('ST')
            p9_pref = set('KR')
        else:
            p1_pref = set('YFWLIV')
            p4_pref = set('AGST')
            p9_pref = set('LMIVF')

        best_score = -1
        best_core = peptide[:9]
        best_pos = 0

        for i in range(len(peptide) - 8):
            core = peptide[i:i+9]
            score = 0
            if core[0] in p1_pref:
                score += 1
            if core[3] in p4_pref:
                score += 1
            if core[8] in p9_pref:
                score += 1

            if score > best_score:
                best_score = score
                best_core = core
                best_pos = i

        return {'core': best_core, 'position': best_pos}

    def predict_from_protein(
        self,
        sequence: str,
        alleles: List[str] = None,
        mhc_class: int = 1,
        top_n: int = 100
    ) -> Dict:
        """
        Predict binding for multiple alleles across a protein.

        Args:
            sequence: Protein sequence
            alleles: List of alleles (default: all available)
            mhc_class: 1 or 2
            top_n: Return top N binders per allele

        Returns:
            Dict with per-allele predictions
        """
        if alleles is None:
            alleles = self.available_alleles(mhc_class)[:5]  # Top 5 by default

        results = {}
        for allele in alleles:
            try:
                if mhc_class == 1:
                    pred = self.predict_mhc1(sequence, allele)
                else:
                    pred = self.predict_mhc2(sequence, allele)

                # Keep only top binders
                pred['peptides'] = pred['peptides'][:top_n]
                results[allele] = pred

            except Exception as e:
                results[allele] = {'error': str(e)}

        return results

    def get_binders(
        self,
        predictions: Dict,
        ic50_threshold: float = 500
    ) -> List[Dict]:
        """Extract binders below IC50 threshold."""
        binders = []
        for pep in predictions.get('peptides', []):
            if pep.get('ic50', float('inf')) < ic50_threshold:
                binders.append(pep)
        return binders


# Convenience functions
_predictor = None

def get_predictor() -> MHCPredictor:
    """Get singleton predictor instance."""
    global _predictor
    if _predictor is None:
        _predictor = MHCPredictor()
    return _predictor


def predict_mhc1(sequence: str, allele: str = 'HLA-A*02:01', **kwargs) -> Dict:
    """Predict MHC-I binding."""
    return get_predictor().predict_mhc1(sequence, allele, **kwargs)


def predict_mhc2(sequence: str, allele: str = 'HLA-DRB1*01:01', **kwargs) -> Dict:
    """Predict MHC-II binding."""
    return get_predictor().predict_mhc2(sequence, allele, **kwargs)


def available_alleles(mhc_class: int = None) -> List[str]:
    """Get available alleles."""
    return get_predictor().available_alleles(mhc_class)
