"""
B-cell linear epitope prediction model (BepiPred v2.0 style).

This module implements a Random Forest-based B-cell linear epitope predictor
using engineered features from amino acid property scales. The approach is
modeled after BepiPred v2.0 which achieved state-of-the-art performance
using Random Forest on biophysical features.

Key features:
- One-hot amino acid encoding (20 features)
- Multiple biophysical property scales
- Sliding window context for each position
- Random Forest classifier (robust, interpretable)
- Pre-trained on synthetic data capturing biological signal

References:
    Jespersen MC et al. BepiPred-2.0: improving sequence-based B-cell epitope
    prediction using conformational epitopes. Nucleic Acids Res. 2017.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from numpy.random import RandomState

try:
    import joblib
    JOBLIB_AVAILABLE = True
except ImportError:
    JOBLIB_AVAILABLE = False

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.preprocessing import StandardScaler
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


# ============================================================================
# Amino Acid Property Scales
# ============================================================================

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INDEX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# Kyte-Doolittle hydrophobicity scale
# Kyte J, Doolittle RF. J Mol Biol. 1982;157(1):105-132. PMID: 7108955
KYTE_DOOLITTLE = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

# Emini surface accessibility scale
# Emini EA et al. J Virol. 1985;55(3):836-839. PMID: 2410642
EMINI = {
    "A": 0.815, "C": 0.394, "D": 1.283, "E": 1.445, "F": 0.695,
    "G": 0.714, "H": 1.180, "I": 0.603, "K": 1.545, "L": 0.603,
    "M": 0.714, "N": 1.296, "P": 1.236, "Q": 1.348, "R": 1.475,
    "S": 1.115, "T": 1.184, "V": 0.606, "W": 0.808, "Y": 1.089,
}

# Karplus-Schulz flexibility scale
# Karplus PA, Schulz GE. Naturwissenschaften. 1985;72:212-213.
KARPLUS_SCHULZ = {
    "A": 1.046, "C": 0.906, "D": 1.068, "E": 1.094, "F": 0.915,
    "G": 1.061, "H": 0.950, "I": 0.927, "K": 1.102, "L": 0.935,
    "M": 0.893, "N": 1.036, "P": 1.049, "Q": 1.037, "R": 1.008,
    "S": 1.046, "T": 0.997, "V": 0.912, "W": 0.904, "Y": 0.929,
}

# Chou-Fasman beta-turn propensity
# Chou PY, Fasman GD. Adv Enzymol. 1978;47:45-148. PMID: 364941
CHOU_FASMAN_TURN = {
    "A": 0.66, "C": 1.19, "D": 1.46, "E": 0.74, "F": 0.60,
    "G": 1.56, "H": 0.95, "I": 0.47, "K": 1.01, "L": 0.59,
    "M": 0.60, "N": 1.56, "P": 1.52, "Q": 0.98, "R": 0.95,
    "S": 1.43, "T": 0.96, "V": 0.50, "W": 0.96, "Y": 1.14,
}

# Kolaskar-Tongaonkar antigenicity scale
# Kolaskar AS, Tongaonkar PC. FEBS Lett. 1990;276(1-2):172-174. PMID: 1702393
KOLASKAR_TONGAONKAR = {
    "A": 1.064, "C": 1.412, "D": 0.866, "E": 0.851, "F": 1.091,
    "G": 0.874, "H": 1.105, "I": 1.152, "K": 0.930, "L": 1.250,
    "M": 0.826, "N": 0.776, "P": 1.064, "Q": 1.015, "R": 0.873,
    "S": 1.012, "T": 0.909, "V": 1.383, "W": 0.893, "Y": 1.161,
}

# Chou-Fasman alpha-helix propensity
CHOU_FASMAN_HELIX = {
    "A": 1.42, "C": 0.70, "D": 1.01, "E": 1.51, "F": 1.13,
    "G": 0.57, "H": 1.00, "I": 1.08, "K": 1.16, "L": 1.21,
    "M": 1.45, "N": 0.67, "P": 0.57, "Q": 1.11, "R": 0.98,
    "S": 0.77, "T": 0.83, "V": 1.06, "W": 1.08, "Y": 0.69,
}

# Chou-Fasman beta-sheet propensity
CHOU_FASMAN_SHEET = {
    "A": 0.83, "C": 1.19, "D": 0.54, "E": 0.37, "F": 1.38,
    "G": 0.75, "H": 0.87, "I": 1.60, "K": 0.74, "L": 1.30,
    "M": 1.05, "N": 0.89, "P": 0.55, "Q": 1.10, "R": 0.93,
    "S": 0.75, "T": 1.19, "V": 1.70, "W": 1.37, "Y": 1.47,
}

# Amino acid charge at pH 7
AA_CHARGE = {
    "A": 0.0, "C": 0.0, "D": -1.0, "E": -1.0, "F": 0.0,
    "G": 0.0, "H": 0.1, "I": 0.0, "K": 1.0, "L": 0.0,  # H is ~10% charged at pH 7
    "M": 0.0, "N": 0.0, "P": 0.0, "Q": 0.0, "R": 1.0,
    "S": 0.0, "T": 0.0, "V": 0.0, "W": 0.0, "Y": 0.0,
}

# Amino acid volume (Zamyatnin, 1984)
AA_VOLUME = {
    "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9,
    "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7,
    "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4,
    "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6,
}

# Polarity (Grantham, 1974)
AA_POLARITY = {
    "A": 8.1, "C": 5.5, "D": 13.0, "E": 12.3, "F": 5.2,
    "G": 9.0, "H": 10.4, "I": 5.2, "K": 11.3, "L": 4.9,
    "M": 5.7, "N": 11.6, "P": 8.0, "Q": 10.5, "R": 10.5,
    "S": 9.2, "T": 8.6, "V": 5.9, "W": 5.4, "Y": 6.2,
}

# Number of features per position (before window aggregation)
N_BASE_FEATURES = 29  # 20 one-hot + 9 scales


def _normalize_scale(scale: Dict[str, float]) -> Dict[str, float]:
    """Normalize a scale to zero mean and unit variance."""
    values = np.array(list(scale.values()))
    mean, std = values.mean(), values.std()
    if std == 0:
        std = 1.0
    return {k: (v - mean) / std for k, v in scale.items()}


# Normalized scales for feature extraction
_NORM_HYDRO = _normalize_scale(KYTE_DOOLITTLE)
_NORM_EMINI = _normalize_scale(EMINI)
_NORM_FLEX = _normalize_scale(KARPLUS_SCHULZ)
_NORM_TURN = _normalize_scale(CHOU_FASMAN_TURN)
_NORM_ANTIG = _normalize_scale(KOLASKAR_TONGAONKAR)
_NORM_HELIX = _normalize_scale(CHOU_FASMAN_HELIX)
_NORM_SHEET = _normalize_scale(CHOU_FASMAN_SHEET)
_NORM_VOL = _normalize_scale(AA_VOLUME)
_NORM_POL = _normalize_scale(AA_POLARITY)


# ============================================================================
# Feature Engineering
# ============================================================================

def one_hot_encode(aa: str) -> np.ndarray:
    """
    One-hot encode a single amino acid.

    Args:
        aa: Single amino acid character.

    Returns:
        20-element array with 1 at the amino acid's index, 0 elsewhere.
    """
    encoding = np.zeros(20, dtype=np.float32)
    if aa in AA_TO_INDEX:
        encoding[AA_TO_INDEX[aa]] = 1.0
    return encoding


def get_residue_features(aa: str) -> np.ndarray:
    """
    Extract all features for a single amino acid residue.

    Features (29 total):
    - One-hot encoding (20)
    - Hydrophobicity (1)
    - Surface accessibility (1)
    - Flexibility (1)
    - Beta-turn propensity (1)
    - Antigenicity (1)
    - Alpha-helix propensity (1)
    - Beta-sheet propensity (1)
    - Charge (1)
    - Volume/size (1)

    Args:
        aa: Single amino acid character.

    Returns:
        Feature vector of length 29.
    """
    features = np.zeros(N_BASE_FEATURES, dtype=np.float32)

    # One-hot encoding (positions 0-19)
    if aa in AA_TO_INDEX:
        features[AA_TO_INDEX[aa]] = 1.0

    # Biophysical properties (positions 20-28)
    features[20] = _NORM_HYDRO.get(aa, 0.0)
    features[21] = _NORM_EMINI.get(aa, 0.0)
    features[22] = _NORM_FLEX.get(aa, 0.0)
    features[23] = _NORM_TURN.get(aa, 0.0)
    features[24] = _NORM_ANTIG.get(aa, 0.0)
    features[25] = _NORM_HELIX.get(aa, 0.0)
    features[26] = _NORM_SHEET.get(aa, 0.0)
    features[27] = AA_CHARGE.get(aa, 0.0)  # Already normalized (-1 to 1)
    features[28] = _NORM_VOL.get(aa, 0.0)

    return features


def extract_window_features(
    sequence: str,
    position: int,
    window_size: int = 9
) -> np.ndarray:
    """
    Extract features for a position with window context.

    The feature vector includes:
    - Features for the center residue (29)
    - Mean features for left neighbors (29)
    - Mean features for right neighbors (29)
    - Positional features (2): relative position in sequence, distance from termini

    Total: 89 features

    Args:
        sequence: Protein sequence.
        position: Center position (0-indexed).
        window_size: Size of context window (must be odd).

    Returns:
        Feature vector of length 89.
    """
    seq_len = len(sequence)
    half_window = window_size // 2

    # Get features for center residue
    center_features = get_residue_features(sequence[position])

    # Get features for left context (average)
    left_features = np.zeros(N_BASE_FEATURES, dtype=np.float32)
    left_count = 0
    for i in range(max(0, position - half_window), position):
        left_features += get_residue_features(sequence[i])
        left_count += 1
    if left_count > 0:
        left_features /= left_count

    # Get features for right context (average)
    right_features = np.zeros(N_BASE_FEATURES, dtype=np.float32)
    right_count = 0
    for i in range(position + 1, min(seq_len, position + half_window + 1)):
        right_features += get_residue_features(sequence[i])
        right_count += 1
    if right_count > 0:
        right_features /= right_count

    # Positional features
    rel_position = position / max(1, seq_len - 1)  # 0 to 1
    dist_from_termini = min(position, seq_len - 1 - position) / max(1, (seq_len - 1) / 2)

    # Combine all features
    return np.concatenate([
        center_features,       # 29 features
        left_features,         # 29 features
        right_features,        # 29 features
        [rel_position, dist_from_termini]  # 2 features
    ])  # Total: 89 features


# ============================================================================
# Synthetic Training Data Generation
# ============================================================================

def generate_synthetic_epitope_sequence(
    length: int,
    rng: RandomState
) -> str:
    """
    Generate a synthetic epitope-like sequence.

    Epitopes tend to be:
    - Hydrophilic (polar/charged residues)
    - Flexible (G, S, N, D)
    - In turns (P, G, N, S)
    - Surface exposed (K, R, E, D, N, Q)

    Args:
        length: Sequence length.
        rng: Random state for reproducibility.

    Returns:
        Generated amino acid sequence.
    """
    # Epitope-favored amino acids with weights
    epitope_aa = "KRENDQSGPT"  # Charged, polar, flexible, turn-forming
    other_aa = "ACFHILMVWY"

    # Higher probability for epitope-like residues
    seq = []
    for _ in range(length):
        if rng.random() < 0.75:  # 75% chance of epitope-favored AA
            seq.append(rng.choice(list(epitope_aa)))
        else:
            seq.append(rng.choice(list(other_aa)))

    return "".join(seq)


def generate_synthetic_nonepitope_sequence(
    length: int,
    rng: RandomState
) -> str:
    """
    Generate a synthetic non-epitope (buried/structured) sequence.

    Non-epitopes tend to be:
    - Hydrophobic (buried)
    - Structured (helices, sheets)
    - Low flexibility

    Args:
        length: Sequence length.
        rng: Random state for reproducibility.

    Returns:
        Generated amino acid sequence.
    """
    # Non-epitope favored amino acids (hydrophobic, structured)
    buried_aa = "AVILMFYW"  # Hydrophobic
    helix_aa = "AELM"      # High helix propensity
    sheet_aa = "VIFY"      # High sheet propensity

    seq = []
    for _ in range(length):
        choice = rng.random()
        if choice < 0.5:  # 50% buried hydrophobic
            seq.append(rng.choice(list(buried_aa)))
        elif choice < 0.75:  # 25% helix-forming
            seq.append(rng.choice(list(helix_aa)))
        else:  # 25% sheet-forming
            seq.append(rng.choice(list(sheet_aa)))

    return "".join(seq)


def generate_training_data(
    n_samples: int = 50000,
    window_size: int = 9,
    seed: int = 42
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate synthetic training data for the B-cell epitope model.

    Creates balanced dataset with:
    - Epitope samples: hydrophilic, flexible, surface-exposed, turns
    - Non-epitope samples: hydrophobic, buried, structured (helices/sheets)

    Args:
        n_samples: Total number of training samples.
        window_size: Context window size.
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (features, labels) arrays.
    """
    rng = RandomState(seed)
    n_positive = n_samples // 2
    n_negative = n_samples - n_positive

    # Preallocate arrays
    n_features = 89  # From extract_window_features
    X = np.zeros((n_samples, n_features), dtype=np.float32)
    y = np.zeros(n_samples, dtype=np.int32)

    idx = 0

    # Generate epitope samples (positive class)
    seq_len_min, seq_len_max = 15, 50
    samples_per_seq = 5  # Extract multiple positions per sequence

    while idx < n_positive:
        seq_len = rng.randint(seq_len_min, seq_len_max)
        seq = generate_synthetic_epitope_sequence(seq_len, rng)

        # Sample positions from this sequence
        n_to_sample = min(samples_per_seq, n_positive - idx, seq_len)
        positions = rng.choice(seq_len, n_to_sample, replace=False)

        for pos in positions:
            X[idx] = extract_window_features(seq, pos, window_size)
            y[idx] = 1
            idx += 1
            if idx >= n_positive:
                break

    # Generate non-epitope samples (negative class)
    while idx < n_samples:
        seq_len = rng.randint(seq_len_min, seq_len_max)
        seq = generate_synthetic_nonepitope_sequence(seq_len, rng)

        # Sample positions from this sequence
        n_to_sample = min(samples_per_seq, n_samples - idx, seq_len)
        positions = rng.choice(seq_len, n_to_sample, replace=False)

        for pos in positions:
            X[idx] = extract_window_features(seq, pos, window_size)
            y[idx] = 0
            idx += 1
            if idx >= n_samples:
                break

    # Shuffle the data
    perm = rng.permutation(n_samples)
    X = X[perm]
    y = y[perm]

    return X, y


# ============================================================================
# B-cell Linear Epitope Predictor
# ============================================================================

class BcellLinearPredictor:
    """
    Random Forest-based B-cell linear epitope predictor.

    This predictor uses engineered features from multiple amino acid property
    scales combined with sliding window context. The model architecture is
    inspired by BepiPred v2.0.

    Features (89 per position):
    - One-hot amino acid encoding (20)
    - Biophysical properties: hydrophobicity, accessibility, flexibility,
      beta-turn, antigenicity, helix propensity, sheet propensity, charge,
      volume (9)
    - Left context average (29)
    - Right context average (29)
    - Positional encoding (2)

    Attributes:
        window_size: Context window size for feature extraction.
        model: Trained RandomForestClassifier.
        scaler: StandardScaler for feature normalization.
        is_trained: Whether the model has been trained.

    Example:
        >>> predictor = BcellLinearPredictor()
        >>> predictor.train()
        >>> probabilities = predictor.predict("MKTAYIAKQRQISFVKSHFSRQLE")
        >>> epitopes = predictor.predict_epitopes("MKTAYIAKQRQISFVKSHFSRQLE")
    """

    def __init__(self, window_size: int = 9):
        """
        Initialize the predictor.

        Args:
            window_size: Context window size (must be odd, default 9).
        """
        if window_size % 2 == 0:
            raise ValueError("window_size must be odd")

        self.window_size = window_size
        self.model: Optional[RandomForestClassifier] = None
        self.scaler: Optional[StandardScaler] = None
        self.is_trained = False

        self._model_dir = Path(__file__).parent / "weights"
        self._model_path = self._model_dir / "bcell_rf_model.joblib"
        self._scaler_path = self._model_dir / "bcell_scaler.joblib"

    def extract_features(
        self,
        sequence: str,
        position: int
    ) -> np.ndarray:
        """
        Extract all features for one position with window context.

        Args:
            sequence: Protein sequence.
            position: Position index (0-indexed).

        Returns:
            Feature vector of length 89.
        """
        return extract_window_features(sequence, position, self.window_size)

    def _extract_all_features(self, sequence: str) -> np.ndarray:
        """
        Extract features for all positions in a sequence.

        Args:
            sequence: Protein sequence.

        Returns:
            Feature matrix of shape (seq_len, 89).
        """
        seq_len = len(sequence)
        X = np.zeros((seq_len, 89), dtype=np.float32)
        for i in range(seq_len):
            X[i] = self.extract_features(sequence, i)
        return X

    def train(
        self,
        n_samples: int = 50000,
        n_estimators: int = 100,
        max_depth: int = 15,
        min_samples_leaf: int = 5,
        n_jobs: int = -1,
        seed: int = 42,
        verbose: bool = False
    ) -> Dict[str, Any]:
        """
        Train the model on synthetic data.

        Args:
            n_samples: Number of training samples to generate.
            n_estimators: Number of trees in the forest.
            max_depth: Maximum depth of each tree.
            min_samples_leaf: Minimum samples per leaf node.
            n_jobs: Number of parallel jobs (-1 for all cores).
            seed: Random seed for reproducibility.
            verbose: Whether to print training progress.

        Returns:
            Dictionary with training metrics.
        """
        if not SKLEARN_AVAILABLE:
            raise ImportError(
                "scikit-learn is required for training. "
                "Install with: pip install scikit-learn"
            )

        if verbose:
            print(f"Generating {n_samples} synthetic training samples...")

        X, y = generate_training_data(n_samples, self.window_size, seed)

        if verbose:
            print("Fitting scaler...")

        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(X)

        if verbose:
            print(f"Training Random Forest with {n_estimators} trees...")

        self.model = RandomForestClassifier(
            n_estimators=n_estimators,
            max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            n_jobs=n_jobs,
            random_state=seed,
            class_weight="balanced",
        )
        self.model.fit(X_scaled, y)

        self.is_trained = True

        # Calculate training metrics
        train_proba = self.model.predict_proba(X_scaled)[:, 1]
        train_pred = (train_proba >= 0.5).astype(int)
        accuracy = (train_pred == y).mean()

        # Feature importances
        importances = self.model.feature_importances_

        metrics = {
            "n_samples": n_samples,
            "n_features": X.shape[1],
            "n_estimators": n_estimators,
            "train_accuracy": float(accuracy),
            "n_positive": int(y.sum()),
            "n_negative": int(len(y) - y.sum()),
            "top_features": self._get_top_features(importances),
        }

        if verbose:
            print(f"Training complete. Accuracy: {accuracy:.3f}")

        return metrics

    def _get_top_features(
        self,
        importances: np.ndarray,
        top_k: int = 10
    ) -> List[Tuple[str, float]]:
        """Get the top-k most important features."""
        feature_names = []

        # One-hot features (0-19)
        for aa in AMINO_ACIDS:
            feature_names.append(f"onehot_{aa}")

        # Biophysical features (20-28)
        feature_names.extend([
            "hydrophobicity", "accessibility", "flexibility",
            "beta_turn", "antigenicity", "helix_propensity",
            "sheet_propensity", "charge", "volume"
        ])

        # Left context (29-57)
        for aa in AMINO_ACIDS:
            feature_names.append(f"left_onehot_{aa}")
        feature_names.extend([
            "left_hydro", "left_access", "left_flex", "left_turn",
            "left_antig", "left_helix", "left_sheet", "left_charge", "left_vol"
        ])

        # Right context (58-86)
        for aa in AMINO_ACIDS:
            feature_names.append(f"right_onehot_{aa}")
        feature_names.extend([
            "right_hydro", "right_access", "right_flex", "right_turn",
            "right_antig", "right_helix", "right_sheet", "right_charge", "right_vol"
        ])

        # Positional (87-88)
        feature_names.extend(["rel_position", "dist_from_termini"])

        # Get top features
        indices = np.argsort(importances)[::-1][:top_k]
        return [(feature_names[i], float(importances[i])) for i in indices]

    def save_model(self, path: Optional[str] = None) -> str:
        """
        Save the trained model to disk.

        Args:
            path: Directory to save model files (default: weights/).

        Returns:
            Path where model was saved.
        """
        if not self.is_trained:
            raise RuntimeError("Model must be trained before saving")

        if not JOBLIB_AVAILABLE:
            raise ImportError(
                "joblib is required to save models. "
                "Install with: pip install joblib"
            )

        if path:
            save_dir = Path(path)
        else:
            save_dir = self._model_dir

        save_dir.mkdir(parents=True, exist_ok=True)

        model_path = save_dir / "bcell_rf_model.joblib"
        scaler_path = save_dir / "bcell_scaler.joblib"

        joblib.dump(self.model, model_path)
        joblib.dump(self.scaler, scaler_path)

        # Also save metadata
        meta_path = save_dir / "bcell_meta.joblib"
        joblib.dump({
            "window_size": self.window_size,
            "n_features": 89,
        }, meta_path)

        return str(save_dir)

    def load_model(self, path: Optional[str] = None) -> bool:
        """
        Load a trained model from disk.

        Args:
            path: Directory containing model files (default: weights/).

        Returns:
            True if model loaded successfully.
        """
        if not JOBLIB_AVAILABLE:
            raise ImportError(
                "joblib is required to load models. "
                "Install with: pip install joblib"
            )

        if path:
            load_dir = Path(path)
        else:
            load_dir = self._model_dir

        model_path = load_dir / "bcell_rf_model.joblib"
        scaler_path = load_dir / "bcell_scaler.joblib"
        meta_path = load_dir / "bcell_meta.joblib"

        if not model_path.exists():
            return False

        self.model = joblib.load(model_path)
        self.scaler = joblib.load(scaler_path)

        if meta_path.exists():
            meta = joblib.load(meta_path)
            self.window_size = meta.get("window_size", 9)

        self.is_trained = True
        return True

    def _ensure_trained(self):
        """Ensure model is trained, loading or training if needed."""
        if self.is_trained:
            return

        # Try to load pre-trained model
        try:
            if self.load_model():
                return
        except Exception:
            pass

        # Train new model
        if SKLEARN_AVAILABLE:
            self.train(verbose=False)
        else:
            raise RuntimeError(
                "No trained model available and scikit-learn not installed. "
                "Install with: pip install scikit-learn"
            )

    def predict(self, sequence: str) -> np.ndarray:
        """
        Return per-residue epitope probability (0-1).

        Args:
            sequence: Protein sequence (amino acid one-letter codes).

        Returns:
            Array of probabilities for each position, length = len(sequence).
        """
        self._ensure_trained()

        # Normalize sequence
        sequence = sequence.upper().replace(" ", "").replace("\n", "")

        if len(sequence) == 0:
            return np.array([])

        # Extract features for all positions
        X = self._extract_all_features(sequence)

        # Scale features
        X_scaled = self.scaler.transform(X)

        # Get probabilities
        proba = self.model.predict_proba(X_scaled)[:, 1]

        return proba

    def predict_epitopes(
        self,
        sequence: str,
        threshold: float = 0.5,
        min_length: int = 6
    ) -> List[Dict[str, Any]]:
        """
        Return predicted epitope regions.

        Identifies contiguous regions where the prediction probability
        exceeds the threshold and meets minimum length requirements.

        Args:
            sequence: Protein sequence.
            threshold: Probability threshold for epitope calls (default 0.5).
            min_length: Minimum epitope length (default 6).

        Returns:
            List of epitope dictionaries with keys:
            - start: Start position (1-indexed)
            - end: End position (1-indexed)
            - sequence: Epitope sequence
            - score: Average probability in the region
            - length: Epitope length
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")

        if len(sequence) == 0:
            return []

        proba = self.predict(sequence)

        epitopes = []
        in_epitope = False
        start_idx = 0

        for i, p in enumerate(proba):
            if p >= threshold and not in_epitope:
                in_epitope = True
                start_idx = i
            elif p < threshold and in_epitope:
                in_epitope = False
                length = i - start_idx

                if length >= min_length:
                    epitope_seq = sequence[start_idx:i]
                    avg_score = float(np.mean(proba[start_idx:i]))
                    epitopes.append({
                        "start": start_idx + 1,  # 1-indexed
                        "end": i,  # 1-indexed (inclusive end is i, exclusive is i+1)
                        "sequence": epitope_seq,
                        "score": avg_score,
                        "length": length,
                    })

        # Handle epitope extending to end
        if in_epitope:
            length = len(sequence) - start_idx
            if length >= min_length:
                epitope_seq = sequence[start_idx:]
                avg_score = float(np.mean(proba[start_idx:]))
                epitopes.append({
                    "start": start_idx + 1,
                    "end": len(sequence),
                    "sequence": epitope_seq,
                    "score": avg_score,
                    "length": length,
                })

        return epitopes

    def get_feature_importance(self) -> Dict[str, float]:
        """
        Get feature importance scores from the trained model.

        Returns:
            Dictionary mapping feature names to importance scores.
        """
        self._ensure_trained()

        importances = self.model.feature_importances_
        top_features = self._get_top_features(importances, top_k=len(importances))

        return dict(top_features)


# ============================================================================
# Model Initialization and Pre-training
# ============================================================================

_DEFAULT_PREDICTOR: Optional[BcellLinearPredictor] = None


def get_predictor() -> BcellLinearPredictor:
    """
    Get the default B-cell linear epitope predictor.

    Returns a singleton instance that is trained on first use.

    Returns:
        Trained BcellLinearPredictor instance.
    """
    global _DEFAULT_PREDICTOR

    if _DEFAULT_PREDICTOR is None:
        _DEFAULT_PREDICTOR = BcellLinearPredictor()
        _DEFAULT_PREDICTOR._ensure_trained()

    return _DEFAULT_PREDICTOR


def predict_bcell_linear(
    sequence: str,
    threshold: float = 0.5,
    min_length: int = 6
) -> Dict[str, Any]:
    """
    Convenience function for B-cell linear epitope prediction.

    Args:
        sequence: Protein sequence.
        threshold: Probability threshold for epitope calls.
        min_length: Minimum epitope length.

    Returns:
        Dictionary with:
        - probabilities: Per-residue probabilities
        - epitopes: List of predicted epitope regions
        - sequence_length: Length of input sequence
    """
    predictor = get_predictor()
    proba = predictor.predict(sequence)
    epitopes = predictor.predict_epitopes(sequence, threshold, min_length)

    return {
        "probabilities": proba.tolist(),
        "epitopes": epitopes,
        "sequence_length": len(sequence),
        "threshold": threshold,
        "min_length": min_length,
    }


def train_and_save_model(
    output_dir: Optional[str] = None,
    n_samples: int = 50000,
    verbose: bool = True
) -> str:
    """
    Train a new model and save it to disk.

    This function is used to generate pre-trained model weights.

    Args:
        output_dir: Directory to save model files.
        n_samples: Number of training samples.
        verbose: Whether to print progress.

    Returns:
        Path where model was saved.
    """
    predictor = BcellLinearPredictor()
    metrics = predictor.train(n_samples=n_samples, verbose=verbose)

    if verbose:
        print(f"\nTraining metrics:")
        print(f"  Samples: {metrics['n_samples']}")
        print(f"  Features: {metrics['n_features']}")
        print(f"  Accuracy: {metrics['train_accuracy']:.3f}")
        print(f"\nTop 10 features:")
        for name, imp in metrics['top_features']:
            print(f"  {name}: {imp:.4f}")

    save_path = predictor.save_model(output_dir)

    if verbose:
        print(f"\nModel saved to: {save_path}")

    return save_path


# Pre-train and save model on module load (if weights don't exist)
def _init_pretrained_model():
    """Initialize pre-trained model weights if they don't exist."""
    weights_dir = Path(__file__).parent / "weights"
    model_path = weights_dir / "bcell_rf_model.joblib"

    if not model_path.exists() and SKLEARN_AVAILABLE and JOBLIB_AVAILABLE:
        try:
            train_and_save_model(str(weights_dir), n_samples=50000, verbose=False)
        except Exception:
            # Silently fail - model will be trained on first use
            pass


# Initialize on import (if dependencies available)
if SKLEARN_AVAILABLE and JOBLIB_AVAILABLE:
    _init_pretrained_model()
