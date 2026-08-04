"""
MHC Class I binding prediction using neural network.

This module provides a neural network-based predictor for MHC Class I
peptide binding affinity. It uses a combination of one-hot encoding
and BLOSUM62-based features to represent peptides.

The model is trained on synthetic data based on known anchor preferences
for common HLA alleles, with pre-trained weights bundled for offline use.

Example:
    >>> from isis.models.mhc1_model import MHC1Predictor
    >>> predictor = MHC1Predictor(allele='HLA-A*02:01')
    >>> ic50_values = predictor.predict(['YLQPRTFLL', 'GILGFVFTL'])
    >>> percentiles = predictor.predict_percentile(['YLQPRTFLL', 'GILGFVFTL'])
"""
from __future__ import annotations

import hashlib
import os
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

try:
    from sklearn.neural_network import MLPRegressor
    from sklearn.preprocessing import StandardScaler
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


# =============================================================================
# Constants
# =============================================================================

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_IDX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
N_AMINO_ACIDS = len(AMINO_ACIDS)

# Supported peptide lengths
MIN_PEPTIDE_LENGTH = 8
MAX_PEPTIDE_LENGTH = 11
CANONICAL_LENGTH = 9  # Pad/truncate to this for model input

# Model hyperparameters
HIDDEN_LAYER_SIZES = (128, 64, 32)
MAX_ITER = 500
RANDOM_STATE = 42

# IC50 bounds for predictions (nM)
MIN_IC50 = 1.0
MAX_IC50 = 50000.0


# =============================================================================
# BLOSUM62 Substitution Matrix
# =============================================================================

# BLOSUM62 matrix for the 20 standard amino acids
# Order: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
BLOSUM62 = {
    'A': {'A': 4, 'C': 0, 'D':-2, 'E':-1, 'F':-2, 'G': 0, 'H':-2, 'I':-1, 'K':-1, 'L':-1,
          'M':-1, 'N':-2, 'P':-1, 'Q':-1, 'R':-1, 'S': 1, 'T': 0, 'V': 0, 'W':-3, 'Y':-2},
    'C': {'A': 0, 'C': 9, 'D':-3, 'E':-4, 'F':-2, 'G':-3, 'H':-3, 'I':-1, 'K':-3, 'L':-1,
          'M':-1, 'N':-3, 'P':-3, 'Q':-3, 'R':-3, 'S':-1, 'T':-1, 'V':-1, 'W':-2, 'Y':-2},
    'D': {'A':-2, 'C':-3, 'D': 6, 'E': 2, 'F':-3, 'G':-1, 'H':-1, 'I':-3, 'K':-1, 'L':-4,
          'M':-3, 'N': 1, 'P':-1, 'Q': 0, 'R':-2, 'S': 0, 'T':-1, 'V':-3, 'W':-4, 'Y':-3},
    'E': {'A':-1, 'C':-4, 'D': 2, 'E': 5, 'F':-3, 'G':-2, 'H': 0, 'I':-3, 'K': 1, 'L':-3,
          'M':-2, 'N': 0, 'P':-1, 'Q': 2, 'R': 0, 'S': 0, 'T':-1, 'V':-2, 'W':-3, 'Y':-2},
    'F': {'A':-2, 'C':-2, 'D':-3, 'E':-3, 'F': 6, 'G':-3, 'H':-1, 'I': 0, 'K':-3, 'L': 0,
          'M': 0, 'N':-3, 'P':-4, 'Q':-3, 'R':-3, 'S':-2, 'T':-2, 'V':-1, 'W': 1, 'Y': 3},
    'G': {'A': 0, 'C':-3, 'D':-1, 'E':-2, 'F':-3, 'G': 6, 'H':-2, 'I':-4, 'K':-2, 'L':-4,
          'M':-3, 'N': 0, 'P':-2, 'Q':-2, 'R':-2, 'S': 0, 'T':-2, 'V':-3, 'W':-2, 'Y':-3},
    'H': {'A':-2, 'C':-3, 'D':-1, 'E': 0, 'F':-1, 'G':-2, 'H': 8, 'I':-3, 'K':-1, 'L':-3,
          'M':-2, 'N': 1, 'P':-2, 'Q': 0, 'R': 0, 'S':-1, 'T':-2, 'V':-3, 'W':-2, 'Y': 2},
    'I': {'A':-1, 'C':-1, 'D':-3, 'E':-3, 'F': 0, 'G':-4, 'H':-3, 'I': 4, 'K':-3, 'L': 2,
          'M': 1, 'N':-3, 'P':-3, 'Q':-3, 'R':-3, 'S':-2, 'T':-1, 'V': 3, 'W':-3, 'Y':-1},
    'K': {'A':-1, 'C':-3, 'D':-1, 'E': 1, 'F':-3, 'G':-2, 'H':-1, 'I':-3, 'K': 5, 'L':-2,
          'M':-1, 'N': 0, 'P':-1, 'Q': 1, 'R': 2, 'S': 0, 'T':-1, 'V':-2, 'W':-3, 'Y':-2},
    'L': {'A':-1, 'C':-1, 'D':-4, 'E':-3, 'F': 0, 'G':-4, 'H':-3, 'I': 2, 'K':-2, 'L': 4,
          'M': 2, 'N':-3, 'P':-3, 'Q':-2, 'R':-2, 'S':-2, 'T':-1, 'V': 1, 'W':-2, 'Y':-1},
    'M': {'A':-1, 'C':-1, 'D':-3, 'E':-2, 'F': 0, 'G':-3, 'H':-2, 'I': 1, 'K':-1, 'L': 2,
          'M': 5, 'N':-2, 'P':-2, 'Q': 0, 'R':-1, 'S':-1, 'T':-1, 'V': 1, 'W':-1, 'Y':-1},
    'N': {'A':-2, 'C':-3, 'D': 1, 'E': 0, 'F':-3, 'G': 0, 'H': 1, 'I':-3, 'K': 0, 'L':-3,
          'M':-2, 'N': 6, 'P':-2, 'Q': 0, 'R': 0, 'S': 1, 'T': 0, 'V':-3, 'W':-4, 'Y':-2},
    'P': {'A':-1, 'C':-3, 'D':-1, 'E':-1, 'F':-4, 'G':-2, 'H':-2, 'I':-3, 'K':-1, 'L':-3,
          'M':-2, 'N':-2, 'P': 7, 'Q':-1, 'R':-2, 'S':-1, 'T':-1, 'V':-2, 'W':-4, 'Y':-3},
    'Q': {'A':-1, 'C':-3, 'D': 0, 'E': 2, 'F':-3, 'G':-2, 'H': 0, 'I':-3, 'K': 1, 'L':-2,
          'M': 0, 'N': 0, 'P':-1, 'Q': 5, 'R': 1, 'S': 0, 'T':-1, 'V':-2, 'W':-2, 'Y':-1},
    'R': {'A':-1, 'C':-3, 'D':-2, 'E': 0, 'F':-3, 'G':-2, 'H': 0, 'I':-3, 'K': 2, 'L':-2,
          'M':-1, 'N': 0, 'P':-2, 'Q': 1, 'R': 5, 'S':-1, 'T':-1, 'V':-3, 'W':-3, 'Y':-2},
    'S': {'A': 1, 'C':-1, 'D': 0, 'E': 0, 'F':-2, 'G': 0, 'H':-1, 'I':-2, 'K': 0, 'L':-2,
          'M':-1, 'N': 1, 'P':-1, 'Q': 0, 'R':-1, 'S': 4, 'T': 1, 'V':-2, 'W':-3, 'Y':-2},
    'T': {'A': 0, 'C':-1, 'D':-1, 'E':-1, 'F':-2, 'G':-2, 'H':-2, 'I':-1, 'K':-1, 'L':-1,
          'M':-1, 'N': 0, 'P':-1, 'Q':-1, 'R':-1, 'S': 1, 'T': 5, 'V': 0, 'W':-2, 'Y':-2},
    'V': {'A': 0, 'C':-1, 'D':-3, 'E':-2, 'F':-1, 'G':-3, 'H':-3, 'I': 3, 'K':-2, 'L': 1,
          'M': 1, 'N':-3, 'P':-2, 'Q':-2, 'R':-3, 'S':-2, 'T': 0, 'V': 4, 'W':-3, 'Y':-1},
    'W': {'A':-3, 'C':-2, 'D':-4, 'E':-3, 'F': 1, 'G':-2, 'H':-2, 'I':-3, 'K':-3, 'L':-2,
          'M':-1, 'N':-4, 'P':-4, 'Q':-2, 'R':-3, 'S':-3, 'T':-2, 'V':-3, 'W':11, 'Y': 2},
    'Y': {'A':-2, 'C':-2, 'D':-3, 'E':-2, 'F': 3, 'G':-3, 'H': 2, 'I':-1, 'K':-2, 'L':-1,
          'M':-1, 'N':-2, 'P':-3, 'Q':-1, 'R':-2, 'S':-2, 'T':-2, 'V':-1, 'W': 2, 'Y': 7},
}

# Convert BLOSUM62 to numpy array for faster encoding
BLOSUM62_MATRIX = np.zeros((N_AMINO_ACIDS, N_AMINO_ACIDS), dtype=np.float32)
for i, aa1 in enumerate(AMINO_ACIDS):
    for j, aa2 in enumerate(AMINO_ACIDS):
        BLOSUM62_MATRIX[i, j] = BLOSUM62[aa1][aa2]


# =============================================================================
# HLA Allele Anchor Preferences
# =============================================================================

# Anchor residue preferences for different HLA alleles
# Position 2 and C-terminus (position 9 for 9-mers) are primary anchors
# Values are log-odds scores (higher = more preferred)
ALLELE_ANCHORS = {
    'HLA-A*02:01': {
        'position_2': {
            'L': 2.0, 'M': 1.8, 'I': 1.2, 'V': 0.8, 'A': 0.3, 'T': 0.1,
            'Q': -0.5, 'F': -0.3, 'S': -0.4, 'N': -0.8, 'G': -1.0,
            'D': -1.5, 'E': -1.5, 'K': -1.2, 'R': -1.0, 'P': -1.5,
            'H': -0.6, 'W': -0.5, 'Y': -0.3, 'C': -0.2
        },
        'position_c': {
            'V': 2.2, 'L': 2.0, 'I': 1.5, 'A': 0.5, 'T': 0.3, 'M': 1.0,
            'F': 0.2, 'Y': -0.2, 'W': -0.5, 'S': -0.3, 'N': -0.8,
            'G': -1.0, 'D': -1.5, 'E': -1.5, 'K': -1.0, 'R': -0.8,
            'P': -1.5, 'H': -0.5, 'Q': -0.5, 'C': -0.2
        },
        'auxiliary': {
            # Positions 1, 3, 6 have some preferences
            1: {'G': 0.3, 'A': 0.2, 'S': 0.2, 'T': 0.1},
            3: {'L': 0.3, 'V': 0.2, 'I': 0.2, 'A': 0.1},
            6: {'V': 0.2, 'I': 0.2, 'L': 0.2, 'T': 0.1},
        }
    },
    'HLA-A*01:01': {
        'position_2': {
            'T': 2.0, 'S': 1.8, 'M': 0.5, 'A': 0.3, 'L': -0.3, 'V': -0.4,
            'I': -0.5, 'Q': -0.3, 'N': -0.2, 'G': -0.8, 'D': -1.0,
            'E': -1.2, 'K': -0.8, 'R': -0.8, 'P': -1.5, 'H': -0.4,
            'F': -0.5, 'W': -0.6, 'Y': -0.4, 'C': -0.3
        },
        'position_c': {
            'Y': 2.5, 'F': 1.2, 'W': 0.8, 'L': -0.5, 'V': -0.8, 'I': -0.6,
            'A': -0.5, 'T': -0.3, 'M': -0.4, 'S': -0.5, 'N': -0.6,
            'G': -1.0, 'D': -1.5, 'E': -1.5, 'K': -0.8, 'R': -0.6,
            'P': -1.5, 'H': -0.3, 'Q': -0.5, 'C': -0.4
        },
        'auxiliary': {
            3: {'D': 0.5, 'E': 0.4, 'A': 0.2, 'S': 0.1},
        }
    },
    'HLA-B*07:02': {
        'position_2': {
            'P': 2.5, 'A': 0.5, 'S': 0.3, 'T': 0.2, 'L': -0.5, 'V': -0.6,
            'I': -0.5, 'M': -0.3, 'Q': -0.4, 'N': -0.5, 'G': -0.6,
            'D': -1.0, 'E': -1.2, 'K': -0.8, 'R': -0.8, 'H': -0.5,
            'F': -0.8, 'W': -0.8, 'Y': -0.6, 'C': -0.4
        },
        'position_c': {
            'L': 2.0, 'I': 1.5, 'V': 1.2, 'M': 1.0, 'F': 0.6, 'A': 0.4,
            'Y': 0.2, 'W': 0.0, 'T': -0.2, 'S': -0.4, 'N': -0.6,
            'G': -1.0, 'D': -1.2, 'E': -1.2, 'K': -0.5, 'R': -0.4,
            'P': -1.5, 'H': -0.4, 'Q': -0.5, 'C': -0.3
        },
        'auxiliary': {}
    },
    'HLA-A*03:01': {
        'position_2': {
            'L': 1.5, 'V': 1.2, 'M': 1.0, 'I': 0.8, 'A': 0.5, 'T': 0.3,
            'S': 0.2, 'Q': -0.3, 'N': -0.4, 'G': -0.8, 'D': -1.2,
            'E': -1.0, 'K': -0.5, 'R': -0.5, 'P': -1.5, 'H': -0.4,
            'F': -0.2, 'W': -0.4, 'Y': -0.3, 'C': -0.2
        },
        'position_c': {
            'K': 2.0, 'R': 1.8, 'Y': 0.5, 'F': 0.3, 'L': -0.3, 'V': -0.5,
            'I': -0.4, 'A': -0.5, 'T': -0.4, 'M': -0.3, 'S': -0.5,
            'N': -0.4, 'G': -1.0, 'D': -1.2, 'E': -1.0, 'P': -1.5,
            'H': 0.3, 'Q': -0.3, 'W': -0.3, 'C': -0.4
        },
        'auxiliary': {}
    },
    'HLA-B*08:01': {
        'position_2': {
            'R': 2.0, 'K': 1.8, 'H': 1.0, 'Q': 0.5, 'L': -0.3, 'V': -0.4,
            'I': -0.4, 'M': -0.3, 'A': -0.2, 'T': -0.3, 'S': -0.4,
            'N': -0.3, 'G': -0.8, 'D': -1.0, 'E': -0.8, 'P': -1.5,
            'F': -0.5, 'W': -0.6, 'Y': -0.4, 'C': -0.4
        },
        'position_c': {
            'L': 2.0, 'K': 1.2, 'R': 1.0, 'I': 0.8, 'V': 0.6, 'M': 0.5,
            'F': 0.3, 'Y': 0.2, 'A': 0.0, 'T': -0.3, 'S': -0.5,
            'N': -0.6, 'G': -1.0, 'D': -1.2, 'E': -1.0, 'P': -1.5,
            'H': -0.2, 'Q': -0.4, 'W': -0.2, 'C': -0.4
        },
        'auxiliary': {}
    },
}


# =============================================================================
# Feature Encoding Functions
# =============================================================================

def one_hot_encode_aa(aa: str) -> np.ndarray:
    """
    One-hot encode a single amino acid.

    Args:
        aa: Single amino acid character.

    Returns:
        One-hot vector of length 20.
    """
    vec = np.zeros(N_AMINO_ACIDS, dtype=np.float32)
    if aa in AA_TO_IDX:
        vec[AA_TO_IDX[aa]] = 1.0
    return vec


def blosum62_encode_aa(aa: str) -> np.ndarray:
    """
    Encode amino acid using its BLOSUM62 row.

    This provides a physicochemical similarity representation
    rather than pure identity.

    Args:
        aa: Single amino acid character.

    Returns:
        BLOSUM62 row vector of length 20.
    """
    if aa in AA_TO_IDX:
        return BLOSUM62_MATRIX[AA_TO_IDX[aa]].copy()
    return np.zeros(N_AMINO_ACIDS, dtype=np.float32)


def encode_peptide(
    peptide: str,
    max_length: int = CANONICAL_LENGTH,
    use_blosum: bool = True
) -> np.ndarray:
    """
    Encode a peptide as feature vector combining one-hot and BLOSUM62.

    For peptides shorter than max_length, padding is added at the
    center to preserve anchor positions (2 and C-terminus).

    For peptides longer than max_length, the middle residues are
    compressed while preserving anchors.

    Args:
        peptide: Amino acid sequence (8-11 residues).
        max_length: Target length for encoding (default 9).
        use_blosum: Include BLOSUM62 features in addition to one-hot.

    Returns:
        Feature vector of shape (max_length * features_per_position,).
    """
    peptide = peptide.upper()
    pep_len = len(peptide)

    # Normalize length to canonical length
    if pep_len < max_length:
        # Pad in the middle to preserve N-term (pos 1-2) and C-term
        n_pad = max_length - pep_len
        mid = pep_len // 2
        # Insert padding (represented as zero vectors) in middle
        normalized = peptide[:mid + 1] + 'X' * n_pad + peptide[mid + 1:]
    elif pep_len > max_length:
        # Keep anchor positions, compress middle
        # Keep positions 1, 2, and last 2
        keep_n = 3  # First 3 positions
        keep_c = 2  # Last 2 positions
        middle_len = max_length - keep_n - keep_c
        # Sample middle positions
        middle_start = keep_n
        middle_end = pep_len - keep_c
        middle_indices = np.linspace(middle_start, middle_end - 1, middle_len, dtype=int)
        normalized = peptide[:keep_n] + ''.join(peptide[i] for i in middle_indices) + peptide[-keep_c:]
    else:
        normalized = peptide

    # Encode each position
    features_per_pos = N_AMINO_ACIDS * (2 if use_blosum else 1)
    features = np.zeros(max_length * features_per_pos, dtype=np.float32)

    for i, aa in enumerate(normalized):
        start_idx = i * features_per_pos

        if aa == 'X':
            # Padding - leave as zeros
            continue

        # One-hot encoding
        onehot = one_hot_encode_aa(aa)
        features[start_idx:start_idx + N_AMINO_ACIDS] = onehot

        # BLOSUM62 encoding
        if use_blosum:
            blosum = blosum62_encode_aa(aa)
            # Normalize BLOSUM scores to [0, 1] range
            blosum_norm = (blosum + 4) / 15  # BLOSUM62 range is roughly -4 to 11
            features[start_idx + N_AMINO_ACIDS:start_idx + 2 * N_AMINO_ACIDS] = blosum_norm

    return features


def encode_peptides(
    peptides: List[str],
    max_length: int = CANONICAL_LENGTH,
    use_blosum: bool = True
) -> np.ndarray:
    """
    Encode multiple peptides as feature matrix.

    Args:
        peptides: List of peptide sequences.
        max_length: Target length for encoding.
        use_blosum: Include BLOSUM62 features.

    Returns:
        Feature matrix of shape (n_peptides, n_features).
    """
    return np.array([encode_peptide(p, max_length, use_blosum) for p in peptides])


# =============================================================================
# Synthetic Training Data Generation
# =============================================================================

def generate_random_peptide(length: int = 9, rng: np.random.Generator = None) -> str:
    """Generate a random peptide of given length."""
    if rng is None:
        rng = np.random.default_rng()
    return ''.join(rng.choice(list(AMINO_ACIDS), size=length))


def compute_anchor_score(peptide: str, allele: str) -> float:
    """
    Compute anchor-based binding score for a peptide.

    Args:
        peptide: Amino acid sequence.
        allele: HLA allele identifier.

    Returns:
        Log-odds score based on anchor residues.
    """
    if allele not in ALLELE_ANCHORS:
        allele = 'HLA-A*02:01'  # Default

    anchors = ALLELE_ANCHORS[allele]
    score = 0.0

    # Position 2 anchor (strong)
    if len(peptide) >= 2:
        aa2 = peptide[1]
        score += anchors['position_2'].get(aa2, -0.5) * 2.0

    # C-terminal anchor (strong)
    aa_c = peptide[-1]
    score += anchors['position_c'].get(aa_c, -0.5) * 2.0

    # Auxiliary positions (weak)
    for pos, prefs in anchors.get('auxiliary', {}).items():
        if pos <= len(peptide):
            aa = peptide[pos - 1]
            score += prefs.get(aa, 0.0) * 0.5

    return score


def score_to_ic50(score: float, baseline: float = 1000.0) -> float:
    """
    Convert anchor score to IC50 value.

    Higher scores (better binding) -> lower IC50.
    Uses exponential transformation.

    Args:
        score: Anchor-based binding score.
        baseline: Baseline IC50 for neutral peptides.

    Returns:
        Predicted IC50 in nanomolar.
    """
    # Score of 0 -> baseline IC50
    # Score of +4 -> ~50 nM (strong binder)
    # Score of -4 -> ~20000 nM (non-binder)
    ic50 = baseline * np.exp(-score * 0.5)
    return np.clip(ic50, MIN_IC50, MAX_IC50)


def ic50_to_log_ic50(ic50: float) -> float:
    """Convert IC50 to log-transformed value for training."""
    return np.log10(np.clip(ic50, MIN_IC50, MAX_IC50))


def log_ic50_to_ic50(log_ic50: float) -> float:
    """Convert log-transformed IC50 back to nM."""
    return np.clip(10 ** log_ic50, MIN_IC50, MAX_IC50)


def generate_training_data(
    allele: str,
    n_samples: int = 10000,
    peptide_lengths: List[int] = None,
    seed: int = RANDOM_STATE
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """
    Generate synthetic training data based on anchor preferences.

    Creates peptides with varying binding affinities by:
    1. Generating random peptides
    2. For some, specifically engineering strong anchor residues
    3. Computing IC50 based on anchor rules + noise

    Args:
        allele: HLA allele identifier.
        n_samples: Number of training samples.
        peptide_lengths: Allowed peptide lengths (default [8,9,10,11]).
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (features, targets, peptides):
        - features: Encoded peptide features, shape (n_samples, n_features)
        - targets: Log10(IC50) values, shape (n_samples,)
        - peptides: List of peptide sequences
    """
    if peptide_lengths is None:
        peptide_lengths = [8, 9, 10, 11]

    rng = np.random.default_rng(seed)
    peptides = []
    ic50_values = []

    if allele not in ALLELE_ANCHORS:
        allele = 'HLA-A*02:01'
    anchors = ALLELE_ANCHORS[allele]

    # Generate three categories of peptides
    n_strong = n_samples // 4      # Strong binders (engineered)
    n_medium = n_samples // 4      # Medium binders
    n_random = n_samples - n_strong - n_medium  # Random (mostly weak)

    # Strong binders - engineer optimal anchors
    strong_p2 = [aa for aa, score in anchors['position_2'].items() if score > 1.0]
    strong_pc = [aa for aa, score in anchors['position_c'].items() if score > 1.0]

    if not strong_p2:
        strong_p2 = ['L', 'M', 'I']
    if not strong_pc:
        strong_pc = ['V', 'L', 'I']

    for _ in range(n_strong):
        length = rng.choice(peptide_lengths)
        peptide = list(generate_random_peptide(length, rng))
        # Insert optimal anchors
        peptide[1] = rng.choice(strong_p2)
        peptide[-1] = rng.choice(strong_pc)
        peptide = ''.join(peptide)

        score = compute_anchor_score(peptide, allele)
        # Add noise
        score += rng.normal(0, 0.3)
        ic50 = score_to_ic50(score, baseline=800)

        peptides.append(peptide)
        ic50_values.append(ic50)

    # Medium binders - one optimal anchor
    for _ in range(n_medium):
        length = rng.choice(peptide_lengths)
        peptide = list(generate_random_peptide(length, rng))
        # Insert one optimal anchor
        if rng.random() > 0.5:
            peptide[1] = rng.choice(strong_p2)
        else:
            peptide[-1] = rng.choice(strong_pc)
        peptide = ''.join(peptide)

        score = compute_anchor_score(peptide, allele)
        score += rng.normal(0, 0.4)
        ic50 = score_to_ic50(score, baseline=1000)

        peptides.append(peptide)
        ic50_values.append(ic50)

    # Random peptides - mostly weak binders
    for _ in range(n_random):
        length = rng.choice(peptide_lengths)
        peptide = generate_random_peptide(length, rng)

        score = compute_anchor_score(peptide, allele)
        score += rng.normal(0, 0.5)
        ic50 = score_to_ic50(score, baseline=1200)

        peptides.append(peptide)
        ic50_values.append(ic50)

    # Shuffle all data
    indices = rng.permutation(len(peptides))
    peptides = [peptides[i] for i in indices]
    ic50_values = [ic50_values[i] for i in indices]

    # Encode features
    features = encode_peptides(peptides, max_length=CANONICAL_LENGTH, use_blosum=True)
    targets = np.array([ic50_to_log_ic50(ic50) for ic50 in ic50_values], dtype=np.float32)

    return features, targets, peptides


# =============================================================================
# MHC1Predictor Class
# =============================================================================

class MHC1Predictor:
    """
    Neural network-based MHC Class I binding predictor.

    Uses an MLPRegressor with combined one-hot and BLOSUM62 features
    to predict peptide-MHC binding affinity. Pre-trained weights are
    bundled for common HLA alleles.

    Attributes:
        allele: HLA allele identifier (e.g., 'HLA-A*02:01').
        model: Trained MLPRegressor or None if not loaded.
        scaler: Feature scaler for normalization.

    Example:
        >>> predictor = MHC1Predictor('HLA-A*02:01')
        >>> ic50 = predictor.predict(['YLQPRTFLL', 'GILGFVFTL'])
        >>> print(f"IC50 values: {ic50}")
    """

    def __init__(
        self,
        allele: str = 'HLA-A*02:01',
        model_dir: Optional[str] = None
    ):
        """
        Initialize the MHC1 predictor.

        Args:
            allele: HLA allele identifier. Supported alleles:
                - HLA-A*02:01 (most common, default)
                - HLA-A*01:01
                - HLA-A*03:01
                - HLA-B*07:02
                - HLA-B*08:01
            model_dir: Directory containing pre-trained weights.
                If None, uses bundled weights in package directory.
        """
        self.allele = self._normalize_allele(allele)
        self.model: Optional[MLPRegressor] = None
        self.scaler: Optional[StandardScaler] = None
        self._reference_ic50s: Optional[np.ndarray] = None

        # Try to load pre-trained model
        if model_dir:
            self._model_dir = Path(model_dir)
        else:
            self._model_dir = Path(__file__).parent / 'weights'

        self._load_or_train_model()

    def _normalize_allele(self, allele: str) -> str:
        """Normalize allele name to standard format."""
        allele = allele.upper().replace(' ', '')

        # Handle common variants
        allele = allele.replace('HLA_', 'HLA-')
        allele = allele.replace('A02:01', 'A*02:01')
        allele = allele.replace('A0201', 'A*02:01')
        allele = allele.replace('B07:02', 'B*07:02')
        allele = allele.replace('B0702', 'B*07:02')

        # Ensure HLA- prefix
        if not allele.startswith('HLA-'):
            if allele.startswith('A*') or allele.startswith('B*'):
                allele = 'HLA-' + allele

        return allele

    def _get_model_path(self) -> Path:
        """Get path for model weights file."""
        safe_name = self.allele.replace('*', '_').replace(':', '_')
        return self._model_dir / f'mhc1_{safe_name}.pkl'

    def _load_or_train_model(self) -> None:
        """Load pre-trained model or train a new one."""
        model_path = self._get_model_path()

        if model_path.exists():
            self._load_model(model_path)
        else:
            # Train new model with synthetic data
            self._train_synthetic_model()
            # Save for future use
            self._save_model(model_path)

    def _load_model(self, path: Path) -> None:
        """Load model from pickle file."""
        with open(path, 'rb') as f:
            data = pickle.load(f)

        self.model = data['model']
        self.scaler = data['scaler']
        self._reference_ic50s = data.get('reference_ic50s')

    def _save_model(self, path: Path) -> None:
        """Save model to pickle file."""
        path.parent.mkdir(parents=True, exist_ok=True)

        data = {
            'model': self.model,
            'scaler': self.scaler,
            'allele': self.allele,
            'reference_ic50s': self._reference_ic50s,
        }

        with open(path, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

    def _train_synthetic_model(self, n_samples: int = 10000) -> None:
        """Train model on synthetic data."""
        if not HAS_SKLEARN:
            raise ImportError(
                "scikit-learn is required for training. "
                "Install with: pip install scikit-learn"
            )

        # Generate training data
        X, y, peptides = generate_training_data(
            self.allele,
            n_samples=n_samples,
            seed=RANDOM_STATE
        )

        # Scale features
        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(X)

        # Train MLPRegressor
        self.model = MLPRegressor(
            hidden_layer_sizes=HIDDEN_LAYER_SIZES,
            activation='relu',
            solver='adam',
            alpha=0.001,  # L2 regularization
            batch_size='auto',
            learning_rate='adaptive',
            learning_rate_init=0.001,
            max_iter=MAX_ITER,
            shuffle=True,
            random_state=RANDOM_STATE,
            early_stopping=True,
            validation_fraction=0.1,
            n_iter_no_change=20,
            verbose=False
        )

        self.model.fit(X_scaled, y)

        # Generate reference distribution for percentile calculation
        n_ref = 5000
        ref_peptides = [generate_random_peptide(9) for _ in range(n_ref)]
        self._reference_ic50s = self.predict(ref_peptides)
        self._reference_ic50s.sort()

    def encode_peptide(self, peptide: str) -> np.ndarray:
        """
        Encode a peptide to feature vector.

        Combines one-hot encoding with BLOSUM62-based features.

        Args:
            peptide: Amino acid sequence (8-11 residues).

        Returns:
            Feature vector of shape (n_features,).

        Raises:
            ValueError: If peptide length is outside supported range.
        """
        peptide = peptide.upper().strip()

        if len(peptide) < MIN_PEPTIDE_LENGTH or len(peptide) > MAX_PEPTIDE_LENGTH:
            raise ValueError(
                f"Peptide length must be {MIN_PEPTIDE_LENGTH}-{MAX_PEPTIDE_LENGTH}, "
                f"got {len(peptide)}"
            )

        # Check for invalid amino acids
        for aa in peptide:
            if aa not in AMINO_ACIDS:
                raise ValueError(f"Invalid amino acid: {aa}")

        return encode_peptide(peptide, CANONICAL_LENGTH, use_blosum=True)

    def predict(self, peptides: Union[str, List[str]]) -> np.ndarray:
        """
        Predict binding affinity (IC50 in nM) for peptides.

        Args:
            peptides: Single peptide string or list of peptides.

        Returns:
            Array of IC50 predictions in nanomolar.
            Lower values indicate stronger binding.
            Typical interpretation:
            - IC50 < 50 nM: Strong binder
            - IC50 < 500 nM: Moderate binder
            - IC50 < 5000 nM: Weak binder
            - IC50 >= 5000 nM: Non-binder
        """
        if isinstance(peptides, str):
            peptides = [peptides]

        if not peptides:
            return np.array([], dtype=np.float32)

        # Encode all peptides
        X = encode_peptides(peptides, CANONICAL_LENGTH, use_blosum=True)

        # Scale features
        if self.scaler is not None:
            X = self.scaler.transform(X)

        # Predict log(IC50)
        if self.model is not None:
            log_ic50 = self.model.predict(X)
        else:
            # Fallback to simple anchor-based prediction
            log_ic50 = np.array([
                ic50_to_log_ic50(score_to_ic50(compute_anchor_score(p, self.allele)))
                for p in peptides
            ])

        # Convert to IC50
        ic50 = np.array([log_ic50_to_ic50(x) for x in log_ic50])

        return ic50

    def predict_percentile(self, peptides: Union[str, List[str]]) -> np.ndarray:
        """
        Predict binding as percentile rank.

        Percentile rank indicates what fraction of random peptides
        have better (lower) IC50 values. Lower percentile = stronger binder.

        Args:
            peptides: Single peptide string or list of peptides.

        Returns:
            Array of percentile ranks (0-100).
            Typical interpretation:
            - Percentile < 0.5%: Strong binder
            - Percentile < 2%: Moderate binder
            - Percentile < 10%: Weak binder
            - Percentile >= 10%: Non-binder
        """
        ic50_values = self.predict(peptides)

        if self._reference_ic50s is None or len(self._reference_ic50s) == 0:
            # Fallback: estimate percentile from IC50 directly
            percentiles = np.zeros_like(ic50_values)
            for i, ic50 in enumerate(ic50_values):
                if ic50 < 50:
                    percentiles[i] = 0.1
                elif ic50 < 500:
                    percentiles[i] = 0.5 + (ic50 - 50) / 450 * 1.5
                elif ic50 < 5000:
                    percentiles[i] = 2.0 + (ic50 - 500) / 4500 * 8.0
                else:
                    percentiles[i] = 10.0 + min(90.0, (ic50 - 5000) / 45000 * 90.0)
            return percentiles

        # Use reference distribution
        percentiles = np.searchsorted(self._reference_ic50s, ic50_values) / len(self._reference_ic50s) * 100
        return percentiles.astype(np.float32)

    def predict_binding_level(self, peptides: Union[str, List[str]]) -> List[str]:
        """
        Classify peptides into binding level categories.

        Args:
            peptides: Single peptide string or list of peptides.

        Returns:
            List of binding level strings:
            - 'Strong': IC50 < 50 nM
            - 'Moderate': 50 <= IC50 < 500 nM
            - 'Weak': 500 <= IC50 < 5000 nM
            - 'Non-binder': IC50 >= 5000 nM
        """
        ic50_values = self.predict(peptides)

        levels = []
        for ic50 in ic50_values:
            if ic50 < 50:
                levels.append('Strong')
            elif ic50 < 500:
                levels.append('Moderate')
            elif ic50 < 5000:
                levels.append('Weak')
            else:
                levels.append('Non-binder')

        return levels

    @classmethod
    def available_alleles(cls) -> List[str]:
        """Return list of supported HLA alleles."""
        return list(ALLELE_ANCHORS.keys())

    def get_anchor_preferences(self) -> Dict[str, Any]:
        """
        Get anchor residue preferences for the current allele.

        Returns:
            Dictionary with position_2, position_c, and auxiliary preferences.
        """
        if self.allele in ALLELE_ANCHORS:
            return ALLELE_ANCHORS[self.allele].copy()
        return ALLELE_ANCHORS['HLA-A*02:01'].copy()


# =============================================================================
# Training Script
# =============================================================================

def train_from_iedb(
    data_file: str,
    allele: str,
    output_dir: Optional[str] = None,
    test_fraction: float = 0.1
) -> Dict[str, float]:
    """
    Train model from IEDB binding affinity data.

    IEDB data can be downloaded from:
    https://www.iedb.org/database_export_v3.php

    Expected file format (CSV or TSV):
    - Column 'Description' or 'Epitope': Peptide sequence
    - Column 'Quantitative measurement': IC50 value
    - Column 'Allele Name' or 'MHC Allele': HLA allele

    Args:
        data_file: Path to IEDB export file.
        allele: HLA allele to train for.
        output_dir: Directory to save trained model.
        test_fraction: Fraction of data for testing.

    Returns:
        Dictionary with training metrics:
        - 'mse': Mean squared error on test set
        - 'r2': R-squared correlation
        - 'n_train': Number of training samples
        - 'n_test': Number of test samples
    """
    if not HAS_SKLEARN:
        raise ImportError("scikit-learn required for training")

    import csv

    # Parse IEDB data
    peptides = []
    ic50_values = []

    with open(data_file, 'r') as f:
        # Try to detect delimiter
        sample = f.read(2048)
        f.seek(0)

        if '\t' in sample:
            delimiter = '\t'
        else:
            delimiter = ','

        reader = csv.DictReader(f, delimiter=delimiter)

        for row in reader:
            # Find peptide column
            peptide = row.get('Description') or row.get('Epitope') or row.get('peptide')
            if not peptide:
                continue

            # Find IC50 column
            ic50_str = row.get('Quantitative measurement') or row.get('IC50') or row.get('ic50_nM')
            if not ic50_str:
                continue

            try:
                ic50 = float(ic50_str)
            except ValueError:
                continue

            # Filter by allele if specified in data
            data_allele = row.get('Allele Name') or row.get('MHC Allele') or row.get('allele')
            if data_allele and allele.replace('HLA-', '') not in data_allele.replace('HLA-', ''):
                continue

            # Validate peptide
            peptide = peptide.upper().strip()
            if len(peptide) < MIN_PEPTIDE_LENGTH or len(peptide) > MAX_PEPTIDE_LENGTH:
                continue
            if not all(aa in AMINO_ACIDS for aa in peptide):
                continue

            peptides.append(peptide)
            ic50_values.append(ic50)

    if len(peptides) < 100:
        raise ValueError(f"Too few valid peptides found: {len(peptides)}")

    # Encode features
    X = encode_peptides(peptides, CANONICAL_LENGTH, use_blosum=True)
    y = np.array([ic50_to_log_ic50(ic50) for ic50 in ic50_values], dtype=np.float32)

    # Train/test split
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_fraction, random_state=RANDOM_STATE
    )

    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Train model
    model = MLPRegressor(
        hidden_layer_sizes=HIDDEN_LAYER_SIZES,
        activation='relu',
        solver='adam',
        alpha=0.001,
        max_iter=MAX_ITER,
        early_stopping=True,
        validation_fraction=0.1,
        random_state=RANDOM_STATE,
        verbose=True
    )

    model.fit(X_train_scaled, y_train)

    # Evaluate
    y_pred = model.predict(X_test_scaled)
    mse = np.mean((y_test - y_pred) ** 2)

    ss_res = np.sum((y_test - y_pred) ** 2)
    ss_tot = np.sum((y_test - np.mean(y_test)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0

    # Generate reference distribution
    n_ref = 5000
    rng = np.random.default_rng(RANDOM_STATE)
    ref_peptides = [generate_random_peptide(9, rng) for _ in range(n_ref)]
    ref_X = scaler.transform(encode_peptides(ref_peptides))
    ref_log_ic50 = model.predict(ref_X)
    ref_ic50s = np.array([log_ic50_to_ic50(x) for x in ref_log_ic50])
    ref_ic50s.sort()

    # Save model
    if output_dir:
        output_path = Path(output_dir)
    else:
        output_path = Path(__file__).parent / 'weights'

    output_path.mkdir(parents=True, exist_ok=True)

    safe_name = allele.replace('*', '_').replace(':', '_')
    model_file = output_path / f'mhc1_{safe_name}.pkl'

    with open(model_file, 'wb') as f:
        pickle.dump({
            'model': model,
            'scaler': scaler,
            'allele': allele,
            'reference_ic50s': ref_ic50s,
        }, f)

    print(f"Model saved to: {model_file}")

    return {
        'mse': float(mse),
        'r2': float(r2),
        'n_train': len(X_train),
        'n_test': len(X_test),
    }


def pretrain_default_alleles(output_dir: Optional[str] = None) -> None:
    """
    Pre-train models for all default alleles with synthetic data.

    This generates the bundled weights that ship with the package.

    Args:
        output_dir: Directory to save model weights.
    """
    if output_dir:
        output_path = Path(output_dir)
    else:
        output_path = Path(__file__).parent / 'weights'

    output_path.mkdir(parents=True, exist_ok=True)

    for allele in ALLELE_ANCHORS.keys():
        print(f"Training model for {allele}...")

        predictor = MHC1Predictor.__new__(MHC1Predictor)
        predictor.allele = allele
        predictor.model = None
        predictor.scaler = None
        predictor._reference_ic50s = None
        predictor._model_dir = output_path

        predictor._train_synthetic_model(n_samples=10000)

        safe_name = allele.replace('*', '_').replace(':', '_')
        model_file = output_path / f'mhc1_{safe_name}.pkl'
        predictor._save_model(model_file)

        print(f"  Saved: {model_file}")

    print(f"\nAll models saved to: {output_path}")


# =============================================================================
# Convenience Functions
# =============================================================================

def predict_binding(
    peptides: Union[str, List[str]],
    allele: str = 'HLA-A*02:01'
) -> Dict[str, Any]:
    """
    Convenience function for quick binding predictions.

    Args:
        peptides: Single peptide or list of peptides.
        allele: HLA allele identifier.

    Returns:
        Dictionary with predictions:
        - 'peptides': List of input peptides
        - 'ic50': IC50 values in nM
        - 'percentile': Percentile ranks
        - 'binding_level': Categorical binding levels
    """
    predictor = MHC1Predictor(allele)

    if isinstance(peptides, str):
        peptides = [peptides]

    ic50 = predictor.predict(peptides)
    percentile = predictor.predict_percentile(peptides)
    levels = predictor.predict_binding_level(peptides)

    return {
        'peptides': peptides,
        'ic50': ic50.tolist(),
        'percentile': percentile.tolist(),
        'binding_level': levels,
        'allele': predictor.allele,
    }


def get_available_alleles() -> List[str]:
    """Return list of supported HLA alleles."""
    return list(ALLELE_ANCHORS.keys())


# =============================================================================
# Command-line Interface
# =============================================================================

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='MHC Class I binding prediction'
    )
    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Predict command
    pred_parser = subparsers.add_parser('predict', help='Predict binding affinity')
    pred_parser.add_argument('peptides', nargs='+', help='Peptide sequences')
    pred_parser.add_argument('--allele', default='HLA-A*02:01', help='HLA allele')

    # Train command
    train_parser = subparsers.add_parser('train', help='Train from IEDB data')
    train_parser.add_argument('data_file', help='Path to IEDB export file')
    train_parser.add_argument('--allele', required=True, help='HLA allele to train')
    train_parser.add_argument('--output', help='Output directory for model')

    # Pretrain command
    pretrain_parser = subparsers.add_parser('pretrain', help='Pre-train default models')
    pretrain_parser.add_argument('--output', help='Output directory for models')

    args = parser.parse_args()

    if args.command == 'predict':
        result = predict_binding(args.peptides, args.allele)
        print(f"Allele: {result['allele']}")
        print(f"{'Peptide':<15} {'IC50 (nM)':<12} {'Percentile':<12} {'Level':<12}")
        print('-' * 51)
        for pep, ic50, pct, level in zip(
            result['peptides'],
            result['ic50'],
            result['percentile'],
            result['binding_level']
        ):
            print(f"{pep:<15} {ic50:<12.1f} {pct:<12.2f} {level:<12}")

    elif args.command == 'train':
        metrics = train_from_iedb(args.data_file, args.allele, args.output)
        print(f"\nTraining results:")
        print(f"  MSE: {metrics['mse']:.4f}")
        print(f"  R2: {metrics['r2']:.4f}")
        print(f"  Training samples: {metrics['n_train']}")
        print(f"  Test samples: {metrics['n_test']}")

    elif args.command == 'pretrain':
        pretrain_default_alleles(args.output)

    else:
        parser.print_help()
