"""
MHC Class II Binding Prediction Model.

This module provides a machine learning-based predictor for MHC Class II
peptide binding. Unlike MHC Class I (which binds fixed 8-11mer peptides),
MHC Class II binds variable length peptides (13-25 amino acids) through
a 9-residue binding core within the peptide.

Key features:
- Handles variable length peptides (13-25 amino acids)
- Identifies optimal 9-mer binding core via sliding window
- Pre-trained models for common HLA-DR alleles
- Returns IC50 predictions and binding core positions

Architecture uses scikit-learn MLPRegressor with features derived from
amino acid properties at anchor positions.
"""
from __future__ import annotations

import base64
import io
import pickle
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler


# =============================================================================
# Constants
# =============================================================================

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INDEX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# Minimum and maximum peptide lengths for MHC-II
MIN_PEPTIDE_LENGTH = 13
MAX_PEPTIDE_LENGTH = 25
CORE_LENGTH = 9


# =============================================================================
# Amino Acid Properties for Feature Encoding
# =============================================================================

# Hydrophobicity (Kyte-Doolittle scale, normalized to [0, 1])
HYDROPHOBICITY = {
    'A': 0.700, 'C': 0.778, 'D': 0.111, 'E': 0.111, 'F': 0.811,
    'G': 0.456, 'H': 0.144, 'I': 1.000, 'K': 0.067, 'L': 0.922,
    'M': 0.711, 'N': 0.111, 'P': 0.322, 'Q': 0.111, 'R': 0.000,
    'S': 0.411, 'T': 0.422, 'V': 0.967, 'W': 0.400, 'Y': 0.356,
}

# Volume (normalized)
VOLUME = {
    'A': 0.167, 'C': 0.287, 'D': 0.313, 'E': 0.380, 'F': 0.593,
    'G': 0.000, 'H': 0.433, 'I': 0.453, 'K': 0.487, 'L': 0.453,
    'M': 0.467, 'N': 0.320, 'P': 0.287, 'Q': 0.387, 'R': 0.547,
    'S': 0.207, 'T': 0.280, 'V': 0.380, 'W': 0.733, 'Y': 0.600,
}

# Polarity
POLARITY = {
    'A': 0.000, 'C': 0.128, 'D': 1.000, 'E': 0.932, 'F': 0.068,
    'G': 0.000, 'H': 0.449, 'I': 0.000, 'K': 0.694, 'L': 0.000,
    'M': 0.068, 'N': 0.593, 'P': 0.221, 'Q': 0.559, 'R': 0.763,
    'S': 0.322, 'T': 0.322, 'V': 0.000, 'W': 0.186, 'Y': 0.373,
}

# Aromaticity (1 for aromatic, 0 otherwise)
AROMATICITY = {
    'A': 0.0, 'C': 0.0, 'D': 0.0, 'E': 0.0, 'F': 1.0,
    'G': 0.0, 'H': 1.0, 'I': 0.0, 'K': 0.0, 'L': 0.0,
    'M': 0.0, 'N': 0.0, 'P': 0.0, 'Q': 0.0, 'R': 0.0,
    'S': 0.0, 'T': 0.0, 'V': 0.0, 'W': 1.0, 'Y': 1.0,
}

# Charge at pH 7 (-1, 0, or 1)
CHARGE = {
    'A': 0.5, 'C': 0.5, 'D': 0.0, 'E': 0.0, 'F': 0.5,
    'G': 0.5, 'H': 0.6, 'I': 0.5, 'K': 1.0, 'L': 0.5,
    'M': 0.5, 'N': 0.5, 'P': 0.5, 'Q': 0.5, 'R': 1.0,
    'S': 0.5, 'T': 0.5, 'V': 0.5, 'W': 0.5, 'Y': 0.5,
}


# =============================================================================
# Allele-Specific Anchor Preferences
# =============================================================================

# Anchor preferences define which amino acids are preferred at key positions
# in the 9-mer binding core. Higher values = stronger preference.

ANCHOR_PREFERENCES = {
    'HLA-DRB1*01:01': {
        # P1: Aromatic/large hydrophobic strongly preferred
        1: {
            'preferred': {'Y': 1.8, 'F': 1.7, 'W': 1.5, 'L': 1.3, 'I': 1.1, 'V': 0.9, 'M': 0.8},
            'tolerated': {'A': 0.3, 'T': 0.2},
            'disfavored': {'D': -1.0, 'E': -1.0, 'K': -0.8, 'R': -0.8, 'P': -1.2, 'G': -0.5},
        },
        # P4: Small residues preferred
        4: {
            'preferred': {'A': 1.0, 'G': 0.9, 'S': 0.8, 'T': 0.7},
            'tolerated': {'V': 0.3, 'N': 0.3, 'Q': 0.2},
            'disfavored': {'W': -0.8, 'F': -0.6, 'Y': -0.5, 'P': -0.7},
        },
        # P6: Small polar preferred
        6: {
            'preferred': {'S': 0.9, 'T': 0.8, 'N': 0.8, 'Q': 0.7, 'A': 0.6},
            'tolerated': {'G': 0.4, 'K': 0.3, 'R': 0.2},
            'disfavored': {'W': -0.6, 'F': -0.5, 'P': -0.8},
        },
        # P9: Hydrophobic preferred
        9: {
            'preferred': {'L': 1.2, 'I': 1.1, 'V': 1.0, 'M': 0.9, 'F': 0.8, 'A': 0.5},
            'tolerated': {'Y': 0.4, 'W': 0.3, 'T': 0.2},
            'disfavored': {'D': -0.8, 'E': -0.7, 'K': -0.5, 'R': -0.5, 'P': -1.0, 'G': -0.3},
        },
    },
    'HLA-DRB1*04:01': {
        # P1: Aromatic strongly preferred
        1: {
            'preferred': {'Y': 1.9, 'F': 1.8, 'W': 1.6, 'L': 1.0, 'I': 0.8},
            'tolerated': {'V': 0.4, 'M': 0.3},
            'disfavored': {'D': -1.2, 'E': -1.1, 'K': -0.9, 'R': -0.8, 'P': -1.3, 'G': -0.6},
        },
        # P4: Polar/small preferred
        4: {
            'preferred': {'T': 1.0, 'S': 0.9, 'Q': 0.8, 'N': 0.7, 'A': 0.6},
            'tolerated': {'G': 0.3, 'M': 0.2},
            'disfavored': {'W': -0.7, 'F': -0.6, 'P': -0.9},
        },
        # P6: Charged/polar preferred
        6: {
            'preferred': {'K': 1.0, 'R': 0.9, 'N': 0.8, 'Q': 0.7, 'H': 0.6},
            'tolerated': {'S': 0.4, 'T': 0.3},
            'disfavored': {'W': -0.5, 'P': -0.9},
        },
        # P9: Basic/polar C-terminus
        9: {
            'preferred': {'R': 1.1, 'K': 1.0, 'Q': 0.8, 'N': 0.7, 'L': 0.5},
            'tolerated': {'H': 0.4, 'Y': 0.3, 'I': 0.3},
            'disfavored': {'D': -0.6, 'E': -0.5, 'P': -1.0},
        },
    },
    'HLA-DRB1*07:01': {
        # P1: Large hydrophobic/aromatic
        1: {
            'preferred': {'L': 1.5, 'I': 1.4, 'V': 1.3, 'F': 1.2, 'Y': 1.1, 'W': 1.0},
            'tolerated': {'M': 0.5, 'A': 0.3},
            'disfavored': {'D': -1.1, 'E': -1.0, 'K': -0.7, 'R': -0.6, 'P': -1.2, 'G': -0.5},
        },
        # P4: Hydrophobic
        4: {
            'preferred': {'M': 1.0, 'L': 0.9, 'I': 0.8, 'V': 0.7, 'F': 0.6},
            'tolerated': {'A': 0.4, 'T': 0.3},
            'disfavored': {'D': -0.6, 'E': -0.5, 'P': -0.8, 'G': -0.4},
        },
        # P6: Small residues
        6: {
            'preferred': {'S': 0.9, 'T': 0.8, 'A': 0.7, 'G': 0.6, 'N': 0.5},
            'tolerated': {'Q': 0.3, 'V': 0.2},
            'disfavored': {'W': -0.6, 'F': -0.5, 'P': -0.8},
        },
        # P9: Aliphatic/small
        9: {
            'preferred': {'L': 1.2, 'I': 1.1, 'V': 1.0, 'A': 0.8, 'M': 0.7},
            'tolerated': {'F': 0.4, 'T': 0.3},
            'disfavored': {'D': -0.7, 'E': -0.6, 'K': -0.4, 'R': -0.4, 'P': -1.0},
        },
    },
}


# =============================================================================
# Feature Encoding
# =============================================================================

def encode_amino_acid(aa: str) -> np.ndarray:
    """
    Encode a single amino acid as a feature vector.

    Features (5 total):
    - Hydrophobicity
    - Volume
    - Polarity
    - Aromaticity
    - Charge

    Args:
        aa: Single letter amino acid code.

    Returns:
        Feature vector of length 5.
    """
    if aa not in AMINO_ACIDS:
        # Return neutral values for unknown amino acids
        return np.array([0.5, 0.5, 0.5, 0.0, 0.5])

    return np.array([
        HYDROPHOBICITY[aa],
        VOLUME[aa],
        POLARITY[aa],
        AROMATICITY[aa],
        CHARGE[aa],
    ])


def encode_core(core: str) -> np.ndarray:
    """
    Encode a 9-mer binding core as a feature vector.

    Includes:
    - Individual position features (9 positions x 5 features = 45)
    - Anchor position one-hot encodings (positions 1, 4, 6, 9)

    Args:
        core: 9-mer peptide sequence.

    Returns:
        Feature vector of length 125 (45 + 20*4 for one-hot anchors).
    """
    if len(core) != CORE_LENGTH:
        raise ValueError(f"Core must be {CORE_LENGTH} residues, got {len(core)}")

    features = []

    # Position-specific features
    for aa in core:
        features.extend(encode_amino_acid(aa))

    # One-hot encoding for anchor positions (1, 4, 6, 9)
    anchor_positions = [0, 3, 5, 8]  # 0-indexed
    for pos in anchor_positions:
        aa = core[pos]
        one_hot = np.zeros(20)
        if aa in AA_TO_INDEX:
            one_hot[AA_TO_INDEX[aa]] = 1.0
        features.extend(one_hot)

    return np.array(features)


def get_anchor_score(core: str, allele: str) -> float:
    """
    Calculate anchor-based score for a binding core.

    Args:
        core: 9-mer binding core sequence.
        allele: HLA allele name.

    Returns:
        Sum of anchor position scores.
    """
    if allele not in ANCHOR_PREFERENCES:
        return 0.0

    prefs = ANCHOR_PREFERENCES[allele]
    score = 0.0

    for pos, pos_prefs in prefs.items():
        aa = core[pos - 1]  # Convert to 0-indexed

        # Check preferred, tolerated, disfavored
        if aa in pos_prefs.get('preferred', {}):
            score += pos_prefs['preferred'][aa]
        elif aa in pos_prefs.get('tolerated', {}):
            score += pos_prefs['tolerated'][aa]
        elif aa in pos_prefs.get('disfavored', {}):
            score += pos_prefs['disfavored'][aa]
        else:
            # Neutral amino acid at this position
            score += 0.0

    return score


# =============================================================================
# Synthetic Training Data Generation
# =============================================================================

def generate_synthetic_binder(allele: str, rng: np.random.Generator) -> Tuple[str, float]:
    """
    Generate a synthetic binding peptide core for an allele.

    Creates sequences that match anchor preferences with added noise.

    Args:
        allele: HLA allele name.
        rng: NumPy random generator.

    Returns:
        Tuple of (9-mer core sequence, binding score).
    """
    if allele not in ANCHOR_PREFERENCES:
        raise ValueError(f"Unknown allele: {allele}")

    prefs = ANCHOR_PREFERENCES[allele]
    core = list('AAAAAAAAA')  # Start with alanines
    base_score = 0.0

    # Fill anchor positions based on preferences
    for pos, pos_prefs in prefs.items():
        preferred = pos_prefs.get('preferred', {})
        if preferred:
            # Weight selection by preference scores
            aas = list(preferred.keys())
            weights = np.array([preferred[aa] for aa in aas])
            weights = np.exp(weights) / np.sum(np.exp(weights))  # Softmax

            selected_aa = rng.choice(aas, p=weights)
            core[pos - 1] = selected_aa
            base_score += preferred[selected_aa]

    # Fill non-anchor positions with random amino acids
    non_anchor_positions = [i for i in range(9) if (i + 1) not in prefs]
    for pos in non_anchor_positions:
        core[pos] = rng.choice(list(AMINO_ACIDS))

    # Add noise to score
    noise = rng.normal(0, 0.3)
    final_score = base_score + noise

    return ''.join(core), final_score


def generate_synthetic_nonbinder(allele: str, rng: np.random.Generator) -> Tuple[str, float]:
    """
    Generate a synthetic non-binding peptide core.

    Creates sequences with disfavored residues at anchor positions.

    Args:
        allele: HLA allele name.
        rng: NumPy random generator.

    Returns:
        Tuple of (9-mer core sequence, binding score).
    """
    if allele not in ANCHOR_PREFERENCES:
        raise ValueError(f"Unknown allele: {allele}")

    prefs = ANCHOR_PREFERENCES[allele]
    core = list('AAAAAAAAA')
    base_score = 0.0

    # Fill anchor positions with disfavored residues
    for pos, pos_prefs in prefs.items():
        disfavored = pos_prefs.get('disfavored', {})
        if disfavored and rng.random() > 0.3:  # 70% chance of disfavored
            aas = list(disfavored.keys())
            selected_aa = rng.choice(aas)
            core[pos - 1] = selected_aa
            base_score += disfavored[selected_aa]
        else:
            # Random amino acid
            core[pos - 1] = rng.choice(list(AMINO_ACIDS))

    # Fill non-anchor positions
    non_anchor_positions = [i for i in range(9) if (i + 1) not in prefs]
    for pos in non_anchor_positions:
        core[pos] = rng.choice(list(AMINO_ACIDS))

    noise = rng.normal(0, 0.3)
    final_score = base_score + noise

    return ''.join(core), final_score


def generate_training_data(
    allele: str,
    n_samples: int = 2000,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate synthetic training data for an allele.

    Creates balanced dataset of binders and non-binders based on
    anchor position preferences.

    Args:
        allele: HLA allele name.
        n_samples: Total number of samples to generate.
        seed: Random seed for reproducibility.

    Returns:
        Tuple of (feature matrix, score array).
    """
    rng = np.random.default_rng(seed)

    cores = []
    scores = []

    n_binders = n_samples // 2
    n_nonbinders = n_samples - n_binders

    # Generate binders
    for _ in range(n_binders):
        core, score = generate_synthetic_binder(allele, rng)
        cores.append(core)
        scores.append(score)

    # Generate non-binders
    for _ in range(n_nonbinders):
        core, score = generate_synthetic_nonbinder(allele, rng)
        cores.append(core)
        scores.append(score)

    # Encode features
    X = np.array([encode_core(core) for core in cores])
    y = np.array(scores)

    return X, y


# =============================================================================
# Score to IC50 Conversion
# =============================================================================

def score_to_ic50(score: float, scale: float = 500.0) -> float:
    """
    Convert binding score to estimated IC50 (nM).

    Higher scores indicate better binding (lower IC50).

    Args:
        score: Model binding score.
        scale: Scaling factor (default 500 nM at score=0).

    Returns:
        Estimated IC50 in nanomolar.
    """
    ic50 = scale * np.exp(-score * 0.6)
    return float(max(1.0, min(50000.0, ic50)))


def ic50_to_binding_level(ic50: float) -> str:
    """
    Classify binding strength based on IC50.

    Args:
        ic50: IC50 value in nanomolar.

    Returns:
        Binding level classification.
    """
    if ic50 < 50:
        return "Strong"
    elif ic50 < 500:
        return "Moderate"
    elif ic50 < 5000:
        return "Weak"
    else:
        return "Non-binder"


# =============================================================================
# Prediction Result Dataclass
# =============================================================================

@dataclass
class MHC2PredictionResult:
    """Result from MHC Class II binding prediction."""
    peptide: str
    core: str
    core_start: int  # 0-indexed position of core within peptide
    score: float
    ic50: float
    binding_level: str
    allele: str
    flanking_n: str = ""  # N-terminal flanking region
    flanking_c: str = ""  # C-terminal flanking region

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'peptide': self.peptide,
            'core': self.core,
            'core_start': self.core_start,
            'score': self.score,
            'ic50': self.ic50,
            'binding_level': self.binding_level,
            'allele': self.allele,
            'flanking_n': self.flanking_n,
            'flanking_c': self.flanking_c,
        }


# =============================================================================
# Pre-trained Model Weights
# =============================================================================

def _create_pretrained_model(allele: str) -> Tuple[MLPRegressor, StandardScaler]:
    """
    Create and train a model for the specified allele.

    Uses synthetic training data based on anchor preferences.

    Args:
        allele: HLA allele name.

    Returns:
        Tuple of (trained MLPRegressor, fitted StandardScaler).
    """
    # Generate training data
    X, y = generate_training_data(allele, n_samples=2000, seed=42)

    # Scale features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Train MLP
    model = MLPRegressor(
        hidden_layer_sizes=(64, 32),
        activation='relu',
        solver='adam',
        max_iter=500,
        random_state=42,
        early_stopping=True,
        validation_fraction=0.1,
        n_iter_no_change=20,
    )
    model.fit(X_scaled, y)

    return model, scaler


# Cache for pre-trained models
_MODEL_CACHE: Dict[str, Tuple[MLPRegressor, StandardScaler]] = {}


def _get_pretrained_model(allele: str) -> Tuple[MLPRegressor, StandardScaler]:
    """
    Get or create a pre-trained model for the allele.

    Models are cached after first creation.

    Args:
        allele: HLA allele name.

    Returns:
        Tuple of (MLPRegressor, StandardScaler).
    """
    if allele not in _MODEL_CACHE:
        _MODEL_CACHE[allele] = _create_pretrained_model(allele)
    return _MODEL_CACHE[allele]


# =============================================================================
# Main Predictor Class
# =============================================================================

class MHC2Predictor:
    """
    MHC Class II binding prediction using neural network model.

    Predicts peptide binding to MHC Class II molecules (HLA-DR alleles).
    Handles variable length peptides (13-25 amino acids) by finding the
    optimal 9-mer binding core through a sliding window approach.

    Example:
        >>> predictor = MHC2Predictor(allele='HLA-DRB1*01:01')
        >>> results = predictor.predict(['AGFKGEQGPKGEPGAA', 'YSFPLGAKVILDMFRG'])
        >>> for r in results:
        ...     print(f"{r['peptide']}: IC50={r['ic50']:.1f} nM, core={r['core']}")

    Attributes:
        allele: HLA-DR allele for predictions.
        model: Trained MLPRegressor model.
        scaler: StandardScaler for feature normalization.
    """

    # Supported alleles with pre-trained models
    SUPPORTED_ALLELES = [
        'HLA-DRB1*01:01',
        'HLA-DRB1*04:01',
        'HLA-DRB1*07:01',
    ]

    def __init__(self, allele: str = 'HLA-DRB1*01:01'):
        """
        Initialize the MHC-II predictor.

        Args:
            allele: HLA-DR allele name. Must be one of SUPPORTED_ALLELES.

        Raises:
            ValueError: If allele is not supported.
        """
        if allele not in self.SUPPORTED_ALLELES:
            raise ValueError(
                f"Unsupported allele: {allele}. "
                f"Supported alleles: {self.SUPPORTED_ALLELES}"
            )

        self.allele = allele
        self.model, self.scaler = _get_pretrained_model(allele)

    def find_binding_core(self, peptide: str) -> Tuple[int, str, float]:
        """
        Find the 9-mer binding core with the best binding score.

        Slides a 9-residue window across the peptide and scores each
        potential binding core. Returns the position and sequence of
        the core with the highest predicted binding affinity.

        Args:
            peptide: Amino acid sequence (13-25 residues).

        Returns:
            Tuple of (core_start_position, core_sequence, core_score).
            Position is 0-indexed within the peptide.

        Raises:
            ValueError: If peptide length is outside valid range.
        """
        peptide = peptide.upper()

        if len(peptide) < MIN_PEPTIDE_LENGTH:
            raise ValueError(
                f"Peptide too short: {len(peptide)} < {MIN_PEPTIDE_LENGTH}"
            )
        if len(peptide) > MAX_PEPTIDE_LENGTH:
            raise ValueError(
                f"Peptide too long: {len(peptide)} > {MAX_PEPTIDE_LENGTH}"
            )

        best_pos = 0
        best_core = peptide[:CORE_LENGTH]
        best_score = float('-inf')

        # Slide window across peptide
        n_positions = len(peptide) - CORE_LENGTH + 1
        for i in range(n_positions):
            core = peptide[i:i + CORE_LENGTH]

            # Skip cores with non-standard amino acids
            if not all(aa in AMINO_ACIDS for aa in core):
                continue

            # Encode and score core
            features = encode_core(core).reshape(1, -1)
            features_scaled = self.scaler.transform(features)
            score = float(self.model.predict(features_scaled)[0])

            # Add anchor preference bonus
            anchor_bonus = get_anchor_score(core, self.allele) * 0.2
            score += anchor_bonus

            if score > best_score:
                best_score = score
                best_pos = i
                best_core = core

        return best_pos, best_core, best_score

    def predict_single(self, peptide: str) -> MHC2PredictionResult:
        """
        Predict MHC-II binding for a single peptide.

        Args:
            peptide: Amino acid sequence (13-25 residues).

        Returns:
            MHC2PredictionResult with binding prediction.
        """
        peptide = peptide.upper()

        # Validate peptide
        if len(peptide) < MIN_PEPTIDE_LENGTH:
            raise ValueError(
                f"Peptide '{peptide}' too short: {len(peptide)} < {MIN_PEPTIDE_LENGTH}"
            )
        if len(peptide) > MAX_PEPTIDE_LENGTH:
            raise ValueError(
                f"Peptide '{peptide}' too long: {len(peptide)} > {MAX_PEPTIDE_LENGTH}"
            )

        # Find best binding core
        core_start, core, score = self.find_binding_core(peptide)

        # Convert score to IC50
        ic50 = score_to_ic50(score)
        binding_level = ic50_to_binding_level(ic50)

        # Extract flanking regions
        flanking_n = peptide[:core_start]
        flanking_c = peptide[core_start + CORE_LENGTH:]

        return MHC2PredictionResult(
            peptide=peptide,
            core=core,
            core_start=core_start,
            score=score,
            ic50=ic50,
            binding_level=binding_level,
            allele=self.allele,
            flanking_n=flanking_n,
            flanking_c=flanking_c,
        )

    def predict(self, peptides: List[str]) -> List[Dict[str, Any]]:
        """
        Predict MHC-II binding for multiple peptides.

        Args:
            peptides: List of amino acid sequences (13-25 residues each).

        Returns:
            List of prediction dictionaries with keys:
                - peptide: Input peptide sequence
                - core: 9-mer binding core
                - core_start: Position of core in peptide (0-indexed)
                - score: Raw binding score
                - ic50: Estimated IC50 in nanomolar
                - binding_level: Classification (Strong/Moderate/Weak/Non-binder)
                - allele: HLA allele used
                - flanking_n: N-terminal flanking region
                - flanking_c: C-terminal flanking region
        """
        results = []

        for peptide in peptides:
            try:
                result = self.predict_single(peptide)
                results.append(result.to_dict())
            except ValueError as e:
                # Include error information in result
                results.append({
                    'peptide': peptide,
                    'core': '',
                    'core_start': -1,
                    'score': float('nan'),
                    'ic50': float('nan'),
                    'binding_level': 'Error',
                    'allele': self.allele,
                    'error': str(e),
                })

        return results

    def predict_from_protein(
        self,
        sequence: str,
        peptide_length: int = 15,
        overlap: int = 1,
    ) -> List[Dict[str, Any]]:
        """
        Generate overlapping peptides from a protein and predict binding.

        Slides a window across the protein sequence to generate peptides,
        then predicts MHC-II binding for each.

        Args:
            sequence: Protein amino acid sequence.
            peptide_length: Length of peptides to generate (13-25).
            overlap: Step size for sliding window (default 1).

        Returns:
            List of prediction dictionaries with additional 'position' key
            indicating the start position in the protein (1-indexed).
        """
        sequence = sequence.upper()

        if peptide_length < MIN_PEPTIDE_LENGTH or peptide_length > MAX_PEPTIDE_LENGTH:
            raise ValueError(
                f"Peptide length must be {MIN_PEPTIDE_LENGTH}-{MAX_PEPTIDE_LENGTH}"
            )

        if len(sequence) < peptide_length:
            raise ValueError(
                f"Protein sequence ({len(sequence)} aa) shorter than "
                f"peptide length ({peptide_length})"
            )

        # Generate peptides
        peptides = []
        positions = []
        step = peptide_length - overlap + 1 if overlap < peptide_length else 1

        for i in range(0, len(sequence) - peptide_length + 1, step):
            peptide = sequence[i:i + peptide_length]
            # Skip peptides with non-standard amino acids
            if all(aa in AMINO_ACIDS for aa in peptide):
                peptides.append(peptide)
                positions.append(i + 1)  # 1-indexed

        # Predict binding
        results = self.predict(peptides)

        # Add position information
        for result, pos in zip(results, positions):
            result['position'] = pos

        return results

    def get_binders(
        self,
        results: List[Dict[str, Any]],
        ic50_threshold: float = 500.0,
    ) -> List[Dict[str, Any]]:
        """
        Filter predictions to return only binders below IC50 threshold.

        Args:
            results: List of prediction dictionaries from predict().
            ic50_threshold: Maximum IC50 in nM for binders (default 500).

        Returns:
            Filtered list of predictions sorted by IC50 (best first).
        """
        binders = [
            r for r in results
            if not np.isnan(r.get('ic50', float('nan'))) and r['ic50'] <= ic50_threshold
        ]
        return sorted(binders, key=lambda x: x['ic50'])

    @classmethod
    def list_alleles(cls) -> List[str]:
        """
        List all supported HLA-DR alleles.

        Returns:
            List of allele names.
        """
        return cls.SUPPORTED_ALLELES.copy()

    def __repr__(self) -> str:
        return f"MHC2Predictor(allele='{self.allele}')"


# =============================================================================
# Convenience Functions
# =============================================================================

def predict_mhc2_binding(
    peptides: List[str],
    allele: str = 'HLA-DRB1*01:01',
) -> List[Dict[str, Any]]:
    """
    Convenience function for MHC-II binding prediction.

    Args:
        peptides: List of peptide sequences (13-25 amino acids each).
        allele: HLA-DR allele name.

    Returns:
        List of prediction dictionaries.

    Example:
        >>> results = predict_mhc2_binding(['AGFKGEQGPKGEPGAA'])
        >>> print(f"IC50: {results[0]['ic50']:.1f} nM")
    """
    predictor = MHC2Predictor(allele=allele)
    return predictor.predict(peptides)


def predict_protein_mhc2(
    sequence: str,
    allele: str = 'HLA-DRB1*01:01',
    peptide_length: int = 15,
    ic50_threshold: float = 500.0,
) -> List[Dict[str, Any]]:
    """
    Predict MHC-II epitopes from a protein sequence.

    Args:
        sequence: Protein amino acid sequence.
        allele: HLA-DR allele name.
        peptide_length: Length of peptides to generate.
        ic50_threshold: IC50 threshold for reporting binders.

    Returns:
        List of predicted binders sorted by IC50.
    """
    predictor = MHC2Predictor(allele=allele)
    results = predictor.predict_from_protein(sequence, peptide_length=peptide_length)
    return predictor.get_binders(results, ic50_threshold=ic50_threshold)
