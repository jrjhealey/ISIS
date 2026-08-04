"""
T-cell epitope prediction methods.

This module provides prediction methods for T-cell epitopes including:
- MHC Class I binding prediction (8-11mer peptides)
- MHC Class II binding prediction (15mer peptides with 9mer core)
- Proteasomal cleavage prediction
- TAP transport prediction
- Combined pipeline scoring

Methods use Position-Specific Scoring Matrices (PSSMs) derived from
published IEDB/SMM data and known anchor residue preferences.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .base import EpitopeMethod, MethodCategory, MethodResult


# Standard amino acids
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Default score for unknown amino acids
UNKNOWN_AA_SCORE = -1.0


# =============================================================================
# MHC Class I PSSMs (9-mer)
# Position-specific scoring matrices based on published anchor preferences
# =============================================================================

def _build_pssm(anchors: Dict[int, Dict[str, float]], length: int = 9) -> Dict[int, Dict[str, float]]:
    """
    Build a complete PSSM from anchor specifications.

    Non-anchor positions get neutral scores based on amino acid frequencies.
    Anchor positions use specified preferences.

    Args:
        anchors: Dict mapping positions to amino acid score dicts.
        length: Peptide length.

    Returns:
        Complete PSSM with scores for all positions and amino acids.
    """
    # Base scores for non-anchor positions (slight preference for common residues)
    base_scores = {
        'A': 0.05, 'C': -0.1, 'D': 0.0, 'E': 0.0, 'F': -0.05,
        'G': 0.05, 'H': 0.0, 'I': -0.05, 'K': 0.1, 'L': -0.05,
        'M': 0.0, 'N': 0.05, 'P': 0.0, 'Q': 0.05, 'R': 0.1,
        'S': 0.05, 'T': 0.05, 'V': -0.05, 'W': -0.1, 'Y': 0.0
    }

    pssm = {}
    for pos in range(1, length + 1):
        if pos in anchors:
            # Use anchor-specific scores
            pssm[pos] = anchors[pos].copy()
            # Fill in missing amino acids with low scores
            for aa in AMINO_ACIDS:
                if aa not in pssm[pos]:
                    pssm[pos][aa] = -0.5
        else:
            # Use base scores for non-anchor positions
            pssm[pos] = base_scores.copy()

    return pssm


# HLA-A*02:01: Anchor positions 2 (L/M) and 9 (V/L)
HLA_A0201_ANCHORS = {
    2: {'L': 1.5, 'M': 1.2, 'I': 0.8, 'V': 0.5, 'A': 0.2, 'T': -0.2, 'Q': -0.3},
    9: {'V': 1.8, 'L': 1.5, 'I': 1.2, 'A': 0.3, 'T': 0.1, 'M': 0.5}
}
HLA_A0201_PSSM_9 = _build_pssm(HLA_A0201_ANCHORS, 9)

# HLA-A*01:01: Anchor positions 2 (T/S) and 9 (Y)
HLA_A0101_ANCHORS = {
    2: {'T': 1.4, 'S': 1.2, 'M': 0.3, 'L': -0.2, 'I': -0.3, 'V': -0.2},
    3: {'D': 0.8, 'E': 0.6, 'A': 0.2},  # Secondary preference
    9: {'Y': 1.8, 'F': 0.8, 'W': 0.5, 'L': -0.3, 'V': -0.4}
}
HLA_A0101_PSSM_9 = _build_pssm(HLA_A0101_ANCHORS, 9)

# HLA-B*07:02: Anchor positions 2 (P) and 9 (L)
HLA_B0702_ANCHORS = {
    2: {'P': 2.0, 'A': 0.3, 'S': 0.1, 'L': -0.3, 'V': -0.4},
    9: {'L': 1.6, 'I': 1.0, 'V': 0.8, 'M': 0.6, 'F': 0.4, 'A': 0.2}
}
HLA_B0702_PSSM_9 = _build_pssm(HLA_B0702_ANCHORS, 9)

# Collected MHC-I PSSMs by allele
MHC1_PSSM = {
    'HLA-A*02:01': {9: HLA_A0201_PSSM_9},
    'HLA-A*01:01': {9: HLA_A0101_PSSM_9},
    'HLA-B*07:02': {9: HLA_B0702_PSSM_9},
}


# =============================================================================
# MHC Class II PSSMs (15-mer with 9-mer core)
# =============================================================================

# HLA-DRB1*01:01: P1 (aromatic/large hydrophobic), P4, P6, P9 anchors
HLA_DRB1_0101_CORE_ANCHORS = {
    1: {'F': 1.8, 'Y': 1.6, 'W': 1.4, 'L': 1.2, 'I': 0.8, 'V': 0.5, 'M': 0.4},
    4: {'M': 0.8, 'L': 0.6, 'I': 0.4, 'V': 0.3, 'T': 0.2},
    6: {'N': 0.6, 'Q': 0.5, 'S': 0.4, 'T': 0.4, 'K': 0.3},
    9: {'L': 1.2, 'I': 0.9, 'V': 0.7, 'M': 0.5, 'A': 0.2}
}
HLA_DRB1_0101_PSSM = _build_pssm(HLA_DRB1_0101_CORE_ANCHORS, 9)

# HLA-DRB1*04:01: P1 (aromatic), P4 (hydrophobic), P6, P9 anchors
HLA_DRB1_0401_CORE_ANCHORS = {
    1: {'Y': 1.8, 'F': 1.6, 'W': 1.4, 'L': 0.8, 'I': 0.5, 'V': 0.3},
    4: {'T': 0.9, 'S': 0.7, 'M': 0.5, 'L': 0.3, 'V': 0.2},
    6: {'K': 0.7, 'R': 0.6, 'N': 0.5, 'Q': 0.4, 'H': 0.3},
    9: {'R': 1.0, 'K': 0.8, 'Q': 0.6, 'N': 0.5, 'L': 0.3}
}
HLA_DRB1_0401_PSSM = _build_pssm(HLA_DRB1_0401_CORE_ANCHORS, 9)

# Collected MHC-II core PSSMs by allele
MHC2_CORE_PSSM = {
    'HLA-DRB1*01:01': HLA_DRB1_0101_PSSM,
    'HLA-DRB1*04:01': HLA_DRB1_0401_PSSM,
}


# =============================================================================
# Proteasomal Cleavage PSSM
# Based on Netchop approach - scores for positions P6 to P6'
# =============================================================================

# Cleavage occurs between P1 and P1'
# Positive values favor cleavage, negative values disfavor
PROTEASOME_CLEAVAGE_PSSM = {
    # P6 to P1 (before cleavage site)
    -6: {'L': 0.2, 'V': 0.1, 'I': 0.1, 'F': 0.1, 'Y': 0.05, 'W': 0.05},
    -5: {'E': 0.1, 'D': 0.1, 'K': 0.1, 'R': 0.1},
    -4: {'A': 0.1, 'G': 0.05, 'S': 0.05, 'T': 0.05},
    -3: {'L': 0.15, 'V': 0.1, 'I': 0.1, 'F': 0.1, 'M': 0.1},
    -2: {'E': 0.2, 'D': 0.15, 'Q': 0.1, 'N': 0.1, 'K': 0.1},
    -1: {  # P1 - C-terminal cleavage preference (strong)
        'L': 1.2, 'Y': 1.0, 'F': 0.9, 'M': 0.8, 'W': 0.7,
        'K': 0.6, 'R': 0.5, 'A': 0.4, 'V': 0.3, 'I': 0.3,
        'H': 0.2, 'Q': 0.1, 'N': 0.0, 'T': 0.0, 'S': -0.1,
        'C': -0.2, 'G': -0.3, 'P': -0.8, 'E': -0.4, 'D': -0.5
    },
    # P1' to P6' (after cleavage site)
    1: {  # P1' - N-terminal of cleaved peptide
        'A': 0.3, 'S': 0.2, 'G': 0.2, 'T': 0.1, 'V': 0.1,
        'P': -0.5, 'E': -0.2, 'D': -0.2, 'K': 0.1, 'R': 0.1
    },
    2: {'L': 0.1, 'V': 0.1, 'I': 0.1, 'A': 0.05},
    3: {'E': 0.1, 'D': 0.1, 'K': 0.05, 'R': 0.05},
    4: {'A': 0.05, 'G': 0.05, 'S': 0.05},
    5: {'L': 0.05, 'V': 0.05, 'I': 0.05},
    6: {'K': 0.05, 'R': 0.05, 'E': 0.05},
}


# =============================================================================
# TAP Transport PSSM
# N-terminal 9 residues contribute most to TAP binding
# =============================================================================

TAP_TRANSPORT_PSSM = {
    1: {  # Position 1 (N-terminus) - important
        'R': 1.0, 'K': 0.8, 'H': 0.5, 'Y': 0.4, 'F': 0.3,
        'W': 0.2, 'L': 0.1, 'I': 0.1, 'V': 0.0, 'M': 0.0,
        'A': -0.1, 'G': -0.2, 'S': -0.1, 'T': -0.1, 'C': -0.2,
        'P': -0.4, 'N': -0.1, 'Q': -0.1, 'E': -0.3, 'D': -0.4
    },
    2: {
        'L': 0.4, 'I': 0.3, 'V': 0.2, 'M': 0.2, 'F': 0.3,
        'Y': 0.2, 'W': 0.2, 'A': 0.1, 'K': 0.1, 'R': 0.1
    },
    3: {
        'Y': 0.3, 'F': 0.2, 'W': 0.2, 'L': 0.1, 'I': 0.1,
        'V': 0.1, 'K': 0.1, 'R': 0.1
    },
    4: {'L': 0.1, 'I': 0.1, 'V': 0.1, 'F': 0.1, 'Y': 0.1},
    5: {'K': 0.1, 'R': 0.1, 'E': 0.05, 'D': 0.05},
    6: {'L': 0.1, 'V': 0.1, 'I': 0.1},
    7: {'A': 0.05, 'G': 0.05, 'S': 0.05},
    8: {'L': 0.1, 'V': 0.1, 'I': 0.1},
    9: {  # C-terminus - very important for TAP
        'L': 1.2, 'F': 1.0, 'Y': 0.9, 'M': 0.8, 'I': 0.7,
        'V': 0.5, 'A': 0.2, 'R': 0.6, 'K': 0.5, 'W': 0.4,
        'H': 0.2, 'T': 0.0, 'S': -0.1, 'C': -0.2, 'N': -0.1,
        'Q': -0.1, 'G': -0.3, 'P': -0.6, 'E': -0.4, 'D': -0.5
    },
}


# =============================================================================
# Available Alleles Registry
# =============================================================================

AVAILABLE_ALLELES = {
    'MHC-I': {
        'HLA-A*02:01': {
            'description': 'Most common Caucasian allele, anchor P2=L/M, P9=V/L',
            'peptide_lengths': [8, 9, 10, 11],
            'population_frequency': 0.28,
        },
        'HLA-A*01:01': {
            'description': 'Common allele, anchor P2=T/S, P9=Y',
            'peptide_lengths': [8, 9, 10, 11],
            'population_frequency': 0.16,
        },
        'HLA-B*07:02': {
            'description': 'Common allele, anchor P2=P, P9=L',
            'peptide_lengths': [8, 9, 10, 11],
            'population_frequency': 0.12,
        },
    },
    'MHC-II': {
        'HLA-DRB1*01:01': {
            'description': 'Common DR1 allele, P1=aromatic/hydrophobic',
            'peptide_length': 15,
            'core_length': 9,
            'population_frequency': 0.10,
        },
        'HLA-DRB1*04:01': {
            'description': 'Common DR4 allele, associated with RA',
            'peptide_length': 15,
            'core_length': 9,
            'population_frequency': 0.08,
        },
    },
}


# =============================================================================
# Scoring Functions
# =============================================================================

def score_peptide_pssm(peptide: str, pssm: Dict[int, Dict[str, float]]) -> float:
    """
    Score a peptide using a PSSM.

    Args:
        peptide: Amino acid sequence.
        pssm: Position-specific scoring matrix (1-indexed positions).

    Returns:
        Sum of position-specific scores.
    """
    score = 0.0
    for i, aa in enumerate(peptide):
        pos = i + 1
        if pos in pssm:
            score += pssm[pos].get(aa, UNKNOWN_AA_SCORE)
    return score


def score_to_ic50(score: float, scale: float = 500.0) -> float:
    """
    Convert PSSM score to approximate IC50 (nM).

    Higher scores mean better binding (lower IC50).
    Uses exponential transformation.

    Args:
        score: PSSM score.
        scale: Scaling factor for IC50 range.

    Returns:
        Estimated IC50 in nanomolar.
    """
    # Transform score to IC50 range (50-50000 nM typical)
    ic50 = scale * np.exp(-score * 0.5)
    return max(1.0, min(50000.0, ic50))


def ic50_to_percentile_rank(ic50: float) -> float:
    """
    Convert IC50 to percentile rank (lower is better binding).

    Args:
        ic50: IC50 value in nanomolar.

    Returns:
        Approximate percentile rank (0-100).
    """
    # Approximate percentile based on IC50 thresholds
    if ic50 < 50:
        return 0.1
    elif ic50 < 500:
        return 1.0 + (ic50 - 50) / 450 * 4.0
    elif ic50 < 5000:
        return 5.0 + (ic50 - 500) / 4500 * 45.0
    else:
        return 50.0 + min(50.0, (ic50 - 5000) / 45000 * 50.0)


# =============================================================================
# T-cell Predictor Class
# =============================================================================

@dataclass
class TcellPredictionResult:
    """Result from T-cell epitope prediction."""
    scores: List[float]
    peptides: List[str]
    threshold: float
    positions: Optional[List[int]] = None
    ic50_values: Optional[List[float]] = None
    percentile_ranks: Optional[List[float]] = None
    metadata: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        result = {
            'scores': self.scores,
            'peptides': self.peptides,
            'threshold': self.threshold,
        }
        if self.positions is not None:
            result['positions'] = self.positions
        if self.ic50_values is not None:
            result['ic50_values'] = self.ic50_values
        if self.percentile_ranks is not None:
            result['percentile_ranks'] = self.percentile_ranks
        if self.metadata is not None:
            result['metadata'] = self.metadata
        return result

    def get_binders(self, use_threshold: Optional[float] = None) -> List[Dict[str, Any]]:
        """
        Get peptides predicted to bind above threshold.

        Args:
            use_threshold: Override the default threshold.

        Returns:
            List of dicts with peptide info for predicted binders.
        """
        thresh = use_threshold if use_threshold is not None else self.threshold
        binders = []
        for i, (score, peptide) in enumerate(zip(self.scores, self.peptides)):
            if score >= thresh:
                binder = {
                    'peptide': peptide,
                    'score': score,
                    'position': self.positions[i] if self.positions else i + 1,
                }
                if self.ic50_values:
                    binder['ic50'] = self.ic50_values[i]
                if self.percentile_ranks:
                    binder['percentile_rank'] = self.percentile_ranks[i]
                binders.append(binder)
        return sorted(binders, key=lambda x: x['score'], reverse=True)


class TcellPredictor(EpitopeMethod):
    """
    T-cell epitope predictor using PSSM-based methods.

    Provides prediction methods for:
    - MHC Class I binding (CD8+ T-cell epitopes)
    - MHC Class II binding (CD4+ T-cell epitopes)
    - Proteasomal cleavage
    - TAP transport
    - Combined pipeline scoring

    Example:
        >>> predictor = TcellPredictor()
        >>> result = predictor.predict_mhc1("MKTAYIAKQRQISFVK", allele='HLA-A*02:01')
        >>> print(result.get_binders())
    """

    name = "T-cell Epitope Predictor"
    key = "tcell"
    category = MethodCategory.TCELL
    description = "PSSM-based T-cell epitope prediction for MHC binding"
    citation = "Based on IEDB/SMM methodologies"
    default_window = 9
    default_threshold = 0.5

    # Thresholds for predictions
    MHC1_THRESHOLD = 0.5       # PSSM score
    MHC1_IC50_THRESHOLD = 500  # nM - strong binder
    MHC2_THRESHOLD = 0.5
    CLEAVAGE_THRESHOLD = 0.8
    TAP_THRESHOLD = 0.5

    def __init__(self):
        """Initialize the T-cell predictor."""
        self.mhc1_pssm = MHC1_PSSM
        self.mhc2_pssm = MHC2_CORE_PSSM
        self.cleavage_pssm = PROTEASOME_CLEAVAGE_PSSM
        self.tap_pssm = TAP_TRANSPORT_PSSM

    def predict(self, sequence: str, **kwargs) -> MethodResult:
        """
        General prediction interface (required by base class).

        Uses the combined pipeline by default.
        """
        allele = kwargs.get('allele', 'HLA-A*02:01')
        result = self.predict_pipeline(sequence, allele=allele)

        # Convert to MethodResult format
        scores_array = np.array(result.scores)
        positions_array = np.array(result.positions) if result.positions else np.arange(1, len(result.scores) + 1)

        epitopes = []
        for binder in result.get_binders():
            epitopes.append({
                'start': binder['position'],
                'end': binder['position'] + len(binder['peptide']) - 1,
                'sequence': binder['peptide'],
                'score': binder['score'],
            })

        return MethodResult(
            scores=scores_array,
            positions=positions_array,
            epitopes=epitopes,
            metadata=result.metadata or {},
        )

    def predict_mhc1(
        self,
        sequence: str,
        allele: str = 'HLA-A*02:01',
        peptide_length: int = 9,
    ) -> TcellPredictionResult:
        """
        Predict MHC Class I binding peptides.

        Generates all possible peptides of the specified length and scores
        them using allele-specific PSSMs.

        Args:
            sequence: Protein sequence (amino acid one-letter codes).
            allele: HLA allele name (e.g., 'HLA-A*02:01').
            peptide_length: Length of peptides to generate (8-11, default 9).

        Returns:
            TcellPredictionResult with scores, peptides, and binding predictions.

        Raises:
            ValueError: If allele is not supported or sequence is too short.
        """
        sequence = self.validate_sequence(sequence)

        if allele not in self.mhc1_pssm:
            available = list(self.mhc1_pssm.keys())
            raise ValueError(f"Unsupported MHC-I allele: {allele}. Available: {available}")

        if peptide_length < 8 or peptide_length > 11:
            raise ValueError(f"Peptide length must be 8-11, got {peptide_length}")

        if len(sequence) < peptide_length:
            raise ValueError(f"Sequence length ({len(sequence)}) must be >= peptide length ({peptide_length})")

        # Get PSSM - use 9-mer PSSM for all lengths with position adjustment
        allele_pssm = self.mhc1_pssm[allele]

        # Use native length PSSM if available, otherwise adapt 9-mer
        if peptide_length in allele_pssm:
            pssm = allele_pssm[peptide_length]
        else:
            # Adapt 9-mer PSSM: keep anchor positions, adjust middle
            pssm = self._adapt_pssm_length(allele_pssm[9], 9, peptide_length)

        # Generate and score all peptides
        peptides = []
        scores = []
        positions = []
        ic50_values = []
        percentile_ranks = []

        n_peptides = len(sequence) - peptide_length + 1
        for i in range(n_peptides):
            peptide = sequence[i:i + peptide_length]
            peptides.append(peptide)
            positions.append(i + 1)

            score = score_peptide_pssm(peptide, pssm)
            scores.append(score)

            ic50 = score_to_ic50(score)
            ic50_values.append(ic50)
            percentile_ranks.append(ic50_to_percentile_rank(ic50))

        return TcellPredictionResult(
            scores=scores,
            peptides=peptides,
            threshold=self.MHC1_THRESHOLD,
            positions=positions,
            ic50_values=ic50_values,
            percentile_ranks=percentile_ranks,
            metadata={
                'allele': allele,
                'peptide_length': peptide_length,
                'method': 'PSSM',
                'ic50_threshold': self.MHC1_IC50_THRESHOLD,
            }
        )

    def _adapt_pssm_length(
        self,
        pssm_9: Dict[int, Dict[str, float]],
        original_length: int,
        target_length: int
    ) -> Dict[int, Dict[str, float]]:
        """
        Adapt a 9-mer PSSM to a different peptide length.

        Keeps anchor positions (2 and C-terminus) and adjusts middle positions.
        """
        new_pssm = {}

        if target_length == 8:
            # Remove middle position
            new_pssm[1] = pssm_9[1]
            new_pssm[2] = pssm_9[2]  # Anchor
            new_pssm[3] = pssm_9[3]
            new_pssm[4] = pssm_9[4]
            new_pssm[5] = pssm_9[6]
            new_pssm[6] = pssm_9[7]
            new_pssm[7] = pssm_9[8]
            new_pssm[8] = pssm_9[9]  # C-terminal anchor
        elif target_length == 10:
            # Duplicate middle position
            new_pssm[1] = pssm_9[1]
            new_pssm[2] = pssm_9[2]  # Anchor
            new_pssm[3] = pssm_9[3]
            new_pssm[4] = pssm_9[4]
            new_pssm[5] = pssm_9[5]
            new_pssm[6] = pssm_9[5]  # Duplicate
            new_pssm[7] = pssm_9[6]
            new_pssm[8] = pssm_9[7]
            new_pssm[9] = pssm_9[8]
            new_pssm[10] = pssm_9[9]  # C-terminal anchor
        elif target_length == 11:
            # Add two middle positions
            new_pssm[1] = pssm_9[1]
            new_pssm[2] = pssm_9[2]  # Anchor
            new_pssm[3] = pssm_9[3]
            new_pssm[4] = pssm_9[4]
            new_pssm[5] = pssm_9[5]
            new_pssm[6] = pssm_9[5]  # Duplicate
            new_pssm[7] = pssm_9[5]  # Duplicate
            new_pssm[8] = pssm_9[6]
            new_pssm[9] = pssm_9[7]
            new_pssm[10] = pssm_9[8]
            new_pssm[11] = pssm_9[9]  # C-terminal anchor
        else:
            # Fallback: just use 9-mer directly
            return pssm_9

        return new_pssm

    def predict_mhc2(
        self,
        sequence: str,
        allele: str = 'HLA-DRB1*01:01',
    ) -> TcellPredictionResult:
        """
        Predict MHC Class II binding peptides.

        Generates 15-mer peptides and identifies the best 9-mer binding core
        within each using allele-specific PSSMs.

        Args:
            sequence: Protein sequence.
            allele: HLA-DR allele name.

        Returns:
            TcellPredictionResult with scores for each 15-mer peptide.
        """
        sequence = self.validate_sequence(sequence)
        peptide_length = 15
        core_length = 9

        if allele not in self.mhc2_pssm:
            available = list(self.mhc2_pssm.keys())
            raise ValueError(f"Unsupported MHC-II allele: {allele}. Available: {available}")

        if len(sequence) < peptide_length:
            raise ValueError(f"Sequence length ({len(sequence)}) must be >= {peptide_length}")

        core_pssm = self.mhc2_pssm[allele]

        peptides = []
        scores = []
        positions = []
        core_positions = []
        best_cores = []

        n_peptides = len(sequence) - peptide_length + 1
        for i in range(n_peptides):
            peptide = sequence[i:i + peptide_length]
            peptides.append(peptide)
            positions.append(i + 1)

            # Find best 9-mer core within the 15-mer
            best_score = float('-inf')
            best_core_pos = 0
            best_core = ""

            for j in range(peptide_length - core_length + 1):
                core = peptide[j:j + core_length]
                core_score = score_peptide_pssm(core, core_pssm)
                if core_score > best_score:
                    best_score = core_score
                    best_core_pos = j
                    best_core = core

            scores.append(best_score)
            core_positions.append(best_core_pos + 1)  # 1-indexed within peptide
            best_cores.append(best_core)

        return TcellPredictionResult(
            scores=scores,
            peptides=peptides,
            threshold=self.MHC2_THRESHOLD,
            positions=positions,
            metadata={
                'allele': allele,
                'peptide_length': peptide_length,
                'core_length': core_length,
                'method': 'PSSM_core',
                'core_positions': core_positions,
                'best_cores': best_cores,
            }
        )

    def predict_cleavage(
        self,
        sequence: str,
    ) -> TcellPredictionResult:
        """
        Predict proteasomal cleavage sites.

        Uses a Netchop-style PSSM to score C-terminal cleavage probability
        at each position based on surrounding amino acid context.

        Args:
            sequence: Protein sequence.

        Returns:
            TcellPredictionResult with cleavage scores for each position.
        """
        sequence = self.validate_sequence(sequence)

        if len(sequence) < 7:
            raise ValueError("Sequence must be at least 7 amino acids for cleavage prediction")

        # Context window: P6-P1 ... cleavage ... P1'-P6'
        context_before = 6
        context_after = 6

        scores = []
        positions = []

        # Score cleavage at each internal position
        for i in range(context_before, len(sequence) - context_after):
            score = 0.0

            # Score positions before cleavage (P6 to P1)
            for j in range(-context_before, 0):
                pos_in_seq = i + j
                aa = sequence[pos_in_seq]
                if j in self.cleavage_pssm:
                    score += self.cleavage_pssm[j].get(aa, 0.0)

            # Score positions after cleavage (P1' to P6')
            for j in range(1, context_after + 1):
                pos_in_seq = i + j
                if pos_in_seq < len(sequence):
                    aa = sequence[pos_in_seq]
                    if j in self.cleavage_pssm:
                        score += self.cleavage_pssm[j].get(aa, 0.0)

            scores.append(score)
            positions.append(i + 1)  # 1-indexed, position of cleavage

        # Generate "peptide" fragments (show context around cleavage site)
        peptides = []
        for pos in positions:
            idx = pos - 1  # 0-indexed
            start = max(0, idx - 3)
            end = min(len(sequence), idx + 4)
            context = sequence[start:idx] + '|' + sequence[idx:end]
            peptides.append(context)

        return TcellPredictionResult(
            scores=scores,
            peptides=peptides,
            threshold=self.CLEAVAGE_THRESHOLD,
            positions=positions,
            metadata={
                'method': 'Netchop_style_PSSM',
                'context_window': f'P{context_before}-P{context_after}',
                'cleavage_symbol': '|',
            }
        )

    def predict_tap(
        self,
        sequence: str,
        peptide_length: int = 9,
    ) -> TcellPredictionResult:
        """
        Predict TAP transporter binding.

        Scores peptide affinity for the TAP transporter complex, which
        translocates peptides from cytoplasm to ER for MHC-I loading.

        Args:
            sequence: Protein sequence.
            peptide_length: Length of peptides to score (default 9).

        Returns:
            TcellPredictionResult with TAP binding scores.
        """
        sequence = self.validate_sequence(sequence)

        if len(sequence) < peptide_length:
            raise ValueError(f"Sequence length ({len(sequence)}) must be >= {peptide_length}")

        peptides = []
        scores = []
        positions = []

        n_peptides = len(sequence) - peptide_length + 1
        for i in range(n_peptides):
            peptide = sequence[i:i + peptide_length]
            peptides.append(peptide)
            positions.append(i + 1)

            # Score using TAP PSSM
            score = 0.0
            for j, aa in enumerate(peptide):
                pos = j + 1
                if pos in self.tap_pssm:
                    score += self.tap_pssm[pos].get(aa, 0.0)
                elif pos <= 9:
                    # For positions > 9, only the C-terminus matters
                    pass

            # C-terminus is position 9 in PSSM, map to actual C-term
            if peptide_length > 9:
                c_term_aa = peptide[-1]
                score += self.tap_pssm[9].get(c_term_aa, 0.0)

            scores.append(score)

        return TcellPredictionResult(
            scores=scores,
            peptides=peptides,
            threshold=self.TAP_THRESHOLD,
            positions=positions,
            metadata={
                'method': 'TAP_PSSM',
                'peptide_length': peptide_length,
            }
        )

    def predict_pipeline(
        self,
        sequence: str,
        allele: str = 'HLA-A*02:01',
        peptide_length: int = 9,
        weights: Optional[Dict[str, float]] = None,
    ) -> TcellPredictionResult:
        """
        Combined T-cell epitope prediction pipeline.

        Integrates MHC-I binding, proteasomal cleavage, and TAP transport
        predictions into a unified score. Identifies peptides that are
        likely to be naturally processed and presented.

        Args:
            sequence: Protein sequence.
            allele: HLA allele for MHC binding prediction.
            peptide_length: Length of peptides to predict.
            weights: Optional weights for combining scores. Keys:
                - 'mhc': MHC binding weight (default 0.6)
                - 'cleavage': Cleavage weight (default 0.2)
                - 'tap': TAP transport weight (default 0.2)

        Returns:
            TcellPredictionResult with combined pipeline scores.
        """
        sequence = self.validate_sequence(sequence)

        # Default weights
        if weights is None:
            weights = {'mhc': 0.6, 'cleavage': 0.2, 'tap': 0.2}

        # Run individual predictions
        mhc_result = self.predict_mhc1(sequence, allele=allele, peptide_length=peptide_length)
        cleavage_result = self.predict_cleavage(sequence)
        tap_result = self.predict_tap(sequence, peptide_length=peptide_length)

        # Normalize scores to [0, 1] range for combination
        def normalize_scores(scores: List[float]) -> np.ndarray:
            arr = np.array(scores)
            if len(arr) == 0:
                return arr
            min_val = np.min(arr)
            max_val = np.max(arr)
            if max_val - min_val < 1e-10:
                return np.zeros_like(arr) + 0.5
            return (arr - min_val) / (max_val - min_val)

        mhc_norm = normalize_scores(mhc_result.scores)
        tap_norm = normalize_scores(tap_result.scores)

        # Map cleavage scores to peptide positions
        # Cleavage at C-terminus of peptide is relevant
        cleavage_scores_mapped = []
        for pos in mhc_result.positions:
            c_term_pos = pos + peptide_length - 1
            # Find cleavage score at C-terminal position
            if c_term_pos in cleavage_result.positions:
                idx = cleavage_result.positions.index(c_term_pos)
                cleavage_scores_mapped.append(cleavage_result.scores[idx])
            else:
                cleavage_scores_mapped.append(0.0)

        cleavage_norm = normalize_scores(cleavage_scores_mapped)

        # Combine scores
        combined_scores = []
        for i in range(len(mhc_result.scores)):
            mhc_score = mhc_norm[i] if i < len(mhc_norm) else 0.0
            tap_score = tap_norm[i] if i < len(tap_norm) else 0.0
            cleavage_score = cleavage_norm[i] if i < len(cleavage_norm) else 0.0

            combined = (
                weights['mhc'] * mhc_score +
                weights['tap'] * tap_score +
                weights['cleavage'] * cleavage_score
            )
            combined_scores.append(combined)

        return TcellPredictionResult(
            scores=combined_scores,
            peptides=mhc_result.peptides,
            threshold=0.5,  # Threshold for normalized combined score
            positions=mhc_result.positions,
            ic50_values=mhc_result.ic50_values,
            percentile_ranks=mhc_result.percentile_ranks,
            metadata={
                'allele': allele,
                'peptide_length': peptide_length,
                'method': 'combined_pipeline',
                'weights': weights,
                'component_scores': {
                    'mhc': mhc_result.scores,
                    'tap': tap_result.scores,
                    'cleavage': cleavage_scores_mapped,
                },
            }
        )

    def get_available_alleles(self, mhc_class: Optional[str] = None) -> Dict[str, Any]:
        """
        Get information about available alleles.

        Args:
            mhc_class: Filter by MHC class ('MHC-I' or 'MHC-II').
                      If None, returns all alleles.

        Returns:
            Dictionary with allele information.
        """
        if mhc_class is None:
            return AVAILABLE_ALLELES
        elif mhc_class.upper() in ['MHC-I', 'MHCI', 'CLASS I', 'CLASS1', '1', 'I']:
            return {'MHC-I': AVAILABLE_ALLELES['MHC-I']}
        elif mhc_class.upper() in ['MHC-II', 'MHCII', 'CLASS II', 'CLASS2', '2', 'II']:
            return {'MHC-II': AVAILABLE_ALLELES['MHC-II']}
        else:
            raise ValueError(f"Unknown MHC class: {mhc_class}")


# =============================================================================
# Convenience Functions
# =============================================================================

def predict_mhc1(
    sequence: str,
    allele: str = 'HLA-A*02:01',
    peptide_length: int = 9,
) -> Dict[str, Any]:
    """
    Convenience function for MHC-I binding prediction.

    Args:
        sequence: Protein sequence.
        allele: HLA allele name.
        peptide_length: Peptide length (8-11).

    Returns:
        Dictionary with scores, peptides, and threshold.
    """
    predictor = TcellPredictor()
    result = predictor.predict_mhc1(sequence, allele=allele, peptide_length=peptide_length)
    return result.to_dict()


def predict_mhc2(
    sequence: str,
    allele: str = 'HLA-DRB1*01:01',
) -> Dict[str, Any]:
    """
    Convenience function for MHC-II binding prediction.

    Args:
        sequence: Protein sequence.
        allele: HLA-DR allele name.

    Returns:
        Dictionary with scores, peptides, and threshold.
    """
    predictor = TcellPredictor()
    result = predictor.predict_mhc2(sequence, allele=allele)
    return result.to_dict()


def predict_cleavage(sequence: str) -> Dict[str, Any]:
    """
    Convenience function for proteasomal cleavage prediction.

    Args:
        sequence: Protein sequence.

    Returns:
        Dictionary with scores, peptides, and threshold.
    """
    predictor = TcellPredictor()
    result = predictor.predict_cleavage(sequence)
    return result.to_dict()


def predict_tap(
    sequence: str,
    peptide_length: int = 9,
) -> Dict[str, Any]:
    """
    Convenience function for TAP transport prediction.

    Args:
        sequence: Protein sequence.
        peptide_length: Peptide length.

    Returns:
        Dictionary with scores, peptides, and threshold.
    """
    predictor = TcellPredictor()
    result = predictor.predict_tap(sequence, peptide_length=peptide_length)
    return result.to_dict()


def predict_tcell_pipeline(
    sequence: str,
    allele: str = 'HLA-A*02:01',
    peptide_length: int = 9,
) -> Dict[str, Any]:
    """
    Convenience function for combined T-cell epitope pipeline.

    Args:
        sequence: Protein sequence.
        allele: HLA allele name.
        peptide_length: Peptide length.

    Returns:
        Dictionary with combined scores, peptides, and threshold.
    """
    predictor = TcellPredictor()
    result = predictor.predict_pipeline(sequence, allele=allele, peptide_length=peptide_length)
    return result.to_dict()
