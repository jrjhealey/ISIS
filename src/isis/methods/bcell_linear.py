"""
B-cell linear epitope prediction methods.

These methods predict linear (continuous) B-cell epitopes using amino acid
property scales with sliding window averaging. Each method uses a different
biochemical property scale derived from experimental data.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np

from .base import EpitopeMethod, MethodCategory, MethodResult
from ..scales import (
    EMINI,
    PARKER,
    CHOU_FASMAN,
    KOLASKAR_TONGAONKAR,
    KARPLUS_SCHULZ,
)


class LinearScaleMethod(EpitopeMethod):
    """
    Base class for linear scale-based prediction methods.

    These methods compute a sliding window average of amino acid
    property values to predict epitopes.
    """

    scale: Dict[str, float] = {}

    def _linear_average(
        self,
        sequence: str,
        window_size: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute sliding window average of scale values."""
        values = np.array([self.scale.get(aa, 0.0) for aa in sequence])

        # Use cumsum for efficient sliding window
        cumsum = np.cumsum(np.insert(values, 0, 0))
        scores = (cumsum[window_size:] - cumsum[:-window_size]) / window_size

        # Center positions (1-indexed)
        center_offset = window_size // 2
        positions = np.arange(
            1 + center_offset,
            len(sequence) - window_size + center_offset + 2
        )

        return scores, positions

    def predict(self, sequence: str, **kwargs) -> MethodResult:
        """
        Predict epitopes using linear sliding window averaging.

        Args:
            sequence: Amino acid sequence.
            window_size: Size of sliding window (default: method-specific).
            threshold: Score threshold for epitope calls (default: method-specific).
            min_epitope_length: Minimum epitope length (default: 6).

        Returns:
            MethodResult with scores, positions, epitopes, and metadata.
        """
        sequence = self.validate_sequence(sequence)

        window_size = kwargs.get("window_size", self.default_window)
        threshold = kwargs.get("threshold", self.default_threshold)
        min_epitope_length = kwargs.get("min_epitope_length", 6)

        if len(sequence) < window_size:
            raise ValueError(
                f"Sequence length ({len(sequence)}) must be >= window size ({window_size})"
            )

        scores, positions = self._linear_average(sequence, window_size)

        # Determine threshold if not specified
        if threshold is None:
            threshold = float(np.mean(scores))

        epitopes = self.extract_epitopes(
            sequence, positions, scores, threshold, min_epitope_length
        )

        return MethodResult(
            scores=scores,
            positions=positions,
            epitopes=epitopes,
            metadata={
                "method": self.key,
                "method_name": self.name,
                "window_size": window_size,
                "threshold": threshold,
                "sequence_length": len(sequence),
            }
        )


class EminiMethod(EpitopeMethod):
    """
    Emini surface accessibility prediction.

    Uses a product-based formula rather than simple averaging to predict
    surface-exposed regions of a protein. Higher scores indicate greater
    likelihood of surface accessibility and potential antigenicity.

    Reference:
        Emini EA et al. J Virol. 1985;55(3):836-839. PMID: 2410642
    """

    name = "Emini Surface Accessibility"
    key = "emini"
    category = MethodCategory.BCELL_LINEAR
    description = "Predicts surface-exposed regions based on accessibility."
    citation = "Emini EA et al. J Virol. 1985;55(3):836-839. PMID: 2410642"
    default_window = 6
    default_threshold = 1.0

    scale = EMINI

    def _emini_score(
        self,
        sequence: str,
        window_size: int
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Compute Emini surface accessibility score (product-based)."""
        values = np.array([self.scale.get(aa, 1.0) for aa in sequence])

        n_windows = len(sequence) - window_size + 1
        scores = np.zeros(n_windows)

        # Sliding window product
        for i in range(n_windows):
            window_vals = values[i:i + window_size]
            product = np.prod(window_vals) * (0.37 ** -6)
            scores[i] = product

        # Normalize by average
        avg = np.mean(scores)
        if avg > 0:
            scores = scores / avg

        center_offset = window_size // 2
        positions = np.arange(1 + center_offset, n_windows + center_offset + 1)

        return scores, positions

    def predict(self, sequence: str, **kwargs) -> MethodResult:
        """
        Predict surface-accessible epitopes using Emini method.

        Args:
            sequence: Amino acid sequence.
            window_size: Size of sliding window (default: 6).
            threshold: Score threshold for epitope calls (default: 1.0).
            min_epitope_length: Minimum epitope length (default: 6).

        Returns:
            MethodResult with scores, positions, epitopes, and metadata.
        """
        sequence = self.validate_sequence(sequence)

        window_size = kwargs.get("window_size", self.default_window)
        threshold = kwargs.get("threshold", self.default_threshold)
        min_epitope_length = kwargs.get("min_epitope_length", 6)

        if len(sequence) < window_size:
            raise ValueError(
                f"Sequence length ({len(sequence)}) must be >= window size ({window_size})"
            )

        scores, positions = self._emini_score(sequence, window_size)

        # Determine threshold if not specified
        if threshold is None:
            threshold = float(np.mean(scores))

        epitopes = self.extract_epitopes(
            sequence, positions, scores, threshold, min_epitope_length
        )

        return MethodResult(
            scores=scores,
            positions=positions,
            epitopes=epitopes,
            metadata={
                "method": self.key,
                "method_name": self.name,
                "window_size": window_size,
                "threshold": threshold,
                "sequence_length": len(sequence),
            }
        )


class ParkerMethod(LinearScaleMethod):
    """
    Parker hydrophilicity prediction.

    Predicts hydrophilic regions that are more likely to be exposed on the
    protein surface and thus potentially antigenic. Based on HPLC retention
    times of synthetic peptides.

    Reference:
        Parker JM et al. Biochemistry. 1986;25(19):5425-5432. PMID: 3539163
    """

    name = "Parker Hydrophilicity"
    key = "parker"
    category = MethodCategory.BCELL_LINEAR
    description = "Predicts hydrophilic regions likely to be antigenic."
    citation = "Parker JM et al. Biochemistry. 1986;25(19):5425-5432. PMID: 3539163"
    default_window = 7
    default_threshold = None  # Uses mean

    scale = PARKER


class ChouFasmanMethod(LinearScaleMethod):
    """
    Chou-Fasman beta-turn prediction.

    Predicts regions with high propensity for beta-turn secondary structure.
    Beta-turns often occur on the protein surface and are associated with
    antigenic regions.

    Reference:
        Chou PY, Fasman GD. Adv Enzymol. 1978;47:45-148. PMID: 364941
    """

    name = "Chou-Fasman Beta-Turn"
    key = "chou-fasman"
    category = MethodCategory.BCELL_LINEAR
    description = "Predicts beta-turn regions associated with epitopes."
    citation = "Chou PY, Fasman GD. Adv Enzymol. 1978;47:45-148. PMID: 364941"
    default_window = 7
    default_threshold = None  # Uses mean

    scale = CHOU_FASMAN


class KolaskarTongaonkarMethod(LinearScaleMethod):
    """
    Kolaskar-Tongaonkar antigenicity prediction.

    Semi-empirical method based on physicochemical properties of amino acids
    and their frequency of occurrence in experimentally known antigenic
    determinants. Reported accuracy of approximately 75%.

    Reference:
        Kolaskar AS, Tongaonkar PC. FEBS Lett. 1990;276(1-2):172-174. PMID: 1702393
    """

    name = "Kolaskar-Tongaonkar Antigenicity"
    key = "kolaskar-tongaonkar"
    category = MethodCategory.BCELL_LINEAR
    description = "Semi-empirical antigenicity prediction (~75% accuracy)."
    citation = "Kolaskar AS, Tongaonkar PC. FEBS Lett. 1990;276(1-2):172-174. PMID: 1702393"
    default_window = 7
    default_threshold = 1.0

    scale = KOLASKAR_TONGAONKAR


class KarplusSchulzMethod(LinearScaleMethod):
    """
    Karplus-Schulz flexibility prediction.

    Predicts flexible regions based on B-factor mobility data from protein
    crystal structures. Flexible regions are more likely to be accessible
    and immunogenic.

    Reference:
        Karplus PA, Schulz GE. Naturwissenschaften. 1985;72:212-213.
    """

    name = "Karplus-Schulz Flexibility"
    key = "karplus-schulz"
    category = MethodCategory.BCELL_LINEAR
    description = "Predicts flexible regions likely to be immunogenic."
    citation = "Karplus PA, Schulz GE. Naturwissenschaften. 1985;72:212-213."
    default_window = 7
    default_threshold = 1.0

    scale = KARPLUS_SCHULZ


# Registry of all B-cell linear methods
BCELL_LINEAR_METHODS: Dict[str, type] = {
    "emini": EminiMethod,
    "parker": ParkerMethod,
    "chou-fasman": ChouFasmanMethod,
    "kolaskar-tongaonkar": KolaskarTongaonkarMethod,
    "karplus-schulz": KarplusSchulzMethod,
}


def get_method(name: str) -> EpitopeMethod:
    """
    Get an instance of a B-cell linear prediction method by name.

    Args:
        name: Method identifier (case-insensitive, underscores converted to hyphens).

    Returns:
        Instance of the requested EpitopeMethod.

    Raises:
        ValueError: If method name is not recognized.
    """
    key = name.lower().replace("_", "-")
    if key not in BCELL_LINEAR_METHODS:
        available = list(BCELL_LINEAR_METHODS.keys())
        raise ValueError(f"Unknown method: {name}. Available: {available}")
    return BCELL_LINEAR_METHODS[key]()


def list_methods() -> Dict[str, Dict[str, Any]]:
    """
    List all available B-cell linear methods with their metadata.

    Returns:
        Dict mapping method keys to their info dicts.
    """
    return {
        key: cls().get_info()
        for key, cls in BCELL_LINEAR_METHODS.items()
    }
