"""
Core B-cell epitope prediction engine.

Uses amino acid property scales with sliding window averaging
to predict immunogenic regions of protein sequences.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np

from .scales import get_scale, get_method_info, SCALES


@dataclass
class Epitope:
    """A predicted epitope region."""
    start: int          # 1-indexed start position
    end: int            # 1-indexed end position (inclusive)
    sequence: str       # Amino acid sequence
    score: float        # Average score for the region

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def __str__(self) -> str:
        return f"{self.start}-{self.end}: {self.sequence} (score={self.score:.3f})"


@dataclass
class Prediction:
    """Result of epitope prediction for a single sequence."""
    method: str
    sequence: str
    sequence_name: str
    window_size: int
    positions: np.ndarray       # 1-indexed center positions
    scores: np.ndarray          # Per-position scores
    threshold: float
    _epitopes: List[Epitope] = field(default_factory=list, repr=False)

    @property
    def epitopes(self) -> List[Epitope]:
        """Extract contiguous above-threshold regions as epitopes."""
        if self._epitopes:
            return self._epitopes

        self._epitopes = extract_epitopes(
            self.sequence,
            self.positions,
            self.scores,
            self.threshold,
            min_length=6
        )
        return self._epitopes

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "method": self.method,
            "sequence_name": self.sequence_name,
            "sequence": self.sequence,
            "window_size": self.window_size,
            "threshold": self.threshold,
            "positions": self.positions.tolist(),
            "scores": self.scores.tolist(),
            "epitopes": [
                {"start": e.start, "end": e.end, "sequence": e.sequence, "score": e.score}
                for e in self.epitopes
            ],
        }

    def score_at(self, position: int) -> Optional[float]:
        """Get score at a 1-indexed position, or None if outside scored region."""
        idx = np.where(self.positions == position)[0]
        if len(idx) == 0:
            return None
        return float(self.scores[idx[0]])


def extract_epitopes(
    sequence: str,
    positions: np.ndarray,
    scores: np.ndarray,
    threshold: float,
    min_length: int = 6
) -> List[Epitope]:
    """Extract contiguous above-threshold regions."""
    epitopes = []
    above = scores >= threshold

    in_epitope = False
    start_idx = 0

    for i, is_above in enumerate(above):
        if is_above and not in_epitope:
            in_epitope = True
            start_idx = i
        elif not is_above and in_epitope:
            in_epitope = False
            start_pos = int(positions[start_idx])
            end_pos = int(positions[i - 1])
            length = end_pos - start_pos + 1

            if length >= min_length:
                peptide = sequence[start_pos - 1:end_pos]
                avg_score = float(np.mean(scores[start_idx:i]))
                epitopes.append(Epitope(start_pos, end_pos, peptide, avg_score))

    # Handle epitope extending to end
    if in_epitope:
        start_pos = int(positions[start_idx])
        end_pos = int(positions[-1])
        length = end_pos - start_pos + 1

        if length >= min_length:
            peptide = sequence[start_pos - 1:end_pos]
            avg_score = float(np.mean(scores[start_idx:]))
            epitopes.append(Epitope(start_pos, end_pos, peptide, avg_score))

    return epitopes


def predict(
    sequence: str,
    method: str = "emini",
    window_size: Optional[int] = None,
    threshold: Optional[float] = None,
    sequence_name: str = "Sequence",
) -> Prediction:
    """
    Predict B-cell epitopes using a propensity scale method.

    Args:
        sequence: Amino acid sequence (one-letter codes)
        method: Prediction method name (emini, parker, chou-fasman, etc.)
        window_size: Sliding window size (default: method-specific)
        threshold: Score threshold for epitope calls (default: method-specific)
        sequence_name: Name/identifier for the sequence

    Returns:
        Prediction object with scores and epitope calls
    """
    sequence = sequence.upper().replace(" ", "").replace("\n", "")
    scale = get_scale(method)
    info = get_method_info(method)

    if window_size is None:
        window_size = info["default_window"]

    if len(sequence) < window_size:
        raise ValueError(f"Sequence length ({len(sequence)}) must be >= window size ({window_size})")

    # Vectorized scoring
    if method.lower() == "emini":
        scores, positions = _emini_score(sequence, scale, window_size)
    else:
        scores, positions = _linear_average(sequence, scale, window_size)

    # Determine threshold
    if threshold is None:
        if info["default_threshold"] is not None:
            threshold = info["default_threshold"]
        else:
            threshold = float(np.mean(scores))

    return Prediction(
        method=method,
        sequence=sequence,
        sequence_name=sequence_name,
        window_size=window_size,
        positions=positions,
        scores=scores,
        threshold=threshold,
    )


def _linear_average(
    sequence: str,
    scale: Dict[str, float],
    window_size: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute sliding window average of scale values."""
    values = np.array([scale.get(aa, 0.0) for aa in sequence])

    # Use cumsum for efficient sliding window
    cumsum = np.cumsum(np.insert(values, 0, 0))
    scores = (cumsum[window_size:] - cumsum[:-window_size]) / window_size

    # Center positions (1-indexed)
    center_offset = window_size // 2
    positions = np.arange(1 + center_offset, len(sequence) - window_size + center_offset + 2)

    return scores, positions


def _emini_score(
    sequence: str,
    scale: Dict[str, float],
    window_size: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute Emini surface accessibility score (product-based)."""
    values = np.array([scale.get(aa, 1.0) for aa in sequence])

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


def predict_all(
    sequence: str,
    methods: Optional[List[str]] = None,
    window_size: Optional[int] = None,
    sequence_name: str = "Sequence",
) -> Dict[str, Prediction]:
    """
    Run all (or specified) prediction methods on a sequence.

    Returns dict mapping method name to Prediction.
    """
    if methods is None:
        methods = list(SCALES.keys())

    return {
        method: predict(sequence, method, window_size, sequence_name=sequence_name)
        for method in methods
    }


def available_methods() -> List[str]:
    """List available prediction methods."""
    return list(SCALES.keys())
