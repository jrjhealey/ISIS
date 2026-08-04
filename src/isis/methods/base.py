"""
Base classes for epitope prediction methods.

This module provides the abstract base class and common data structures
for all epitope prediction methods in the ISIS library.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

import numpy as np


class MethodCategory(Enum):
    """Categories of epitope prediction methods."""
    BCELL_LINEAR = "bcell_linear"
    BCELL_CONFORMATIONAL = "bcell_conformational"
    TCELL = "tcell"
    INNATE = "innate"
    STRUCTURAL = "structural"


@dataclass
class MethodResult:
    """
    Result from an epitope prediction method.

    Attributes:
        scores: Array of prediction scores for each position.
        positions: Array of sequence positions (1-indexed) corresponding to scores.
        epitopes: List of predicted epitope regions as dicts with
                  'start', 'end', 'sequence', and 'score' keys.
        metadata: Additional method-specific information.
    """
    scores: np.ndarray
    positions: np.ndarray
    epitopes: List[Dict[str, Any]] = field(default_factory=list)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Ensure arrays are numpy arrays."""
        if not isinstance(self.scores, np.ndarray):
            self.scores = np.array(self.scores)
        if not isinstance(self.positions, np.ndarray):
            self.positions = np.array(self.positions)

    def to_dict(self) -> Dict[str, Any]:
        """Convert result to a serializable dictionary."""
        return {
            "scores": self.scores.tolist(),
            "positions": self.positions.tolist(),
            "epitopes": self.epitopes,
            "metadata": self.metadata,
        }

    @property
    def max_score(self) -> float:
        """Maximum score in the prediction."""
        return float(np.max(self.scores)) if len(self.scores) > 0 else 0.0

    @property
    def min_score(self) -> float:
        """Minimum score in the prediction."""
        return float(np.min(self.scores)) if len(self.scores) > 0 else 0.0

    @property
    def mean_score(self) -> float:
        """Mean score across all positions."""
        return float(np.mean(self.scores)) if len(self.scores) > 0 else 0.0


class EpitopeMethod(ABC):
    """
    Abstract base class for epitope prediction methods.

    All prediction methods must inherit from this class and implement
    the predict() method. Methods should also set the class attributes
    for proper categorization and metadata.

    Attributes:
        name: Human-readable method name.
        key: Lowercase identifier used in method lookups.
        category: MethodCategory enum value.
        description: Brief description of what the method predicts.
        citation: Publication reference for the method.
        default_window: Default sliding window size.
        default_threshold: Default score threshold for epitope calls.
    """

    name: str = "Base Method"
    key: str = "base"
    category: MethodCategory = MethodCategory.BCELL_LINEAR
    description: str = "Base epitope prediction method."
    citation: str = ""
    default_window: int = 7
    default_threshold: Optional[float] = None

    @abstractmethod
    def predict(self, sequence: str, **kwargs) -> MethodResult:
        """
        Predict epitopes in a protein sequence.

        Args:
            sequence: Amino acid sequence (one-letter codes).
            **kwargs: Method-specific parameters. Common options include:
                - window_size: Size of sliding window.
                - threshold: Score threshold for epitope calls.
                - min_epitope_length: Minimum length for epitope regions.

        Returns:
            MethodResult containing scores, positions, epitopes, and metadata.
        """
        pass

    def get_info(self) -> Dict[str, Any]:
        """Get method metadata as a dictionary."""
        return {
            "name": self.name,
            "key": self.key,
            "category": self.category.value,
            "description": self.description,
            "citation": self.citation,
            "default_window": self.default_window,
            "default_threshold": self.default_threshold,
        }

    def validate_sequence(self, sequence: str) -> str:
        """
        Validate and normalize an amino acid sequence.

        Args:
            sequence: Input sequence string.

        Returns:
            Normalized uppercase sequence with whitespace removed.

        Raises:
            ValueError: If sequence is empty after normalization.
        """
        normalized = sequence.upper().replace(" ", "").replace("\n", "")
        if not normalized:
            raise ValueError("Sequence cannot be empty")
        return normalized

    def extract_epitopes(
        self,
        sequence: str,
        positions: np.ndarray,
        scores: np.ndarray,
        threshold: float,
        min_length: int = 6
    ) -> List[Dict[str, Any]]:
        """
        Extract contiguous above-threshold regions as epitopes.

        Args:
            sequence: The protein sequence.
            positions: Array of 1-indexed positions.
            scores: Array of scores for each position.
            threshold: Score threshold for epitope calls.
            min_length: Minimum length for an epitope region.

        Returns:
            List of epitope dicts with 'start', 'end', 'sequence', 'score' keys.
        """
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
                    epitopes.append({
                        "start": start_pos,
                        "end": end_pos,
                        "sequence": peptide,
                        "score": avg_score
                    })

        # Handle epitope extending to end
        if in_epitope:
            start_pos = int(positions[start_idx])
            end_pos = int(positions[-1])
            length = end_pos - start_pos + 1

            if length >= min_length:
                peptide = sequence[start_pos - 1:end_pos]
                avg_score = float(np.mean(scores[start_idx:]))
                epitopes.append({
                    "start": start_pos,
                    "end": end_pos,
                    "sequence": peptide,
                    "score": avg_score
                })

        return epitopes
