"""
Conformational B-cell epitope prediction methods.

These methods REQUIRE 3D structure data (coordinates, SASA, contact numbers)
obtained from structure analysis tools like ChimeraX. They CANNOT work on
sequence alone - the structural context is essential for accurate prediction.

Methods implemented:
- DiscoTope-style: Combines propensity scores with contact numbers and SASA
- ElliPro-style: Uses protrusion index based on distance from protein centroid
- SEPPA-style: Patch-based scoring of surface-exposed residue clusters

References:
- DiscoTope: Kringelum JV et al. PLoS Comput Biol. 2012;8(12):e1002829
- ElliPro: Ponomarenko J et al. BMC Bioinformatics. 2008;9:514
- SEPPA: Sun J et al. Nucleic Acids Res. 2009;37:W612-W616
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
from numpy.typing import ArrayLike


# DiscoTope propensity scale from Kringelum et al. 2012
# Values represent log-odds of residue being in an epitope vs non-epitope
DISCOTOPE_PROPENSITY = {
    'A': -0.02, 'C': -0.41, 'D': 0.25, 'E': 0.21, 'F': -0.44,
    'G': -0.07, 'H': 0.10, 'I': -0.50, 'K': 0.13, 'L': -0.55,
    'M': -0.34, 'N': 0.15, 'P': 0.02, 'Q': 0.07, 'R': 0.09,
    'S': 0.10, 'T': 0.04, 'V': -0.45, 'W': -0.20, 'Y': 0.04
}


@dataclass
class ConformationalEpitope:
    """A predicted conformational epitope region.

    Unlike linear epitopes, conformational epitopes may consist of
    residues that are distant in sequence but close in 3D space.
    """
    residue_indices: List[int]  # 0-indexed residue positions
    residues: str               # Amino acid sequence (may be non-contiguous)
    score: float                # Average score for the epitope
    centroid: Optional[np.ndarray] = None  # 3D centroid of the epitope

    @property
    def size(self) -> int:
        return len(self.residue_indices)

    def __str__(self) -> str:
        if len(self.residue_indices) <= 6:
            positions = ",".join(str(i+1) for i in self.residue_indices)
        else:
            positions = f"{self.residue_indices[0]+1}..{self.residue_indices[-1]+1}"
        return f"Epitope({positions}): {self.residues[:20]}{'...' if len(self.residues) > 20 else ''} (score={self.score:.3f})"


@dataclass
class ConformationalPrediction:
    """Result of structure-based epitope prediction."""
    method: str
    sequence: str
    scores: np.ndarray          # Per-residue scores (0-indexed)
    threshold: float
    epitopes: List[ConformationalEpitope] = field(default_factory=list)
    metadata: Dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            "method": self.method,
            "sequence": self.sequence,
            "scores": self.scores.tolist(),
            "threshold": self.threshold,
            "epitopes": [
                {
                    "residue_indices": e.residue_indices,
                    "residues": e.residues,
                    "score": e.score,
                    "centroid": e.centroid.tolist() if e.centroid is not None else None
                }
                for e in self.epitopes
            ],
            "metadata": self.metadata
        }


# =============================================================================
# Helper Functions
# =============================================================================

def calculate_protrusion_index(coords: ArrayLike) -> np.ndarray:
    """
    Calculate the protrusion index for each residue based on CA coordinates.

    The protrusion index (PI) measures how much a residue protrudes from the
    protein surface. It is calculated as the normalized distance from the
    protein centroid:

        PI = distance_from_centroid / max_distance

    Higher PI values indicate more protruding residues, which are more likely
    to be part of conformational epitopes.

    Args:
        coords: Nx3 array of CA (alpha-carbon) coordinates.
                Each row is [x, y, z] for one residue.

    Returns:
        Array of protrusion index values (0-1) for each residue.

    Raises:
        ValueError: If coords is empty or not Nx3 shape.

    Example:
        >>> coords = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]])
        >>> pi = calculate_protrusion_index(coords)
        >>> pi[-1] > pi[1]  # Last residue protrudes more than middle
        True
    """
    coords = np.asarray(coords, dtype=np.float64)

    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"coords must be Nx3 array, got shape {coords.shape}")

    if len(coords) == 0:
        return np.array([])

    if len(coords) == 1:
        return np.array([0.0])

    # Rows may be all-NaN for residues with no resolved CA position (e.g.
    # missing density in the model) - exclude them from the centroid/max
    # distance calculation so a few missing residues don't poison every
    # other residue's protrusion index. Their own output stays NaN.
    valid = ~np.isnan(coords).any(axis=1)
    if not np.any(valid):
        return np.full(len(coords), np.nan)

    # Calculate centroid over resolved residues only
    centroid = np.mean(coords[valid], axis=0)

    # Calculate distances from centroid (NaN rows propagate to NaN distances)
    distances = np.linalg.norm(coords - centroid, axis=1)

    # Normalize by maximum distance among resolved residues
    max_distance = np.max(distances[valid])
    if max_distance > 0:
        protrusion_index = distances / max_distance
    else:
        protrusion_index = np.where(valid, 0.0, np.nan)

    return protrusion_index


def find_surface_patches(
    coords: ArrayLike,
    sasa: ArrayLike,
    threshold: float = 20.0,
    distance_cutoff: float = 10.0,
    min_patch_size: int = 3
) -> List[List[int]]:
    """
    Find clusters of nearby surface-exposed residues (patches).

    Surface patches are identified by:
    1. Selecting residues with SASA above a threshold
    2. Clustering them by spatial proximity (distance_cutoff)
    3. Filtering patches by minimum size

    This implements a simplified version of the patch detection used in
    SEPPA and similar methods.

    Args:
        coords: Nx3 array of CA coordinates.
        sasa: Array of SASA values per residue (in Angstrom^2).
        threshold: Minimum SASA to consider a residue surface-exposed.
                   Default 20.0 A^2 (roughly 10% relative accessibility).
        distance_cutoff: Maximum distance (Angstroms) between residues
                         in the same patch. Default 10.0 A.
        min_patch_size: Minimum number of residues to form a valid patch.

    Returns:
        List of patches, where each patch is a list of 0-indexed residue
        indices belonging to that cluster.

    Example:
        >>> coords = np.array([[0,0,0], [1,1,0], [10,10,10], [11,11,10]])
        >>> sasa = np.array([50.0, 45.0, 60.0, 55.0])  # All surface-exposed
        >>> patches = find_surface_patches(coords, sasa, threshold=20.0)
        >>> len(patches)  # Should find 2 patches
        2
    """
    coords = np.asarray(coords, dtype=np.float64)
    sasa = np.asarray(sasa, dtype=np.float64)

    if len(coords) != len(sasa):
        raise ValueError(f"coords ({len(coords)}) and sasa ({len(sasa)}) must have same length")

    if len(coords) == 0:
        return []

    # Find surface-exposed residues
    surface_mask = sasa >= threshold
    surface_indices = np.where(surface_mask)[0]

    if len(surface_indices) < min_patch_size:
        return []

    # Get coordinates of surface residues
    surface_coords = coords[surface_indices]

    # Build patches using simple connectivity-based clustering
    # (Union-Find / single-linkage clustering)
    n_surface = len(surface_indices)
    parent = list(range(n_surface))

    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]

    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    # Calculate pairwise distances and cluster nearby residues
    for i in range(n_surface):
        for j in range(i + 1, n_surface):
            dist = np.linalg.norm(surface_coords[i] - surface_coords[j])
            if dist <= distance_cutoff:
                union(i, j)

    # Group by cluster
    clusters: Dict[int, List[int]] = {}
    for i in range(n_surface):
        root = find(i)
        if root not in clusters:
            clusters[root] = []
        clusters[root].append(int(surface_indices[i]))

    # Filter by minimum size and sort by size (largest first)
    patches = [
        sorted(indices) for indices in clusters.values()
        if len(indices) >= min_patch_size
    ]
    patches.sort(key=len, reverse=True)

    return patches


# =============================================================================
# Base Predictor Class
# =============================================================================

class ConformationalPredictor(ABC):
    """
    Abstract base class for conformational B-cell epitope predictors.

    All conformational predictors require 3D structure data and cannot
    operate on sequence alone. Subclasses must implement the
    `predict_with_structure` method.
    """

    name: str = "base"
    description: str = "Base conformational epitope predictor"
    default_threshold: float = 0.5

    @abstractmethod
    def predict_with_structure(
        self,
        sequence: str,
        sasa_values: ArrayLike,
        contact_numbers: ArrayLike,
        coordinates: ArrayLike,
        threshold: Optional[float] = None
    ) -> Dict:
        """
        Predict conformational B-cell epitopes using structure data.

        This method REQUIRES 3D structure information. It cannot work on
        sequence alone - pass structural features from ChimeraX or similar.

        Args:
            sequence: Amino acid sequence (one-letter codes, must match
                      the length of structural arrays).
            sasa_values: Per-residue solvent accessible surface area (SASA)
                         in Angstrom^2. Array of length N.
            contact_numbers: Number of neighboring atoms/residues within
                            a cutoff distance. Array of length N.
            coordinates: Nx3 array of CA coordinates for each residue.
            threshold: Score threshold for epitope calling. If None,
                      uses method-specific default.

        Returns:
            Dictionary containing:
                - 'scores': numpy array of per-residue prediction scores
                - 'epitopes': list of ConformationalEpitope objects
                - 'threshold': the threshold value used
                - 'method': name of the prediction method

        Raises:
            ValueError: If array lengths don't match sequence length.
        """
        pass

    def _validate_inputs(
        self,
        sequence: str,
        sasa_values: ArrayLike,
        contact_numbers: ArrayLike,
        coordinates: ArrayLike
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Validate and convert inputs to numpy arrays."""
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        n_residues = len(sequence)

        sasa = np.asarray(sasa_values, dtype=np.float64)
        contacts = np.asarray(contact_numbers, dtype=np.float64)
        coords = np.asarray(coordinates, dtype=np.float64)

        if len(sasa) != n_residues:
            raise ValueError(
                f"SASA array length ({len(sasa)}) must match "
                f"sequence length ({n_residues})"
            )

        if len(contacts) != n_residues:
            raise ValueError(
                f"Contact numbers length ({len(contacts)}) must match "
                f"sequence length ({n_residues})"
            )

        if coords.shape[0] != n_residues:
            raise ValueError(
                f"Coordinates array rows ({coords.shape[0]}) must match "
                f"sequence length ({n_residues})"
            )

        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError(
                f"Coordinates must be Nx3 array, got shape {coords.shape}"
            )

        return sasa, contacts, coords

    def _extract_epitopes_from_scores(
        self,
        sequence: str,
        scores: np.ndarray,
        threshold: float,
        coordinates: Optional[np.ndarray] = None,
        min_size: int = 3,
        max_gap: int = 2
    ) -> List[ConformationalEpitope]:
        """
        Extract epitopes from per-residue scores.

        Identifies contiguous regions above threshold, allowing small gaps.

        Args:
            sequence: Amino acid sequence
            scores: Per-residue prediction scores
            threshold: Score threshold for epitope calling
            coordinates: Optional coordinates for centroid calculation
            min_size: Minimum epitope size
            max_gap: Maximum allowed gap between above-threshold residues
        """
        epitopes = []
        n = len(scores)

        above_threshold = scores >= threshold

        # Find contiguous regions with allowed gaps
        i = 0
        while i < n:
            if not above_threshold[i]:
                i += 1
                continue

            # Start of potential epitope
            start = i
            end = i
            gap_count = 0

            while end < n - 1:
                if above_threshold[end + 1]:
                    end += 1
                    gap_count = 0
                elif gap_count < max_gap:
                    end += 1
                    gap_count += 1
                else:
                    # End of region (excluding trailing gaps)
                    while end > start and not above_threshold[end]:
                        end -= 1
                    break

            # Handle trailing gaps
            while end > start and not above_threshold[end]:
                end -= 1

            # Create epitope if large enough
            indices = list(range(start, end + 1))
            actual_above = sum(above_threshold[start:end+1])

            if actual_above >= min_size:
                residues = sequence[start:end + 1]
                avg_score = float(np.mean(scores[start:end + 1]))

                centroid = None
                if coordinates is not None:
                    centroid = np.mean(coordinates[start:end + 1], axis=0)

                epitopes.append(ConformationalEpitope(
                    residue_indices=indices,
                    residues=residues,
                    score=avg_score,
                    centroid=centroid
                ))

            i = end + 1

        return epitopes


# =============================================================================
# DiscoTope-style Predictor
# =============================================================================

class DiscoTopePredictor(ConformationalPredictor):
    """
    DiscoTope-style conformational epitope predictor.

    Combines three features to predict epitope residues:
    1. Amino acid propensity (log-odds from epitope vs non-epitope)
    2. Contact number (negative contribution - buried residues less likely)
    3. Solvent accessible surface area (positive contribution)

    The final score is:
        score = propensity + alpha * contact_number + beta * log(SASA + 1)

    Where alpha = -0.02 and beta = 0.35 (from DiscoTope 2.0 paper).

    Reference:
        Kringelum JV et al. (2012) Reliable B Cell Epitope Predictions:
        Impacts of Method Development and Improved Benchmarking.
        PLoS Comput Biol 8(12): e1002829

    IMPORTANT: This method REQUIRES 3D structure data. It cannot predict
    epitopes from sequence alone.
    """

    name = "discotope"
    description = "DiscoTope-style scoring combining propensity, contacts, and SASA"
    default_threshold = -7.7  # DiscoTope 2.0 default threshold

    # Coefficients from DiscoTope 2.0
    ALPHA = -0.02  # Contact number weight (negative - more contacts = lower score)
    BETA = 0.35    # SASA weight (positive - more exposed = higher score)

    def __init__(
        self,
        alpha: float = ALPHA,
        beta: float = BETA,
        propensity_scale: Optional[Dict[str, float]] = None
    ):
        """
        Initialize DiscoTope predictor.

        Args:
            alpha: Weight for contact number term (default -0.02)
            beta: Weight for log(SASA) term (default 0.35)
            propensity_scale: Custom propensity scale, or None to use default
        """
        self.alpha = alpha
        self.beta = beta
        self.propensity_scale = propensity_scale or DISCOTOPE_PROPENSITY

    def predict_with_structure(
        self,
        sequence: str,
        sasa_values: ArrayLike,
        contact_numbers: ArrayLike,
        coordinates: ArrayLike,
        threshold: Optional[float] = None
    ) -> Dict:
        """
        Predict conformational epitopes using DiscoTope-style scoring.

        REQUIRES 3D STRUCTURE DATA - cannot work on sequence alone.

        The score for each residue is calculated as:
            score = propensity + alpha * contact_number + beta * log(SASA + 1)

        Args:
            sequence: Amino acid sequence (one-letter codes)
            sasa_values: Per-residue SASA in Angstrom^2 (from ChimeraX)
            contact_numbers: Per-residue contact counts (from ChimeraX)
            coordinates: Nx3 array of CA coordinates (from ChimeraX)
            threshold: Score threshold (default -7.7 from DiscoTope 2.0)

        Returns:
            Dictionary with 'scores', 'epitopes', 'threshold', 'method'
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        sasa, contacts, coords = self._validate_inputs(
            sequence, sasa_values, contact_numbers, coordinates
        )

        if threshold is None:
            threshold = self.default_threshold

        # Calculate propensity scores
        propensity_scores = np.array([
            self.propensity_scale.get(aa, 0.0) for aa in sequence
        ])

        # Calculate combined score
        # Score = propensity + alpha * contacts + beta * log(SASA + 1)
        scores = (
            propensity_scores +
            self.alpha * contacts +
            self.beta * np.log(sasa + 1)
        )

        # Extract epitopes
        epitopes = self._extract_epitopes_from_scores(
            sequence, scores, threshold, coords
        )

        return {
            'scores': scores,
            'epitopes': epitopes,
            'threshold': threshold,
            'method': self.name,
            'propensity_scores': propensity_scores,
            'contact_contribution': self.alpha * contacts,
            'sasa_contribution': self.beta * np.log(sasa + 1)
        }


# =============================================================================
# ElliPro-style Predictor
# =============================================================================

class ElliProPredictor(ConformationalPredictor):
    """
    ElliPro-style conformational epitope predictor based on protrusion.

    Uses the protrusion index (PI) as the primary feature. The PI measures
    how far each residue protrudes from the protein surface:

        PI = distance_from_centroid / max_distance_from_centroid

    Residues with higher PI values are more exposed and more likely to be
    part of conformational epitopes.

    The method optionally incorporates SASA as a secondary filter to ensure
    predicted epitopes are actually solvent-accessible.

    Reference:
        Ponomarenko J et al. (2008) ElliPro: a new structure-based tool
        for the prediction of antibody epitopes.
        BMC Bioinformatics 9:514

    IMPORTANT: This method REQUIRES 3D structure data (coordinates).
    It cannot predict epitopes from sequence alone.
    """

    name = "ellipro"
    description = "ElliPro-style protrusion index based prediction"
    default_threshold = 0.5  # Protrusion index is 0-1, default = top 50%

    def __init__(
        self,
        min_sasa: float = 5.0,
        combine_with_sasa: bool = True
    ):
        """
        Initialize ElliPro predictor.

        Args:
            min_sasa: Minimum SASA (A^2) to consider a residue surface-exposed.
            combine_with_sasa: If True, multiply PI by normalized SASA.
        """
        self.min_sasa = min_sasa
        self.combine_with_sasa = combine_with_sasa

    def predict_with_structure(
        self,
        sequence: str,
        sasa_values: ArrayLike,
        contact_numbers: ArrayLike,
        coordinates: ArrayLike,
        threshold: Optional[float] = None
    ) -> Dict:
        """
        Predict conformational epitopes using protrusion index.

        REQUIRES 3D STRUCTURE DATA - cannot work on sequence alone.

        The protrusion index (PI) is calculated as:
            PI = distance_from_centroid / max_distance

        Optionally combined with normalized SASA for more accurate results.

        Args:
            sequence: Amino acid sequence (one-letter codes)
            sasa_values: Per-residue SASA in Angstrom^2 (from ChimeraX)
            contact_numbers: Per-residue contact counts (not used but required
                            for interface consistency)
            coordinates: Nx3 array of CA coordinates (from ChimeraX)
            threshold: PI threshold (default 0.5, range 0-1)

        Returns:
            Dictionary with 'scores', 'epitopes', 'threshold', 'method'
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        sasa, contacts, coords = self._validate_inputs(
            sequence, sasa_values, contact_numbers, coordinates
        )

        if threshold is None:
            threshold = self.default_threshold

        # Calculate protrusion index
        protrusion_index = calculate_protrusion_index(coords)

        # Optionally combine with SASA
        if self.combine_with_sasa:
            # Normalize SASA (0-1 scale)
            max_sasa = np.max(sasa)
            if max_sasa > 0:
                normalized_sasa = sasa / max_sasa
            else:
                normalized_sasa = np.zeros_like(sasa)

            # Combined score = PI * normalized_SASA
            # This ensures high scores only for residues that are both
            # protruding AND solvent-accessible
            scores = protrusion_index * normalized_sasa

            # Apply minimum SASA filter
            scores[sasa < self.min_sasa] = 0.0
        else:
            scores = protrusion_index

        # Extract epitopes
        epitopes = self._extract_epitopes_from_scores(
            sequence, scores, threshold, coords
        )

        return {
            'scores': scores,
            'epitopes': epitopes,
            'threshold': threshold,
            'method': self.name,
            'protrusion_index': protrusion_index,
            'sasa_filter_applied': self.combine_with_sasa
        }


# =============================================================================
# SEPPA-style Predictor
# =============================================================================

class SEPPAPredictor(ConformationalPredictor):
    """
    SEPPA-style patch-based conformational epitope predictor.

    This method predicts epitopes based on surface patches (clusters of
    nearby surface-exposed residues). The approach:

    1. Identify surface-exposed residues (SASA above threshold)
    2. Cluster nearby surface residues into patches
    3. Score each patch by average propensity
    4. Map patch scores back to individual residues

    Residues belonging to high-scoring patches are predicted as epitopes.

    Reference:
        Sun J et al. (2009) SEPPA: a computational server for spatial
        epitope prediction of protein antigens.
        Nucleic Acids Res 37:W612-W616

    IMPORTANT: This method REQUIRES 3D structure data (coordinates and SASA).
    It cannot predict epitopes from sequence alone.
    """

    name = "seppa"
    description = "SEPPA-style patch-based surface epitope prediction"
    default_threshold = 0.0  # Average propensity; positive = likely epitope

    def __init__(
        self,
        sasa_cutoff: float = 20.0,
        distance_cutoff: float = 10.0,
        min_patch_size: int = 3,
        propensity_scale: Optional[Dict[str, float]] = None
    ):
        """
        Initialize SEPPA predictor.

        Args:
            sasa_cutoff: Minimum SASA (A^2) to consider a residue surface-exposed.
            distance_cutoff: Maximum CA-CA distance (A) for residues in same patch.
            min_patch_size: Minimum number of residues to form a valid patch.
            propensity_scale: Custom propensity scale, or None to use DiscoTope scale.
        """
        self.sasa_cutoff = sasa_cutoff
        self.distance_cutoff = distance_cutoff
        self.min_patch_size = min_patch_size
        self.propensity_scale = propensity_scale or DISCOTOPE_PROPENSITY

    def predict_with_structure(
        self,
        sequence: str,
        sasa_values: ArrayLike,
        contact_numbers: ArrayLike,
        coordinates: ArrayLike,
        threshold: Optional[float] = None
    ) -> Dict:
        """
        Predict conformational epitopes using patch-based scoring.

        REQUIRES 3D STRUCTURE DATA - cannot work on sequence alone.

        1. Finds surface patches (clusters of exposed residues)
        2. Scores patches by average propensity
        3. Maps patch scores to individual residues

        Args:
            sequence: Amino acid sequence (one-letter codes)
            sasa_values: Per-residue SASA in Angstrom^2 (from ChimeraX)
            contact_numbers: Per-residue contact counts (not used directly)
            coordinates: Nx3 array of CA coordinates (from ChimeraX)
            threshold: Propensity threshold (default 0.0)

        Returns:
            Dictionary with 'scores', 'epitopes', 'threshold', 'method', 'patches'
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        sasa, contacts, coords = self._validate_inputs(
            sequence, sasa_values, contact_numbers, coordinates
        )

        if threshold is None:
            threshold = self.default_threshold

        n_residues = len(sequence)

        # Find surface patches
        patches = find_surface_patches(
            coords, sasa,
            threshold=self.sasa_cutoff,
            distance_cutoff=self.distance_cutoff,
            min_patch_size=self.min_patch_size
        )

        # Score each patch by average propensity
        patch_scores = []
        for patch_indices in patches:
            propensities = [
                self.propensity_scale.get(sequence[i], 0.0)
                for i in patch_indices
            ]
            patch_score = np.mean(propensities)
            patch_scores.append(patch_score)

        # Map patch scores to residues
        # Each residue gets the score of its highest-scoring patch
        residue_scores = np.full(n_residues, np.nan)
        residue_patch_membership = [[] for _ in range(n_residues)]

        for patch_idx, (patch_indices, patch_score) in enumerate(zip(patches, patch_scores)):
            for res_idx in patch_indices:
                residue_patch_membership[res_idx].append(patch_idx)
                if np.isnan(residue_scores[res_idx]):
                    residue_scores[res_idx] = patch_score
                else:
                    residue_scores[res_idx] = max(residue_scores[res_idx], patch_score)

        # Residues not in any patch get minimum score
        min_score = min(patch_scores) if patch_scores else -1.0
        residue_scores = np.where(np.isnan(residue_scores), min_score - 0.5, residue_scores)

        # Create epitopes from high-scoring patches
        epitopes = []
        for patch_idx, (patch_indices, patch_score) in enumerate(zip(patches, patch_scores)):
            if patch_score >= threshold:
                residues = "".join(sequence[i] for i in patch_indices)
                centroid = np.mean(coords[patch_indices], axis=0)

                epitopes.append(ConformationalEpitope(
                    residue_indices=patch_indices,
                    residues=residues,
                    score=patch_score,
                    centroid=centroid
                ))

        # Sort epitopes by score (highest first)
        epitopes.sort(key=lambda e: e.score, reverse=True)

        return {
            'scores': residue_scores,
            'epitopes': epitopes,
            'threshold': threshold,
            'method': self.name,
            'patches': [
                {
                    'indices': patch_indices,
                    'score': score,
                    'residues': "".join(sequence[i] for i in patch_indices)
                }
                for patch_indices, score in zip(patches, patch_scores)
            ],
            'n_patches': len(patches)
        }


# =============================================================================
# Convenience Functions
# =============================================================================

def predict_discotope(
    sequence: str,
    sasa_values: ArrayLike,
    contact_numbers: ArrayLike,
    coordinates: ArrayLike,
    threshold: Optional[float] = None
) -> Dict:
    """
    Convenience function for DiscoTope-style prediction.

    REQUIRES 3D STRUCTURE DATA - cannot work on sequence alone.

    See DiscoTopePredictor for full documentation.
    """
    predictor = DiscoTopePredictor()
    return predictor.predict_with_structure(
        sequence, sasa_values, contact_numbers, coordinates, threshold
    )


def predict_ellipro(
    sequence: str,
    sasa_values: ArrayLike,
    contact_numbers: ArrayLike,
    coordinates: ArrayLike,
    threshold: Optional[float] = None
) -> Dict:
    """
    Convenience function for ElliPro-style prediction.

    REQUIRES 3D STRUCTURE DATA - cannot work on sequence alone.

    See ElliProPredictor for full documentation.
    """
    predictor = ElliProPredictor()
    return predictor.predict_with_structure(
        sequence, sasa_values, contact_numbers, coordinates, threshold
    )


def predict_seppa(
    sequence: str,
    sasa_values: ArrayLike,
    contact_numbers: ArrayLike,
    coordinates: ArrayLike,
    threshold: Optional[float] = None
) -> Dict:
    """
    Convenience function for SEPPA-style prediction.

    REQUIRES 3D STRUCTURE DATA - cannot work on sequence alone.

    See SEPPAPredictor for full documentation.
    """
    predictor = SEPPAPredictor()
    return predictor.predict_with_structure(
        sequence, sasa_values, contact_numbers, coordinates, threshold
    )


# =============================================================================
# Combined Prediction
# =============================================================================

def predict_all_conformational(
    sequence: str,
    sasa_values: ArrayLike,
    contact_numbers: ArrayLike,
    coordinates: ArrayLike
) -> Dict[str, Dict]:
    """
    Run all conformational epitope prediction methods.

    REQUIRES 3D STRUCTURE DATA - cannot work on sequence alone.

    Args:
        sequence: Amino acid sequence
        sasa_values: Per-residue SASA values
        contact_numbers: Per-residue contact counts
        coordinates: Nx3 array of CA coordinates

    Returns:
        Dictionary mapping method name to prediction results.
    """
    return {
        'discotope': predict_discotope(sequence, sasa_values, contact_numbers, coordinates),
        'ellipro': predict_ellipro(sequence, sasa_values, contact_numbers, coordinates),
        'seppa': predict_seppa(sequence, sasa_values, contact_numbers, coordinates)
    }
