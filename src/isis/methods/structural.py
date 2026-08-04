"""
Structural analysis utilities for epitope prediction.

These utilities process structural data extracted from ChimeraX.
They do NOT call ChimeraX directly - the ChimeraX plugin calls
these functions with pre-extracted coordinate and property data.

All functions are pure numpy operations for efficient computation.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray


# Feature names for structural epitope prediction
STRUCTURAL_FEATURES = ['sasa', 'contacts', 'bfactors', 'depth', 'protrusion']

# Default weights for composite scoring
# Positive weights favor epitope prediction, negative disfavor
DEFAULT_WEIGHTS = {
    'sasa': 0.35,       # High surface accessibility -> epitope
    'contacts': -0.25,  # Many contacts -> buried -> not epitope
    'bfactors': 0.25,   # High flexibility -> epitope
    'depth': -0.15,     # Deep residues -> not epitope
}

# Secondary structure accessibility propensities
# Based on empirical observations of epitope locations
SS_PROPENSITY = {
    'C': 1.0,   # Coil/loop - most accessible
    'T': 0.95,  # Turn - highly accessible
    'S': 0.90,  # Bend - accessible
    'G': 0.80,  # 3-10 helix - moderate
    'H': 0.60,  # Alpha helix - less accessible
    'E': 0.55,  # Beta strand - less accessible
    'B': 0.50,  # Beta bridge - less accessible
    'I': 0.65,  # Pi helix - moderate
    ' ': 1.0,   # Undefined treated as coil
    '-': 1.0,   # Gap treated as coil
}


def calculate_residue_depth(
    ca_coords: NDArray[np.floating],
    surface_points: NDArray[np.floating]
) -> NDArray[np.floating]:
    """
    Calculate minimum distance from each CA to nearest surface point.

    Residue depth measures how buried a residue is within the protein.
    Surface residues have low depth values, core residues have high values.

    Args:
        ca_coords: Nx3 array of CA (alpha carbon) coordinates
        surface_points: Mx3 array of molecular surface points

    Returns:
        depths: N array of depth values in Angstroms

    Raises:
        ValueError: If input arrays have invalid shapes

    Example:
        >>> ca_coords = np.array([[0, 0, 0], [1, 1, 1], [5, 5, 5]])
        >>> surface = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
        >>> depths = calculate_residue_depth(ca_coords, surface)
    """
    ca_coords = np.asarray(ca_coords, dtype=np.float64)
    surface_points = np.asarray(surface_points, dtype=np.float64)

    if ca_coords.ndim != 2 or ca_coords.shape[1] != 3:
        raise ValueError(f"ca_coords must be Nx3 array, got shape {ca_coords.shape}")
    if surface_points.ndim != 2 or surface_points.shape[1] != 3:
        raise ValueError(
            f"surface_points must be Mx3 array, got shape {surface_points.shape}"
        )

    n_residues = ca_coords.shape[0]
    depths = np.zeros(n_residues, dtype=np.float64)

    # Compute pairwise distances efficiently using broadcasting
    # For large arrays, process in chunks to avoid memory issues
    chunk_size = 1000
    n_surface = surface_points.shape[0]

    for i in range(n_residues):
        min_dist = np.inf
        for j in range(0, n_surface, chunk_size):
            end_j = min(j + chunk_size, n_surface)
            # Distance from CA i to surface points j:end_j
            diff = surface_points[j:end_j] - ca_coords[i]
            dists = np.sqrt(np.sum(diff ** 2, axis=1))
            chunk_min = np.min(dists)
            if chunk_min < min_dist:
                min_dist = chunk_min
        depths[i] = min_dist

    return depths


def calculate_contact_number(
    ca_coords: NDArray[np.floating],
    cutoff: float = 8.0
) -> NDArray[np.integer]:
    """
    Count neighboring residues within cutoff distance.

    Contact number is a measure of local packing. Buried residues
    typically have more contacts than surface-exposed residues.

    Args:
        ca_coords: Nx3 array of CA coordinates
        cutoff: distance threshold in Angstroms (default 8.0)

    Returns:
        contacts: N array of contact counts (excluding self)

    Raises:
        ValueError: If cutoff is not positive or coords have invalid shape

    Example:
        >>> coords = np.array([[0, 0, 0], [3, 0, 0], [10, 0, 0]])
        >>> contacts = calculate_contact_number(coords, cutoff=5.0)
        >>> # First two residues contact each other, third is isolated
    """
    ca_coords = np.asarray(ca_coords, dtype=np.float64)

    if cutoff <= 0:
        raise ValueError(f"cutoff must be positive, got {cutoff}")
    if ca_coords.ndim != 2 or ca_coords.shape[1] != 3:
        raise ValueError(f"ca_coords must be Nx3 array, got shape {ca_coords.shape}")

    n_residues = ca_coords.shape[0]

    # Compute pairwise distance matrix
    # diff[i, j, k] = ca_coords[i, k] - ca_coords[j, k]
    diff = ca_coords[:, np.newaxis, :] - ca_coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=2))

    # Count neighbors within cutoff (excluding self - diagonal)
    contact_matrix = (dist_matrix < cutoff) & (dist_matrix > 0)
    contacts = np.sum(contact_matrix, axis=1)

    return contacts.astype(np.int32)


def calculate_local_density(
    ca_coords: NDArray[np.floating],
    radius: float = 10.0
) -> NDArray[np.floating]:
    """
    Calculate local packing density around each residue.

    Density is computed as the number of CA atoms within the specified
    radius, divided by the volume of the sphere (4/3 * pi * r^3).

    Args:
        ca_coords: Nx3 array of CA coordinates
        radius: sphere radius in Angstroms (default 10.0)

    Returns:
        densities: N array of density values (atoms per cubic Angstrom)

    Raises:
        ValueError: If radius is not positive or coords have invalid shape

    Example:
        >>> coords = np.random.randn(100, 3) * 5  # Random cloud
        >>> densities = calculate_local_density(coords, radius=8.0)
    """
    ca_coords = np.asarray(ca_coords, dtype=np.float64)

    if radius <= 0:
        raise ValueError(f"radius must be positive, got {radius}")
    if ca_coords.ndim != 2 or ca_coords.shape[1] != 3:
        raise ValueError(f"ca_coords must be Nx3 array, got shape {ca_coords.shape}")

    # Get contact counts within radius
    contacts = calculate_contact_number(ca_coords, cutoff=radius)

    # Volume of sphere
    sphere_volume = (4.0 / 3.0) * np.pi * (radius ** 3)

    # Density = (neighbors + self) / volume
    densities = (contacts.astype(np.float64) + 1) / sphere_volume

    return densities


def score_ss_accessibility(ss_codes: str) -> NDArray[np.floating]:
    """
    Score accessibility based on secondary structure.

    Loops and turns are more likely to be surface-exposed and thus
    more likely to be epitopes. Alpha helices and beta sheets are
    often buried in the protein core.

    Args:
        ss_codes: string of secondary structure codes
            H = alpha helix
            E = beta strand
            C = coil/loop
            T = turn
            S = bend
            G = 3-10 helix
            B = beta bridge
            I = pi helix

    Returns:
        scores: accessibility propensity per position (0-1 scale)

    Example:
        >>> ss = "CCHHHHHHCCEEEECCC"
        >>> scores = score_ss_accessibility(ss)
        >>> # Coil positions will have higher scores
    """
    if not isinstance(ss_codes, str):
        raise TypeError(f"ss_codes must be a string, got {type(ss_codes)}")

    n_residues = len(ss_codes)
    scores = np.zeros(n_residues, dtype=np.float64)

    for i, code in enumerate(ss_codes):
        # Use propensity if known, default to coil propensity
        scores[i] = SS_PROPENSITY.get(code.upper(), SS_PROPENSITY['C'])

    return scores


def normalize_bfactors(
    bfactors: NDArray[np.floating]
) -> NDArray[np.floating]:
    """
    Normalize B-factors to 0-1 range using z-score method.

    B-factors (temperature factors) measure atomic mobility/flexibility.
    Higher values indicate more flexible regions, which are often found
    in epitopes. This function normalizes using z-scores, then scales
    to 0-1 range using a sigmoid-like transformation.

    Args:
        bfactors: N array of B-factor values

    Returns:
        normalized: N array of normalized values in [0, 1] range
            Higher values = more flexible = more likely epitope

    Raises:
        ValueError: If input array is empty

    Example:
        >>> bfactors = np.array([10.0, 25.0, 50.0, 30.0, 15.0])
        >>> norm = normalize_bfactors(bfactors)
        >>> # Highest B-factor (50) will have value closest to 1
    """
    bfactors = np.asarray(bfactors, dtype=np.float64).flatten()

    if len(bfactors) == 0:
        raise ValueError("bfactors array cannot be empty")

    # Handle constant B-factors
    std = np.std(bfactors)
    if std < 1e-10:
        # All B-factors are the same - return 0.5 for all
        return np.full_like(bfactors, 0.5)

    # Z-score normalization
    mean = np.mean(bfactors)
    z_scores = (bfactors - mean) / std

    # Transform to 0-1 using sigmoid
    # tanh maps (-inf, inf) to (-1, 1), then scale to (0, 1)
    normalized = (np.tanh(z_scores / 2) + 1) / 2

    return normalized


def detect_surface_patches(
    ca_coords: NDArray[np.floating],
    sasa_values: NDArray[np.floating],
    min_sasa: float = 20.0,
    cluster_distance: float = 8.0
) -> Tuple[List[List[int]], NDArray[np.floating]]:
    """
    Cluster surface-exposed residues into patches.

    Surface patches are contiguous regions of solvent-exposed residues.
    These patches often correspond to potential epitope sites.

    Args:
        ca_coords: Nx3 array of CA coordinates
        sasa_values: N array of SASA values per residue (Angstrom^2)
        min_sasa: minimum SASA to consider surface-exposed (default 20.0)
        cluster_distance: maximum distance to join residues into patch (default 8.0)

    Returns:
        patches: list of lists, each containing residue indices (0-based) in a patch
        patch_scores: array of average SASA values per patch

    Raises:
        ValueError: If arrays have inconsistent lengths or invalid shapes

    Example:
        >>> coords = np.random.randn(50, 3) * 10
        >>> sasa = np.random.rand(50) * 100
        >>> patches, scores = detect_surface_patches(coords, sasa)
    """
    ca_coords = np.asarray(ca_coords, dtype=np.float64)
    sasa_values = np.asarray(sasa_values, dtype=np.float64).flatten()

    if ca_coords.ndim != 2 or ca_coords.shape[1] != 3:
        raise ValueError(f"ca_coords must be Nx3 array, got shape {ca_coords.shape}")
    if len(sasa_values) != ca_coords.shape[0]:
        raise ValueError(
            f"Length mismatch: {len(sasa_values)} SASA values, "
            f"{ca_coords.shape[0]} coordinates"
        )

    # Find surface-exposed residues
    surface_mask = sasa_values >= min_sasa
    surface_indices = np.where(surface_mask)[0]

    if len(surface_indices) == 0:
        return [], np.array([], dtype=np.float64)

    # Single-linkage clustering of surface residues
    # Build adjacency based on distance
    surface_coords = ca_coords[surface_indices]
    n_surface = len(surface_indices)

    # Distance matrix for surface residues only
    diff = surface_coords[:, np.newaxis, :] - surface_coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=2))

    # Adjacency: residues within cluster_distance
    adjacency = dist_matrix < cluster_distance

    # Connected components via iterative expansion
    visited = np.zeros(n_surface, dtype=bool)
    patches = []

    for start in range(n_surface):
        if visited[start]:
            continue

        # BFS to find connected component
        component = []
        queue = [start]
        visited[start] = True

        while queue:
            current = queue.pop(0)
            component.append(surface_indices[current])

            # Add unvisited neighbors
            neighbors = np.where(adjacency[current] & ~visited)[0]
            for neighbor in neighbors:
                visited[neighbor] = True
                queue.append(neighbor)

        patches.append(sorted(component))

    # Calculate average SASA per patch
    patch_scores = np.array([
        np.mean(sasa_values[patch]) for patch in patches
    ], dtype=np.float64)

    # Sort patches by score (descending)
    if len(patches) > 1:
        sorted_indices = np.argsort(patch_scores)[::-1]
        patches = [patches[i] for i in sorted_indices]
        patch_scores = patch_scores[sorted_indices]

    return patches, patch_scores


def calculate_structural_epitope_score(
    sasa: NDArray[np.floating],
    contacts: NDArray[np.floating],
    bfactors: NDArray[np.floating],
    depth: Optional[NDArray[np.floating]] = None,
    weights: Optional[Dict[str, float]] = None
) -> NDArray[np.floating]:
    """
    Combine structural features into single epitope likelihood score.

    This function normalizes each feature to [0, 1] and combines them
    using weighted sum. Features with positive weights contribute to
    higher epitope likelihood, negative weights reduce it.

    Default interpretation:
        - SASA (+): Higher surface accessibility -> more likely epitope
        - contacts (-): More contacts -> buried -> less likely epitope
        - bfactors (+): Higher flexibility -> more likely epitope
        - depth (-): Greater depth -> buried -> less likely epitope

    Args:
        sasa: N array of solvent accessible surface area values
        contacts: N array of contact counts
        bfactors: N array of B-factor values
        depth: Optional N array of residue depth values
        weights: Optional dict of feature weights. Missing features get
            default weights. Set weight to 0 to exclude a feature.

    Returns:
        scores: N array of composite epitope scores (0-1 range)

    Raises:
        ValueError: If arrays have inconsistent lengths

    Example:
        >>> n = 100
        >>> sasa = np.random.rand(n) * 100
        >>> contacts = np.random.randint(0, 20, n)
        >>> bfactors = np.random.rand(n) * 50
        >>> scores = calculate_structural_epitope_score(sasa, contacts, bfactors)
    """
    # Use default weights, updated with any provided weights
    effective_weights = DEFAULT_WEIGHTS.copy()
    if weights is not None:
        effective_weights.update(weights)

    # Convert inputs
    sasa = np.asarray(sasa, dtype=np.float64).flatten()
    contacts = np.asarray(contacts, dtype=np.float64).flatten()
    bfactors = np.asarray(bfactors, dtype=np.float64).flatten()

    n_residues = len(sasa)

    # Validate array lengths
    if len(contacts) != n_residues:
        raise ValueError(
            f"Length mismatch: sasa has {n_residues}, contacts has {len(contacts)}"
        )
    if len(bfactors) != n_residues:
        raise ValueError(
            f"Length mismatch: sasa has {n_residues}, bfactors has {len(bfactors)}"
        )
    if depth is not None:
        depth = np.asarray(depth, dtype=np.float64).flatten()
        if len(depth) != n_residues:
            raise ValueError(
                f"Length mismatch: sasa has {n_residues}, depth has {len(depth)}"
            )

    def normalize_feature(arr: np.ndarray) -> np.ndarray:
        """Normalize array to [0, 1] using min-max scaling."""
        arr_min = np.min(arr)
        arr_max = np.max(arr)
        if arr_max - arr_min < 1e-10:
            return np.full_like(arr, 0.5)
        return (arr - arr_min) / (arr_max - arr_min)

    # Normalize features
    sasa_norm = normalize_feature(sasa)
    contacts_norm = normalize_feature(contacts)
    bfactors_norm = normalize_feature(bfactors)

    # Compute weighted combination
    scores = np.zeros(n_residues, dtype=np.float64)
    total_abs_weight = 0.0

    # SASA contribution (positive weight = high SASA -> high score)
    w = effective_weights.get('sasa', 0)
    if abs(w) > 1e-10:
        scores += w * sasa_norm
        total_abs_weight += abs(w)

    # Contacts contribution (negative weight = high contacts -> low score)
    # Invert so that low contacts get high contribution
    w = effective_weights.get('contacts', 0)
    if abs(w) > 1e-10:
        if w < 0:
            scores += abs(w) * (1 - contacts_norm)
        else:
            scores += w * contacts_norm
        total_abs_weight += abs(w)

    # B-factors contribution
    w = effective_weights.get('bfactors', 0)
    if abs(w) > 1e-10:
        scores += w * bfactors_norm
        total_abs_weight += abs(w)

    # Depth contribution (if provided)
    if depth is not None:
        w = effective_weights.get('depth', 0)
        if abs(w) > 1e-10:
            depth_norm = normalize_feature(depth)
            if w < 0:
                scores += abs(w) * (1 - depth_norm)
            else:
                scores += w * depth_norm
            total_abs_weight += abs(w)

    # Normalize by total absolute weight to keep scores in [0, 1]
    if total_abs_weight > 1e-10:
        scores = scores / total_abs_weight

    # Clamp to [0, 1] for safety
    scores = np.clip(scores, 0.0, 1.0)

    return scores
