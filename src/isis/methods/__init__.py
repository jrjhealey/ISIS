"""
ISIS epitope prediction methods.

This module provides a unified interface for all epitope prediction methods
organized by category. The new method architecture supports B-cell linear,
B-cell conformational, T-cell, innate immunity, and structural methods.

Example:
    >>> from isis.methods import get_method, MethodCategory
    >>> method = get_method("emini")
    >>> result = method.predict("MKTAYIAKQRQISFVKSHFSRQLE")
    >>> print(result.epitopes)

    >>> from isis.methods import list_all_methods
    >>> methods = list_all_methods()
    >>> for key, info in methods.items():
    ...     print(f"{key}: {info['description']}")
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from .base import (
    EpitopeMethod,
    MethodCategory,
    MethodResult,
)
from .bcell_linear import (
    LinearScaleMethod,
    EminiMethod,
    ParkerMethod,
    ChouFasmanMethod,
    KolaskarTongaonkarMethod,
    KarplusSchulzMethod,
    BCELL_LINEAR_METHODS,
    get_method as get_bcell_linear_method,
    list_methods as list_bcell_linear_methods,
)
from .bcell_conformational import (
    ConformationalPredictor,
    ConformationalEpitope,
    ConformationalPrediction,
    DiscoTopePredictor,
    ElliProPredictor,
    SEPPAPredictor,
    calculate_protrusion_index,
    find_surface_patches,
    predict_discotope,
    predict_ellipro,
    predict_seppa,
    predict_all_conformational,
    DISCOTOPE_PROPENSITY,
)
from .innate import (
    InnatePredictor,
    GlycoSite,
    SignalPeptideResult,
    MotifMatch,
    predict_all_innate,
    NLinkedGlycosylationMethod,
    OLinkedGlycosylationMethod,
    SignalPeptideMethod,
    O_GLYCO_POSITIVE,
    O_GLYCO_NEGATIVE,
    HYDROPHOBICITY,
)

__all__ = [
    # Base classes and types
    "EpitopeMethod",
    "MethodCategory",
    "MethodResult",
    # B-cell linear methods
    "LinearScaleMethod",
    "EminiMethod",
    "ParkerMethod",
    "ChouFasmanMethod",
    "KolaskarTongaonkarMethod",
    "KarplusSchulzMethod",
    "BCELL_LINEAR_METHODS",
    # B-cell conformational methods
    "ConformationalPredictor",
    "ConformationalEpitope",
    "ConformationalPrediction",
    "DiscoTopePredictor",
    "ElliProPredictor",
    "SEPPAPredictor",
    "calculate_protrusion_index",
    "find_surface_patches",
    "predict_discotope",
    "predict_ellipro",
    "predict_seppa",
    "predict_all_conformational",
    "DISCOTOPE_PROPENSITY",
    "BCELL_CONFORMATIONAL_METHODS",
    # Lookup functions
    "get_method",
    "get_method_by_category",
    "list_all_methods",
    "list_methods_by_category",
    "available_methods",
]


# Registry of all methods across categories
_ALL_METHODS: Dict[str, type] = {
    **BCELL_LINEAR_METHODS,
    # Future: add other categories here
    # **TCELL_METHODS,
    # **INNATE_METHODS,
    # **STRUCTURAL_METHODS,
}

# Note: B-cell conformational methods (DiscoTope, ElliPro, SEPPA) are NOT in
# the registry because they require 3D structure data and use a different
# interface (predict_with_structure instead of predict). Use them directly:
#   from isis.methods import DiscoTopePredictor, ElliProPredictor, SEPPAPredictor
# or via convenience functions:
#   from isis.methods import predict_discotope, predict_ellipro, predict_seppa
BCELL_CONFORMATIONAL_METHODS = {
    "discotope": DiscoTopePredictor,
    "ellipro": ElliProPredictor,
    "seppa": SEPPAPredictor,
}

# Category to method registry mapping
_CATEGORY_REGISTRIES: Dict[MethodCategory, Dict[str, type]] = {
    MethodCategory.BCELL_LINEAR: BCELL_LINEAR_METHODS,
    # Future: add other categories here
    # MethodCategory.BCELL_CONFORMATIONAL: BCELL_CONFORMATIONAL_METHODS,
    # MethodCategory.TCELL: TCELL_METHODS,
    # MethodCategory.INNATE: INNATE_METHODS,
    # MethodCategory.STRUCTURAL: STRUCTURAL_METHODS,
}


def get_method(name: str) -> EpitopeMethod:
    """
    Get an instance of an epitope prediction method by name.

    This is the primary entry point for obtaining method instances.
    The method is looked up across all registered categories.

    Args:
        name: Method identifier (case-insensitive, underscores converted to hyphens).

    Returns:
        Instance of the requested EpitopeMethod.

    Raises:
        ValueError: If method name is not recognized.

    Example:
        >>> method = get_method("emini")
        >>> result = method.predict("ACDEFGHIKLMNPQRSTVWY")
    """
    key = name.lower().replace("_", "-")
    if key not in _ALL_METHODS:
        available = list(_ALL_METHODS.keys())
        raise ValueError(f"Unknown method: {name}. Available: {available}")
    return _ALL_METHODS[key]()


def get_method_by_category(
    name: str,
    category: MethodCategory
) -> EpitopeMethod:
    """
    Get a method instance, restricted to a specific category.

    Args:
        name: Method identifier.
        category: MethodCategory to search within.

    Returns:
        Instance of the requested EpitopeMethod.

    Raises:
        ValueError: If method not found in the specified category.
    """
    key = name.lower().replace("_", "-")
    registry = _CATEGORY_REGISTRIES.get(category, {})

    if key not in registry:
        available = list(registry.keys())
        raise ValueError(
            f"Unknown method '{name}' in category {category.value}. "
            f"Available: {available}"
        )
    return registry[key]()


def list_all_methods() -> Dict[str, Dict[str, Any]]:
    """
    List all available methods with their metadata.

    Returns:
        Dict mapping method keys to their info dicts containing:
        - name: Human-readable method name
        - key: Method identifier
        - category: Method category
        - description: Brief description
        - citation: Publication reference
        - default_window: Default sliding window size
        - default_threshold: Default score threshold
    """
    return {
        key: cls().get_info()
        for key, cls in _ALL_METHODS.items()
    }


def list_methods_by_category(
    category: Optional[MethodCategory] = None
) -> Dict[str, Dict[str, Any]]:
    """
    List methods filtered by category.

    Args:
        category: MethodCategory to filter by, or None for all methods.

    Returns:
        Dict mapping method keys to their info dicts.
    """
    if category is None:
        return list_all_methods()

    registry = _CATEGORY_REGISTRIES.get(category, {})
    return {
        key: cls().get_info()
        for key, cls in registry.items()
    }


def available_methods(
    category: Optional[MethodCategory] = None
) -> List[str]:
    """
    Get list of available method names.

    Args:
        category: Optional category to filter by.

    Returns:
        List of method key strings.
    """
    if category is None:
        return list(_ALL_METHODS.keys())
    registry = _CATEGORY_REGISTRIES.get(category, {})
    return list(registry.keys())
