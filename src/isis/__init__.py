"""
ISIS - In Silico Immunogenicity Studies

A library for predicting B-cell epitopes using amino acid property scales.

Example:
    >>> from isis import predict, predict_all
    >>> result = predict("MKTAYIAKQRQISFVKSHFSRQLE", method="emini")
    >>> print(result.epitopes)
"""

from .predictor import (
    predict,
    predict_all,
    available_methods,
    Prediction,
    Epitope,
)
from .scales import (
    SCALES,
    METHOD_INFO,
    get_scale,
    get_method_info,
)

__version__ = "2.0.0"
__all__ = [
    "predict",
    "predict_all",
    "available_methods",
    "Prediction",
    "Epitope",
    "SCALES",
    "METHOD_INFO",
    "get_scale",
    "get_method_info",
]
