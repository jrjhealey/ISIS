"""
Machine learning models for epitope prediction.

This module provides trained ML models for various epitope prediction tasks:
- B-cell linear epitope prediction (BepiPred v2.0 style Random Forest)
- MHC class II binding prediction
- Signal peptide prediction
"""
# Conditional imports for models with optional dependencies

# B-cell linear epitope model
try:
    from .bcell_model import (
        BcellLinearPredictor,
        get_predictor,
        predict_bcell_linear,
        train_and_save_model,
        # Feature extraction utilities
        get_residue_features,
        extract_window_features,
        one_hot_encode,
        # Training data generation
        generate_training_data,
        # Amino acid scales
        AMINO_ACIDS,
        KYTE_DOOLITTLE,
        EMINI,
        KARPLUS_SCHULZ,
        CHOU_FASMAN_TURN,
        KOLASKAR_TONGAONKAR,
        AA_CHARGE,
        AA_VOLUME,
    )
    _has_bcell = True
except ImportError:
    _has_bcell = False

# MHC class II model
try:
    from .mhc2_model import MHC2Predictor
    _has_mhc2 = True
except ImportError:
    _has_mhc2 = False

# Signal peptide model
try:
    from .signal_model import (
        SignalPeptidePredictor,
        SignalPeptideResult,
        predict_signal_peptide,
    )
    _has_signal = True
except ImportError:
    _has_signal = False

__all__ = []

if _has_bcell:
    __all__.extend([
        # Main predictor class
        "BcellLinearPredictor",
        # Convenience functions
        "get_predictor",
        "predict_bcell_linear",
        "train_and_save_model",
        # Feature extraction
        "get_residue_features",
        "extract_window_features",
        "one_hot_encode",
        # Training data
        "generate_training_data",
        # Scales
        "AMINO_ACIDS",
        "KYTE_DOOLITTLE",
        "EMINI",
        "KARPLUS_SCHULZ",
        "CHOU_FASMAN_TURN",
        "KOLASKAR_TONGAONKAR",
        "AA_CHARGE",
        "AA_VOLUME",
    ])

if _has_mhc2:
    __all__.append("MHC2Predictor")

if _has_signal:
    __all__.extend(["SignalPeptidePredictor", "SignalPeptideResult", "predict_signal_peptide"])
