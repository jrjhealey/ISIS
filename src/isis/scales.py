"""
Amino acid property scales for B-cell epitope prediction.

These scales are published, peer-reviewed values from immunology literature.
Each scale assigns a numeric property value to each of the 20 standard amino acids.
"""

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

EMINI = {
    "A": 0.815, "C": 0.394, "D": 1.283, "E": 1.445, "F": 0.695,
    "G": 0.714, "H": 1.180, "I": 0.603, "K": 1.545, "L": 0.603,
    "M": 0.714, "N": 1.296, "P": 1.236, "Q": 1.348, "R": 1.475,
    "S": 1.115, "T": 1.184, "V": 0.606, "W": 0.808, "Y": 1.089,
}
"""
Emini surface accessibility scale.
Emini EA et al. J Virol. 1985;55(3):836-839. PMID: 2410642
Higher values indicate greater surface accessibility.
"""

PARKER = {
    "A": 2.1, "C": 1.4, "D": 10.0, "E": 7.8, "F": -9.2,
    "G": 5.7, "H": 2.1, "I": -8.0, "K": 5.7, "L": -9.2,
    "M": -4.2, "N": 7.0, "P": 2.1, "Q": 6.0, "R": 4.2,
    "S": 6.5, "T": 5.2, "V": -3.7, "W": -10.0, "Y": -1.9,
}
"""
Parker hydrophilicity scale.
Parker JM et al. Biochemistry. 1986;25(19):5425-5432. PMID: 3539163
Higher values indicate greater hydrophilicity.
"""

CHOU_FASMAN = {
    "A": 0.66, "C": 1.19, "D": 1.46, "E": 0.74, "F": 0.60,
    "G": 1.56, "H": 0.95, "I": 0.47, "K": 1.01, "L": 0.59,
    "M": 0.60, "N": 1.56, "P": 1.52, "Q": 0.98, "R": 0.95,
    "S": 1.43, "T": 0.96, "V": 0.50, "W": 0.96, "Y": 1.14,
}
"""
Chou-Fasman beta-turn propensity scale.
Chou PY, Fasman GD. Adv Enzymol. 1978;47:45-148. PMID: 364941
Higher values indicate greater propensity for beta-turns.
"""

KOLASKAR_TONGAONKAR = {
    "A": 1.064, "C": 1.412, "D": 0.866, "E": 0.851, "F": 1.091,
    "G": 0.874, "H": 1.105, "I": 1.152, "K": 0.930, "L": 1.250,
    "M": 0.826, "N": 0.776, "P": 1.064, "Q": 1.015, "R": 0.873,
    "S": 1.012, "T": 0.909, "V": 1.383, "W": 0.893, "Y": 1.161,
}
"""
Kolaskar-Tongaonkar antigenicity scale.
Kolaskar AS, Tongaonkar PC. FEBS Lett. 1990;276(1-2):172-174. PMID: 1702393
Values derived from experimentally determined antigenic determinants.
"""

KARPLUS_SCHULZ = {
    "A": 1.046, "C": 0.906, "D": 1.068, "E": 1.094, "F": 0.915,
    "G": 1.061, "H": 0.950, "I": 0.927, "K": 1.102, "L": 0.935,
    "M": 0.893, "N": 1.036, "P": 1.049, "Q": 1.037, "R": 1.008,
    "S": 1.046, "T": 0.997, "V": 0.912, "W": 0.904, "Y": 0.929,
}
"""
Karplus-Schulz flexibility scale.
Karplus PA, Schulz GE. Naturwissenschaften. 1985;72:212-213.
Based on B-factor mobility of alpha-carbons.
"""

SCALES = {
    "emini": EMINI,
    "parker": PARKER,
    "chou-fasman": CHOU_FASMAN,
    "kolaskar-tongaonkar": KOLASKAR_TONGAONKAR,
    "karplus-schulz": KARPLUS_SCHULZ,
}

METHOD_INFO = {
    "emini": {
        "name": "Emini Surface Accessibility",
        "default_window": 6,
        "default_threshold": 1.0,
        "description": "Predicts surface-exposed regions based on accessibility.",
    },
    "parker": {
        "name": "Parker Hydrophilicity",
        "default_window": 7,
        "default_threshold": None,  # Uses average
        "description": "Predicts hydrophilic regions likely to be antigenic.",
    },
    "chou-fasman": {
        "name": "Chou-Fasman Beta-Turn",
        "default_window": 7,
        "default_threshold": None,
        "description": "Predicts beta-turn regions associated with epitopes.",
    },
    "kolaskar-tongaonkar": {
        "name": "Kolaskar-Tongaonkar Antigenicity",
        "default_window": 7,
        "default_threshold": 1.0,
        "description": "Semi-empirical antigenicity prediction (~75% accuracy).",
    },
    "karplus-schulz": {
        "name": "Karplus-Schulz Flexibility",
        "default_window": 7,
        "default_threshold": 1.0,
        "description": "Predicts flexible regions likely to be immunogenic.",
    },
}


def get_scale(method: str) -> dict[str, float]:
    """Get amino acid scale by method name (case-insensitive)."""
    key = method.lower().replace("_", "-")
    if key not in SCALES:
        raise ValueError(f"Unknown method: {method}. Available: {list(SCALES.keys())}")
    return SCALES[key]


def get_method_info(method: str) -> dict:
    """Get method metadata by name."""
    key = method.lower().replace("_", "-")
    if key not in METHOD_INFO:
        raise ValueError(f"Unknown method: {method}")
    return METHOD_INFO[key]
