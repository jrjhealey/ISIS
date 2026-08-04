"""
Signal Peptide Prediction Model

A lightweight scikit-learn based signal peptide predictor that identifies
signal peptides and their cleavage sites in protein sequences.

Signal peptide structure:
- N-region (1-5 residues): positively charged (K, R)
- H-region (7-15 residues): hydrophobic core (L, I, V, F, A)
- C-region (3-7 residues): small/neutral residues (A, G, S), ends with A-X-A pattern

Total length: typically 15-30 residues
"""

import base64
import pickle
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.preprocessing import StandardScaler
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False


# Amino acid property scales
HYDROPHOBICITY = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2,
    'X': 0.0, 'B': -3.5, 'Z': -3.5, 'U': 2.5, 'O': 0.0
}

# Charge at pH 7
CHARGE = {
    'K': 1.0, 'R': 1.0, 'H': 0.1,
    'D': -1.0, 'E': -1.0,
    'A': 0.0, 'N': 0.0, 'C': 0.0, 'Q': 0.0, 'G': 0.0,
    'I': 0.0, 'L': 0.0, 'M': 0.0, 'F': 0.0, 'P': 0.0,
    'S': 0.0, 'T': 0.0, 'W': 0.0, 'Y': 0.0, 'V': 0.0,
    'X': 0.0, 'B': -0.5, 'Z': -0.5, 'U': 0.0, 'O': 0.0
}

# Residue size (approximate Van der Waals volume)
SIZE = {
    'G': 60.1, 'A': 88.6, 'S': 89.0, 'C': 108.5, 'D': 111.1,
    'P': 112.7, 'N': 114.1, 'T': 116.1, 'E': 138.4, 'V': 140.0,
    'Q': 143.8, 'H': 153.2, 'M': 162.9, 'I': 166.7, 'L': 166.7,
    'K': 168.6, 'R': 173.4, 'F': 189.9, 'Y': 193.6, 'W': 227.8,
    'X': 140.0, 'B': 112.6, 'Z': 141.1, 'U': 108.5, 'O': 168.6
}

# Small residues for cleavage site prediction
SMALL_RESIDUES = {'A', 'G', 'S', 'C', 'T'}
HYDROPHOBIC_RESIDUES = {'L', 'I', 'V', 'F', 'A', 'M', 'W'}
POSITIVE_RESIDUES = {'K', 'R'}


@dataclass
class SignalPeptideResult:
    """Result of signal peptide prediction."""
    is_signal_peptide: bool
    probability: float
    cleavage_site: Optional[int]
    n_region: Optional[Tuple[int, int]]
    h_region: Optional[Tuple[int, int]]
    c_region: Optional[Tuple[int, int]]
    signal_sequence: Optional[str]
    mature_start: Optional[int]

    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'is_signal_peptide': self.is_signal_peptide,
            'probability': self.probability,
            'cleavage_site': self.cleavage_site,
            'n_region': self.n_region,
            'h_region': self.h_region,
            'c_region': self.c_region,
            'signal_sequence': self.signal_sequence,
            'mature_start': self.mature_start
        }


class SignalPeptidePredictor:
    """
    Lightweight signal peptide predictor using scikit-learn.

    Uses a two-stage approach:
    1. Binary classification: Does the sequence have a signal peptide?
    2. Cleavage site prediction: Where is the cleavage site?

    Features are based on signal peptide biology:
    - N-region: positively charged
    - H-region: hydrophobic core
    - C-region: small residues with A-X-A pattern

    Example:
        >>> predictor = SignalPeptidePredictor()
        >>> result = predictor.predict("MKFLIVAALVFMFYSYMAFS...")
        >>> print(result.is_signal_peptide)
        True
        >>> print(result.cleavage_site)
        21
    """

    # Feature dimension for signal peptide detection
    N_FEATURES = 52

    def __init__(self, pretrained: bool = True):
        """
        Initialize the signal peptide predictor.

        Args:
            pretrained: If True, load pre-trained model weights.
                       If False, create untrained models (for training).
        """
        if not HAS_SKLEARN:
            raise ImportError(
                "scikit-learn is required for SignalPeptidePredictor. "
                "Install with: pip install scikit-learn"
            )

        self.classifier = None
        self.cleavage_model = None
        self.scaler = None
        self._is_trained = False

        if pretrained:
            self._load_pretrained_weights()

    def extract_features(self, sequence: str, max_length: int = 70) -> np.ndarray:
        """
        Extract features for signal peptide detection.

        Features based on signal peptide biology:
        - N-region (1-5): positive charge density
        - H-region (6-20): hydrophobicity profile
        - C-region (15-30): small residue content, A-X-A patterns
        - Overall: length, composition, gradients

        Args:
            sequence: Protein sequence (first 50-70 residues recommended)
            max_length: Maximum sequence length to consider

        Returns:
            Feature vector of shape (N_FEATURES,)
        """
        seq = sequence[:max_length].upper()
        n = len(seq)

        if n < 15:
            # Too short for reliable signal peptide
            return np.zeros(self.N_FEATURES)

        features = []

        # === N-region features (residues 1-5) ===
        n_region = seq[:5] if n >= 5 else seq

        # Positive charge count in N-region
        n_positive = sum(1 for aa in n_region if aa in POSITIVE_RESIDUES)
        features.append(n_positive / max(len(n_region), 1))

        # Total charge in N-region
        n_charge = sum(CHARGE.get(aa, 0) for aa in n_region)
        features.append(n_charge)

        # Met at position 1 (common for signal peptides)
        features.append(1.0 if seq[0] == 'M' else 0.0)

        # === H-region features (residues 6-20) ===
        h_region = seq[5:20] if n >= 20 else seq[5:n]
        h_len = len(h_region)

        if h_len > 0:
            # Mean hydrophobicity
            h_hydro = np.mean([HYDROPHOBICITY.get(aa, 0) for aa in h_region])
            features.append(h_hydro)

            # Max hydrophobicity
            features.append(max(HYDROPHOBICITY.get(aa, 0) for aa in h_region))

            # Hydrophobic residue fraction
            h_hydrophobic_count = sum(1 for aa in h_region if aa in HYDROPHOBIC_RESIDUES)
            features.append(h_hydrophobic_count / h_len)

            # Hydrophobicity variance (should be low for consistent H-region)
            h_values = [HYDROPHOBICITY.get(aa, 0) for aa in h_region]
            features.append(np.std(h_values))

            # Length of potential H-region
            features.append(h_len / 15.0)
        else:
            features.extend([0.0, 0.0, 0.0, 0.0, 0.0])

        # === C-region features (residues 15-30) ===
        c_region = seq[15:30] if n >= 30 else seq[15:n] if n > 15 else ""
        c_len = len(c_region)

        if c_len > 0:
            # Small residue fraction
            c_small = sum(1 for aa in c_region if aa in SMALL_RESIDUES)
            features.append(c_small / c_len)

            # Mean hydrophobicity (should decrease from H-region)
            c_hydro = np.mean([HYDROPHOBICITY.get(aa, 0) for aa in c_region])
            features.append(c_hydro)

            # Hydrophobicity gradient (H to C should decrease)
            if h_len > 0:
                features.append(h_hydro - c_hydro)  # Positive for signal peptides
            else:
                features.append(0.0)

            # A-X-A pattern detection in C-region
            axa_score = self._score_axa_patterns(c_region)
            features.append(axa_score)
        else:
            features.extend([0.0, 0.0, 0.0, 0.0])

        # === Position-specific hydrophobicity profile (10 bins) ===
        profile_length = min(n, 50)
        hydro_profile = [HYDROPHOBICITY.get(aa, 0) for aa in seq[:profile_length]]

        # Bin into 10 windows
        bin_size = max(profile_length // 10, 1)
        for i in range(10):
            start = i * bin_size
            end = min(start + bin_size, profile_length)
            if start < profile_length:
                features.append(np.mean(hydro_profile[start:end]))
            else:
                features.append(0.0)

        # === Position-specific charge profile (5 bins) ===
        charge_profile = [CHARGE.get(aa, 0) for aa in seq[:profile_length]]
        bin_size = max(profile_length // 5, 1)
        for i in range(5):
            start = i * bin_size
            end = min(start + bin_size, profile_length)
            if start < profile_length:
                features.append(np.sum(charge_profile[start:end]))
            else:
                features.append(0.0)

        # === Composition features ===
        # Overall amino acid composition for key residues
        for aa in ['L', 'A', 'V', 'I', 'F', 'K', 'R', 'G', 'S']:
            count = seq[:40].count(aa)
            features.append(count / min(n, 40))

        # === Structural features ===
        # Charged residues in first 10
        first_10 = seq[:10]
        features.append(sum(1 for aa in first_10 if aa in POSITIVE_RESIDUES) / 10)
        features.append(sum(1 for aa in first_10 if aa in {'D', 'E'}) / 10)

        # Hydrophobic stretch detection
        max_hydro_stretch = self._find_hydrophobic_stretch(seq[:40])
        features.append(max_hydro_stretch[0] / 40)  # Start position
        features.append(max_hydro_stretch[1] / 20)  # Length
        features.append(max_hydro_stretch[2])  # Mean hydrophobicity

        # === Cleavage site signatures ===
        # Score potential cleavage sites in typical range (15-35)
        cleavage_scores = []
        for pos in range(15, min(36, n)):
            score = self._score_cleavage_site(seq, pos)
            cleavage_scores.append(score)

        if cleavage_scores:
            features.append(max(cleavage_scores))  # Best cleavage score
            features.append(np.mean(cleavage_scores))  # Mean cleavage score
        else:
            features.append(0.0)
            features.append(0.0)

        # Sequence length
        features.append(min(n, 70) / 70)

        # Pad or truncate to fixed size
        features = np.array(features[:self.N_FEATURES])
        if len(features) < self.N_FEATURES:
            features = np.pad(features, (0, self.N_FEATURES - len(features)))

        return features.astype(np.float32)

    def _score_axa_patterns(self, region: str) -> float:
        """
        Score A-X-A type cleavage patterns.

        The -3, -1 rule: small residues at positions -3 and -1 relative to cleavage.
        Common patterns: AXA, GXA, SXA, AXS, AXG
        """
        if len(region) < 3:
            return 0.0

        max_score = 0.0
        for i in range(len(region) - 2):
            pos_minus3 = region[i]
            pos_minus1 = region[i + 2]

            score = 0.0

            # Position -3 scoring
            if pos_minus3 == 'A':
                score += 1.0
            elif pos_minus3 in SMALL_RESIDUES:
                score += 0.7
            elif pos_minus3 in {'V', 'L', 'I'}:
                score += 0.3

            # Position -1 scoring (most important)
            if pos_minus1 == 'A':
                score += 1.2
            elif pos_minus1 in SMALL_RESIDUES:
                score += 0.8

            # Bonus for classic A-X-A
            if pos_minus3 == 'A' and pos_minus1 == 'A':
                score += 0.5

            max_score = max(max_score, score)

        return max_score / 2.7  # Normalize to [0, 1]

    def _find_hydrophobic_stretch(self, seq: str) -> Tuple[int, int, float]:
        """
        Find the longest hydrophobic stretch (H-region candidate).

        Returns: (start_position, length, mean_hydrophobicity)
        """
        if not seq:
            return (0, 0, 0.0)

        best_start = 0
        best_length = 0
        best_hydro = 0.0

        # Sliding window approach
        min_hydro_threshold = 1.0  # Minimum mean hydrophobicity for H-region

        for window_size in range(7, 16):
            for start in range(len(seq) - window_size + 1):
                window = seq[start:start + window_size]
                hydro_values = [HYDROPHOBICITY.get(aa, 0) for aa in window]
                mean_hydro = np.mean(hydro_values)

                if mean_hydro >= min_hydro_threshold:
                    # Score by length * hydrophobicity
                    score = window_size * mean_hydro
                    if score > best_length * best_hydro:
                        best_start = start
                        best_length = window_size
                        best_hydro = mean_hydro

        return (best_start, best_length, best_hydro)

    def _score_cleavage_site(self, seq: str, position: int) -> float:
        """
        Score a potential cleavage site at the given position.

        Based on the -3, -1 rule: small residues at positions -3 and -1.
        Also considers surrounding context.
        """
        if position < 3 or position >= len(seq):
            return 0.0

        score = 0.0

        # Position -1 (most important)
        aa_minus1 = seq[position - 1]
        if aa_minus1 == 'A':
            score += 2.0
        elif aa_minus1 in {'G', 'S', 'C'}:
            score += 1.5
        elif aa_minus1 in {'T', 'V'}:
            score += 0.5

        # Position -3
        aa_minus3 = seq[position - 3]
        if aa_minus3 == 'A':
            score += 1.5
        elif aa_minus3 in {'G', 'S', 'C', 'V', 'L', 'I'}:
            score += 1.0
        elif aa_minus3 in {'T', 'F'}:
            score += 0.3

        # Position -2 (often helix-breaking)
        if position >= 2:
            aa_minus2 = seq[position - 2]
            if aa_minus2 in {'P', 'G'}:
                score -= 0.5  # Helix breakers not favored at -2

        # Bonus for proper context (hydrophobic before, less hydrophobic after)
        if position >= 10 and position < len(seq) - 5:
            before = seq[position-10:position-3]
            after = seq[position:position+5]

            hydro_before = np.mean([HYDROPHOBICITY.get(aa, 0) for aa in before])
            hydro_after = np.mean([HYDROPHOBICITY.get(aa, 0) for aa in after])

            if hydro_before > 1.5 and hydro_after < hydro_before:
                score += 1.0

        return score

    def _extract_cleavage_features(self, seq: str, position: int) -> np.ndarray:
        """
        Extract features for cleavage site prediction at a specific position.
        """
        features = []

        # Position in sequence
        features.append(position / 40.0)

        # Residues around the cleavage site
        for offset in [-5, -4, -3, -2, -1, 0, 1, 2]:
            pos = position + offset
            if 0 <= pos < len(seq):
                aa = seq[pos]
                features.append(HYDROPHOBICITY.get(aa, 0) / 4.5)
                features.append(1.0 if aa in SMALL_RESIDUES else 0.0)
                features.append(CHARGE.get(aa, 0))
            else:
                features.extend([0.0, 0.0, 0.0])

        # Cleavage site score
        features.append(self._score_cleavage_site(seq, position) / 5.0)

        # A-X-A pattern
        if position >= 3:
            features.append(self._score_axa_patterns(seq[position-6:position]) if position >= 6 else 0.0)
        else:
            features.append(0.0)

        return np.array(features, dtype=np.float32)

    def predict(self, sequence: str) -> SignalPeptideResult:
        """
        Predict signal peptide presence and cleavage site.

        Args:
            sequence: Protein sequence (full or first 50-70 residues)

        Returns:
            SignalPeptideResult with prediction details
        """
        if not self._is_trained:
            raise RuntimeError("Model is not trained. Call train() or use pretrained=True.")

        # Clean sequence
        seq = sequence.strip().upper()[:70]

        if len(seq) < 15:
            return SignalPeptideResult(
                is_signal_peptide=False,
                probability=0.0,
                cleavage_site=None,
                n_region=None,
                h_region=None,
                c_region=None,
                signal_sequence=None,
                mature_start=None
            )

        # Stage 1: Classify if sequence has signal peptide
        features = self.extract_features(seq)
        features_scaled = self.scaler.transform(features.reshape(1, -1))

        proba = self.classifier.predict_proba(features_scaled)[0]
        is_signal = proba[1] > 0.5
        probability = float(proba[1])

        if not is_signal:
            return SignalPeptideResult(
                is_signal_peptide=False,
                probability=probability,
                cleavage_site=None,
                n_region=None,
                h_region=None,
                c_region=None,
                signal_sequence=None,
                mature_start=None
            )

        # Stage 2: Predict cleavage site
        cleavage_site = self._predict_cleavage_site(seq)

        # Determine regions
        n_region, h_region, c_region = self._determine_regions(seq, cleavage_site)

        return SignalPeptideResult(
            is_signal_peptide=True,
            probability=probability,
            cleavage_site=cleavage_site,
            n_region=n_region,
            h_region=h_region,
            c_region=c_region,
            signal_sequence=seq[:cleavage_site] if cleavage_site else None,
            mature_start=cleavage_site
        )

    def _predict_cleavage_site(self, seq: str) -> int:
        """
        Predict the most likely cleavage site position.

        Uses a combination of:
        1. Position weight matrix scoring
        2. Machine learning model
        3. Hydrophobicity gradient
        """
        n = len(seq)
        best_position = 20  # Default
        best_score = -float('inf')

        # Search in typical signal peptide length range
        for pos in range(15, min(36, n)):
            # Base score from cleavage site rules
            rule_score = self._score_cleavage_site(seq, pos)

            # ML score if model available
            if self.cleavage_model is not None:
                cleavage_features = self._extract_cleavage_features(seq, pos)
                ml_score = self.cleavage_model.predict_proba(
                    cleavage_features.reshape(1, -1)
                )[0, 1]
            else:
                ml_score = 0.5

            # Combined score
            total_score = 0.6 * rule_score + 0.4 * ml_score * 5

            if total_score > best_score:
                best_score = total_score
                best_position = pos

        return best_position

    def _determine_regions(
        self, seq: str, cleavage_site: int
    ) -> Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]:
        """
        Determine N, H, and C regions based on sequence and cleavage site.
        """
        # N-region: typically first 1-5 residues
        n_end = min(5, cleavage_site - 10) if cleavage_site > 15 else 3
        n_region = (1, n_end)

        # H-region: hydrophobic core
        h_start, h_length, _ = self._find_hydrophobic_stretch(seq[:cleavage_site])
        if h_length > 0:
            h_region = (h_start + 1, h_start + h_length)
        else:
            h_region = (n_end + 1, cleavage_site - 5)

        # C-region: from end of H-region to cleavage site
        c_start = h_region[1] + 1 if h_region else n_end + 10
        c_region = (c_start, cleavage_site)

        return n_region, h_region, c_region

    def train(
        self,
        positive_sequences: List[str],
        negative_sequences: List[str],
        cleavage_sites: Optional[List[int]] = None
    ) -> Dict[str, float]:
        """
        Train the signal peptide predictor.

        Args:
            positive_sequences: Sequences with signal peptides
            negative_sequences: Sequences without signal peptides
            cleavage_sites: Known cleavage sites for positive sequences

        Returns:
            Dictionary with training metrics
        """
        # Prepare training data
        X_pos = np.array([self.extract_features(seq) for seq in positive_sequences])
        X_neg = np.array([self.extract_features(seq) for seq in negative_sequences])

        X = np.vstack([X_pos, X_neg])
        y = np.array([1] * len(positive_sequences) + [0] * len(negative_sequences))

        # Scale features
        self.scaler = StandardScaler()
        X_scaled = self.scaler.fit_transform(X)

        # Train classifier
        self.classifier = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            min_samples_split=5,
            min_samples_leaf=2,
            random_state=42,
            class_weight='balanced'
        )
        self.classifier.fit(X_scaled, y)

        # Train cleavage site model if sites provided
        if cleavage_sites:
            self._train_cleavage_model(positive_sequences, cleavage_sites)

        self._is_trained = True

        # Calculate training accuracy
        train_pred = self.classifier.predict(X_scaled)
        accuracy = np.mean(train_pred == y)

        return {'accuracy': accuracy, 'n_samples': len(y)}

    def _train_cleavage_model(
        self, sequences: List[str], cleavage_sites: List[int]
    ) -> None:
        """Train the cleavage site prediction model."""
        X_cleavage = []
        y_cleavage = []

        for seq, true_site in zip(sequences, cleavage_sites):
            # Add positive example at true cleavage site
            X_cleavage.append(self._extract_cleavage_features(seq, true_site))
            y_cleavage.append(1)

            # Add negative examples at other positions
            for pos in range(15, min(36, len(seq))):
                if abs(pos - true_site) > 2:  # Not too close to true site
                    X_cleavage.append(self._extract_cleavage_features(seq, pos))
                    y_cleavage.append(0)

        X_cleavage = np.array(X_cleavage)
        y_cleavage = np.array(y_cleavage)

        self.cleavage_model = GradientBoostingClassifier(
            n_estimators=50,
            max_depth=5,
            random_state=42
        )
        self.cleavage_model.fit(X_cleavage, y_cleavage)

    def _load_pretrained_weights(self) -> None:
        """Load pre-trained model weights."""
        # Generate synthetic training data and train
        positive_seqs, negative_seqs, cleavage_sites = self._generate_training_data()
        self.train(positive_seqs, negative_seqs, cleavage_sites)

    def _generate_training_data(
        self, n_positive: int = 500, n_negative: int = 500
    ) -> Tuple[List[str], List[str], List[int]]:
        """
        Generate synthetic training data following signal peptide rules.

        Based on biological constraints:
        - N-region: 1-5 residues, positively charged
        - H-region: 7-15 residues, highly hydrophobic
        - C-region: 3-7 residues, small/neutral with A-X-A pattern
        """
        random.seed(42)  # Reproducibility

        positive_seqs = []
        cleavage_sites = []

        for _ in range(n_positive):
            seq, site = self._generate_signal_peptide()
            positive_seqs.append(seq)
            cleavage_sites.append(site)

        negative_seqs = [self._generate_non_signal() for _ in range(n_negative)]

        return positive_seqs, negative_seqs, cleavage_sites

    def _generate_signal_peptide(self) -> Tuple[str, int]:
        """Generate a synthetic signal peptide with mature protein start."""
        # N-region: 1-5 residues, starts with M, followed by positive residues
        n_length = random.randint(1, 5)
        n_region = ['M'] + random.choices(['K', 'R', 'K', 'R', 'S', 'A'], k=n_length - 1)

        # H-region: 7-15 hydrophobic residues
        h_length = random.randint(8, 15)
        h_residues = ['L', 'I', 'V', 'A', 'F', 'L', 'L', 'V', 'A', 'I']  # Weighted
        h_region = random.choices(h_residues, k=h_length)

        # C-region: 3-7 residues ending with A-X-A pattern
        c_length = random.randint(3, 7)
        c_residues = ['A', 'G', 'S', 'T', 'A', 'A', 'G', 'S']  # Weighted toward A
        c_region = random.choices(c_residues, k=c_length - 3)

        # Ensure A-X-A pattern at end
        x_residue = random.choice(['L', 'A', 'V', 'G', 'S', 'T', 'F', 'Q', 'N'])
        c_region.extend(['A', x_residue, 'A'])

        # Mature protein start (first ~30 residues)
        all_residues = list('ACDEFGHIKLMNPQRSTVWY')
        mature_length = random.randint(25, 40)
        mature_region = random.choices(all_residues, k=mature_length)

        signal_seq = ''.join(n_region + h_region + c_region)
        full_seq = signal_seq + ''.join(mature_region)
        cleavage_site = len(signal_seq)

        return full_seq, cleavage_site

    def _generate_non_signal(self) -> str:
        """Generate a non-signal peptide sequence."""
        all_residues = list('ACDEFGHIKLMNPQRSTVWY')

        # Various non-signal patterns
        pattern = random.choice(['cytoplasmic', 'membrane', 'random'])

        if pattern == 'cytoplasmic':
            # Start with M, then mixed composition
            seq = ['M'] + random.choices(all_residues, k=random.randint(50, 70))
            # Add some charged residues in N-terminus
            for i in range(1, min(10, len(seq))):
                if random.random() < 0.3:
                    seq[i] = random.choice(['D', 'E', 'K', 'R'])

        elif pattern == 'membrane':
            # Transmembrane-like: hydrophobic stretch but not at right position
            seq = ['M']
            # Some random residues
            seq.extend(random.choices(all_residues, k=random.randint(15, 25)))
            # Hydrophobic stretch (too late for signal peptide)
            seq.extend(random.choices(['L', 'I', 'V', 'A', 'F'], k=random.randint(15, 22)))
            # More random
            seq.extend(random.choices(all_residues, k=random.randint(10, 20)))

        else:  # random
            length = random.randint(50, 70)
            seq = ['M'] + random.choices(all_residues, k=length - 1)

        return ''.join(seq)

    def save(self, filepath: str) -> None:
        """Save model weights to file."""
        try:
            import joblib
            model_data = {
                'classifier': self.classifier,
                'cleavage_model': self.cleavage_model,
                'scaler': self.scaler,
                'is_trained': self._is_trained
            }
            joblib.dump(model_data, filepath)
        except ImportError:
            # Fallback to pickle
            model_data = {
                'classifier': self.classifier,
                'cleavage_model': self.cleavage_model,
                'scaler': self.scaler,
                'is_trained': self._is_trained
            }
            with open(filepath, 'wb') as f:
                pickle.dump(model_data, f)

    def load(self, filepath: str) -> None:
        """Load model weights from file."""
        try:
            import joblib
            model_data = joblib.load(filepath)
        except ImportError:
            with open(filepath, 'rb') as f:
                model_data = pickle.load(f)

        self.classifier = model_data['classifier']
        self.cleavage_model = model_data['cleavage_model']
        self.scaler = model_data['scaler']
        self._is_trained = model_data['is_trained']


# Convenience function for quick predictions
def predict_signal_peptide(sequence: str) -> SignalPeptideResult:
    """
    Quick signal peptide prediction.

    Args:
        sequence: Protein sequence

    Returns:
        SignalPeptideResult with prediction

    Example:
        >>> result = predict_signal_peptide("MKTAYIAKQRQISFVK...")
        >>> if result.is_signal_peptide:
        ...     print(f"Cleavage at position {result.cleavage_site}")
    """
    predictor = SignalPeptidePredictor(pretrained=True)
    return predictor.predict(sequence)


# Position weight matrix for cleavage site (based on SignalP data patterns)
CLEAVAGE_PWM = {
    -3: {'A': 0.45, 'V': 0.15, 'S': 0.10, 'G': 0.08, 'L': 0.07, 'I': 0.05, 'T': 0.05, 'C': 0.03, 'F': 0.02},
    -1: {'A': 0.55, 'G': 0.15, 'S': 0.12, 'C': 0.06, 'T': 0.05, 'V': 0.04, 'L': 0.02, 'Q': 0.01},
    1:  {'A': 0.20, 'Q': 0.12, 'S': 0.10, 'D': 0.09, 'E': 0.08, 'G': 0.07, 'T': 0.07, 'N': 0.06, 'K': 0.05},
}


def score_cleavage_position(sequence: str, position: int) -> float:
    """
    Score a potential cleavage site using the position weight matrix.

    Args:
        sequence: Protein sequence
        position: Potential cleavage site (1-indexed, cleavage after this position)

    Returns:
        Log-odds score for cleavage at this position
    """
    score = 0.0

    for offset, weights in CLEAVAGE_PWM.items():
        pos = position + offset
        if 0 <= pos < len(sequence):
            aa = sequence[pos]
            prob = weights.get(aa, 0.01)
            background = 0.05  # Background frequency
            score += np.log(prob / background)

    return score
