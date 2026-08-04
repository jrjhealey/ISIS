"""
Innate immunity prediction methods.

Predicts glycosylation sites, signal peptides, TLR ligand motifs,
and other features relevant to innate immune recognition.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .base import EpitopeMethod, MethodCategory, MethodResult


# Amino acid property sets
HYDROPHOBIC = set("AILMFVWPY")
BASIC = set("KRH")
POLAR = set("STNQ")
ACIDIC = set("DE")
SMALL = set("AGST")

# Hydrophobicity scale (Kyte-Doolittle)
HYDROPHOBICITY = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

# O-Glycosylation propensity scales
O_GLYCO_POSITIVE = {"P": 0.3, "T": 0.2, "S": 0.2, "A": 0.1, "V": 0.1}
O_GLYCO_NEGATIVE = {"E": -0.2, "D": -0.2, "K": -0.1, "R": -0.1}

# Signal peptide cleavage site preferences (small, uncharged at -1 and -3)
SP_CLEAVAGE_RESIDUES = set("AGSTCVN")


@dataclass
class GlycoSite:
    """A predicted glycosylation site."""
    position: int  # 1-indexed position of N/S/T
    score: float
    motif: str
    site_type: str  # "canonical" or "non-canonical" for N-glyco; "high"/"medium"/"low" for O-glyco

    def to_dict(self) -> Dict[str, Any]:
        return {
            "position": self.position,
            "score": self.score,
            "motif": self.motif,
            "site_type": self.site_type,
        }


@dataclass
class SignalPeptideResult:
    """Result of signal peptide prediction."""
    is_signal_peptide: bool
    cleavage_position: Optional[int]  # 1-indexed position after cleavage
    confidence: float
    n_region_score: float
    h_region_score: float
    c_region_score: float
    details: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "is_signal_peptide": self.is_signal_peptide,
            "cleavage_position": self.cleavage_position,
            "confidence": self.confidence,
            "n_region_score": self.n_region_score,
            "h_region_score": self.h_region_score,
            "c_region_score": self.c_region_score,
            "details": self.details,
        }


@dataclass
class MotifMatch:
    """A matched sequence motif."""
    position: int  # 1-indexed start position
    end: int  # 1-indexed end position
    sequence: str
    motif_name: str
    score: float
    annotation: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "position": self.position,
            "end": self.end,
            "sequence": self.sequence,
            "motif_name": self.motif_name,
            "score": self.score,
            "annotation": self.annotation,
        }


class InnatePredictor:
    """
    Predictor for innate immunity features.

    Predicts N-glycosylation sites, O-glycosylation sites, signal peptides,
    TLR ligand motifs, and immunostimulatory motifs.
    """

    def __init__(self):
        """Initialize the innate immunity predictor."""
        # Flagellin conserved regions (TLR5)
        # N-terminal D0 domain conserved region
        self.flagellin_n_pattern = re.compile(r"[AILMV]{2,4}[ST][AILMV]{1,3}N")
        # C-terminal D0 domain conserved region
        self.flagellin_c_pattern = re.compile(r"[QN][AILMV]{2,3}[ST][AILMV]{2}")

        # LPS-binding motifs (basic/hydrophobic)
        self.lps_binding_pattern = re.compile(r"[KR]{2,}[AILMFVWY]{2,}|[AILMFVWY]{2,}[KR]{2,}")

        # CpG motifs (for nucleic acid sequences)
        self.cpg_pattern = re.compile(r"[ATGC]*CG[ATGC]*", re.IGNORECASE)

    def predict_n_glyco(self, sequence: str) -> Dict[str, Any]:
        """
        Predict N-glycosylation sites.

        N-glycosylation occurs at Asn (N) in the sequon N-X-S/T where X != P.
        Also detects non-canonical weak sites with relaxed requirements.

        Args:
            sequence: Amino acid sequence (one-letter codes).

        Returns:
            Dictionary with keys:
                - positions: List of 1-indexed N positions with potential glycosylation.
                - scores: List of scores (1.0 for canonical, <1.0 for non-canonical).
                - motifs: List of matched sequon strings.
                - annotations: List of site type descriptions.
                - sites: List of GlycoSite objects.
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        sites: List[GlycoSite] = []

        for i, aa in enumerate(sequence):
            if aa != "N":
                continue

            # Need at least 2 more residues for the sequon
            if i + 2 >= len(sequence):
                continue

            x_residue = sequence[i + 1]
            third_residue = sequence[i + 2]

            # Canonical sequon: N-X-S/T where X != P
            if x_residue != "P" and third_residue in "ST":
                motif = sequence[i:i + 3]
                score = 1.0

                # Adjust score based on favorable/unfavorable factors
                # Threonine at +2 slightly more favorable than Serine
                if third_residue == "T":
                    score *= 1.05

                # Acidic or proline near the site reduces efficiency
                if i > 0 and sequence[i - 1] in "DEP":
                    score *= 0.8
                if i + 3 < len(sequence) and sequence[i + 3] in "DEP":
                    score *= 0.9

                # Cap at 1.0
                score = min(score, 1.0)

                sites.append(GlycoSite(
                    position=i + 1,  # 1-indexed
                    score=score,
                    motif=motif,
                    site_type="canonical"
                ))

            # Non-canonical sites: N-X-C (cysteine variant, rare)
            elif x_residue != "P" and third_residue == "C":
                motif = sequence[i:i + 3]
                sites.append(GlycoSite(
                    position=i + 1,
                    score=0.3,  # Weak site
                    motif=motif,
                    site_type="non-canonical-NXC"
                ))

            # N-P-S/T is blocked but may have weak activity in some contexts
            elif x_residue == "P" and third_residue in "ST":
                motif = sequence[i:i + 3]
                sites.append(GlycoSite(
                    position=i + 1,
                    score=0.1,  # Very weak, usually blocked
                    motif=motif,
                    site_type="non-canonical-blocked"
                ))

        return {
            "positions": [s.position for s in sites],
            "scores": [s.score for s in sites],
            "motifs": [s.motif for s in sites],
            "annotations": [s.site_type for s in sites],
            "sites": [s.to_dict() for s in sites],
        }

    def predict_o_glyco(self, sequence: str) -> Dict[str, Any]:
        """
        Predict O-glycosylation sites.

        O-glycosylation primarily targets Serine (S) and Threonine (T) residues.
        Prediction uses propensity scoring based on surrounding residues.

        Favorable factors:
            - Proline at -1 or +3 position
            - High local S/T density
            - Presence of other favorable residues (A, V)

        Args:
            sequence: Amino acid sequence (one-letter codes).

        Returns:
            Dictionary with keys:
                - positions: List of 1-indexed S/T positions with predicted O-glycosylation.
                - scores: List of propensity scores.
                - motifs: List of local sequence contexts.
                - annotations: List of confidence levels.
                - sites: List of GlycoSite objects.
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        sites: List[GlycoSite] = []

        for i, aa in enumerate(sequence):
            if aa not in "ST":
                continue

            score = 0.5  # Base score for S/T

            # Get flanking region (5 residues on each side)
            start = max(0, i - 5)
            end = min(len(sequence), i + 6)
            flanking = sequence[start:end]

            # Count local S/T density (excluding current position)
            local_st_count = sum(1 for c in flanking if c in "ST") - 1
            st_density_bonus = local_st_count * 0.05
            score += st_density_bonus

            # Proline at -1 position strongly increases probability
            if i > 0 and sequence[i - 1] == "P":
                score += 0.25

            # Proline at +3 position increases probability
            if i + 3 < len(sequence) and sequence[i + 3] == "P":
                score += 0.2

            # Sum propensity contributions from flanking residues
            for j, res in enumerate(flanking):
                if j == i - start:  # Skip the central S/T
                    continue
                if res in O_GLYCO_POSITIVE:
                    score += O_GLYCO_POSITIVE[res] * 0.5
                if res in O_GLYCO_NEGATIVE:
                    score += O_GLYCO_NEGATIVE[res] * 0.5

            # Threonine is slightly more favorable than Serine
            if aa == "T":
                score *= 1.1

            # Clamp score to [0, 1]
            score = max(0.0, min(1.0, score))

            # Determine confidence level
            if score >= 0.7:
                site_type = "high"
            elif score >= 0.5:
                site_type = "medium"
            else:
                site_type = "low"

            # Get motif context (5 residues centered on S/T)
            motif_start = max(0, i - 2)
            motif_end = min(len(sequence), i + 3)
            motif = sequence[motif_start:motif_end]

            sites.append(GlycoSite(
                position=i + 1,  # 1-indexed
                score=round(score, 3),
                motif=motif,
                site_type=site_type
            ))

        # Filter to return only sites with reasonable scores
        significant_sites = [s for s in sites if s.score >= 0.4]

        return {
            "positions": [s.position for s in significant_sites],
            "scores": [s.score for s in significant_sites],
            "motifs": [s.motif for s in significant_sites],
            "annotations": [s.site_type for s in significant_sites],
            "sites": [s.to_dict() for s in significant_sites],
            "all_sites": [s.to_dict() for s in sites],  # Include all for reference
        }

    def predict_signal_peptide(self, sequence: str) -> Dict[str, Any]:
        """
        Predict signal peptide presence and cleavage site.

        Signal peptides typically have three regions:
            - N-region (1-5 aa): Positively charged (K, R)
            - H-region (7-15 aa): Hydrophobic core
            - C-region (3-7 aa): Polar, with cleavage site following A-X-A rule

        Args:
            sequence: Amino acid sequence (one-letter codes).

        Returns:
            Dictionary with keys:
                - is_signal_peptide: Boolean indicating likely signal peptide.
                - cleavage_position: 1-indexed position after the cleavage site.
                - confidence: Confidence score (0-1).
                - n_region_score: Score for N-region (positive charges).
                - h_region_score: Score for H-region (hydrophobicity).
                - c_region_score: Score for C-region (cleavage site).
                - details: Additional analysis details.
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")

        if len(sequence) < 15:
            return SignalPeptideResult(
                is_signal_peptide=False,
                cleavage_position=None,
                confidence=0.0,
                n_region_score=0.0,
                h_region_score=0.0,
                c_region_score=0.0,
                details={"error": "Sequence too short for signal peptide analysis"}
            ).to_dict()

        # Analyze N-region (first 5 residues)
        n_region = sequence[:5]
        n_basic_count = sum(1 for aa in n_region if aa in BASIC)
        n_region_score = min(1.0, n_basic_count * 0.4)

        # Bonus for positive charge at position 2
        if len(sequence) > 1 and sequence[1] in BASIC:
            n_region_score = min(1.0, n_region_score + 0.2)

        # Analyze H-region (residues 6-20)
        h_region = sequence[5:20] if len(sequence) >= 20 else sequence[5:]
        h_hydro_scores = [HYDROPHOBICITY.get(aa, 0) for aa in h_region]
        avg_hydrophobicity = np.mean(h_hydro_scores) if h_hydro_scores else 0

        # Normalize hydrophobicity to 0-1 scale
        # Kyte-Doolittle ranges from -4.5 to 4.5, typical SP H-region is 1.5-3.0
        h_region_score = max(0, min(1.0, (avg_hydrophobicity + 1) / 4))

        # Find potential cleavage sites (A-X-A motif, positions 15-30)
        best_cleavage_pos = None
        best_c_score = 0.0

        for i in range(14, min(30, len(sequence) - 1)):
            # Check -3 and -1 positions relative to cleavage (A-X-A rule)
            if i < 2:
                continue

            pos_minus_3 = sequence[i - 2]  # -3 position
            pos_minus_1 = sequence[i]      # -1 position

            c_score = 0.0

            # Small/uncharged at -1 (Alanine most common)
            if pos_minus_1 == "A":
                c_score += 0.4
            elif pos_minus_1 in SP_CLEAVAGE_RESIDUES:
                c_score += 0.25

            # Small/uncharged at -3
            if pos_minus_3 == "A":
                c_score += 0.3
            elif pos_minus_3 in SP_CLEAVAGE_RESIDUES:
                c_score += 0.2

            # Avoid charged residues at -1
            if pos_minus_1 in BASIC or pos_minus_1 in ACIDIC:
                c_score -= 0.3

            # Check -2 position (often helix-breaking)
            if i >= 1:
                pos_minus_2 = sequence[i - 1]
                if pos_minus_2 in "GP":
                    c_score += 0.1

            if c_score > best_c_score:
                best_c_score = c_score
                best_cleavage_pos = i + 2  # 1-indexed position after cleavage

        c_region_score = max(0, min(1.0, best_c_score))

        # Calculate overall confidence
        confidence = (n_region_score * 0.2 + h_region_score * 0.5 + c_region_score * 0.3)

        # Determine if it's likely a signal peptide
        is_signal_peptide = (
            confidence >= 0.5 and
            h_region_score >= 0.3 and
            best_cleavage_pos is not None
        )

        return SignalPeptideResult(
            is_signal_peptide=is_signal_peptide,
            cleavage_position=best_cleavage_pos if is_signal_peptide else None,
            confidence=round(confidence, 3),
            n_region_score=round(n_region_score, 3),
            h_region_score=round(h_region_score, 3),
            c_region_score=round(c_region_score, 3),
            details={
                "n_region": n_region,
                "h_region": h_region,
                "avg_hydrophobicity": round(avg_hydrophobicity, 3),
                "all_cleavage_candidates": best_cleavage_pos,
            }
        ).to_dict()

    def predict_tlr_motifs(
        self,
        sequence: str,
        nucleic_acid: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Predict Toll-like receptor (TLR) ligand motifs.

        Detects motifs recognized by various TLRs:
            - TLR5: Flagellin conserved regions
            - TLR4: LPS-binding regions (basic/hydrophobic patches)
            - TLR9: CpG motifs (if nucleic acid sequence provided)

        Args:
            sequence: Amino acid sequence (one-letter codes).
            nucleic_acid: Optional nucleic acid sequence for CpG detection.

        Returns:
            Dictionary with keys:
                - positions: List of 1-indexed start positions.
                - scores: List of match scores.
                - motifs: List of matched sequences.
                - annotations: List of TLR/motif type annotations.
                - matches: List of MotifMatch objects.
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        matches: List[MotifMatch] = []

        # TLR5: Flagellin N-terminal motifs
        for match in self.flagellin_n_pattern.finditer(sequence):
            matches.append(MotifMatch(
                position=match.start() + 1,
                end=match.end(),
                sequence=match.group(),
                motif_name="flagellin_n_terminal",
                score=0.8,
                annotation="TLR5 ligand - flagellin N-terminal conserved region"
            ))

        # TLR5: Flagellin C-terminal motifs
        for match in self.flagellin_c_pattern.finditer(sequence):
            matches.append(MotifMatch(
                position=match.start() + 1,
                end=match.end(),
                sequence=match.group(),
                motif_name="flagellin_c_terminal",
                score=0.8,
                annotation="TLR5 ligand - flagellin C-terminal conserved region"
            ))

        # TLR4-associated: LPS-binding regions (basic + hydrophobic patches)
        for match in self.lps_binding_pattern.finditer(sequence):
            matches.append(MotifMatch(
                position=match.start() + 1,
                end=match.end(),
                sequence=match.group(),
                motif_name="lps_binding",
                score=0.6,
                annotation="Potential LPS-binding region (basic/hydrophobic patch)"
            ))

        # Detect extended basic patches (potential LPS binding)
        basic_patch_pattern = re.compile(r"[KRH]{3,}")
        for match in basic_patch_pattern.finditer(sequence):
            if len(match.group()) >= 3:
                matches.append(MotifMatch(
                    position=match.start() + 1,
                    end=match.end(),
                    sequence=match.group(),
                    motif_name="basic_patch",
                    score=0.5 + 0.1 * min(len(match.group()) - 3, 3),
                    annotation="Basic patch - potential LPS/nucleic acid binding"
                ))

        # Detect hydrophobic patches
        hydro_patch_pattern = re.compile(r"[AILMFVWY]{4,}")
        for match in hydro_patch_pattern.finditer(sequence):
            matches.append(MotifMatch(
                position=match.start() + 1,
                end=match.end(),
                sequence=match.group(),
                motif_name="hydrophobic_patch",
                score=0.4,
                annotation="Hydrophobic patch - potential membrane interaction"
            ))

        # TLR9: CpG motifs (nucleic acid)
        cpg_matches: List[MotifMatch] = []
        if nucleic_acid:
            nucleic_acid = nucleic_acid.upper().replace(" ", "").replace("\n", "")
            # Optimal CpG motifs: GTCGTT, GACGTT, etc.
            cpg_optimal = re.compile(r"[GAT]TCGTT|G[AT]CGTT|GTCG[AT]T", re.IGNORECASE)
            for match in cpg_optimal.finditer(nucleic_acid):
                cpg_matches.append(MotifMatch(
                    position=match.start() + 1,
                    end=match.end(),
                    sequence=match.group(),
                    motif_name="cpg_optimal",
                    score=0.9,
                    annotation="TLR9 ligand - optimal CpG motif"
                ))

            # General CpG dinucleotides
            cpg_basic = re.compile(r"CG")
            for match in cpg_basic.finditer(nucleic_acid):
                # Check if not already covered by optimal
                already_covered = any(
                    m.position <= match.start() + 1 <= m.end
                    for m in cpg_matches
                )
                if not already_covered:
                    cpg_matches.append(MotifMatch(
                        position=match.start() + 1,
                        end=match.end(),
                        sequence=match.group(),
                        motif_name="cpg_basic",
                        score=0.4,
                        annotation="TLR9 ligand - CpG dinucleotide"
                    ))

        return {
            "positions": [m.position for m in matches],
            "scores": [m.score for m in matches],
            "motifs": [m.sequence for m in matches],
            "annotations": [m.annotation for m in matches],
            "matches": [m.to_dict() for m in matches],
            "cpg_matches": [m.to_dict() for m in cpg_matches] if nucleic_acid else [],
        }

    def predict_immunostim(self, sequence: str) -> Dict[str, Any]:
        """
        Predict immunostimulatory motifs.

        Detects motifs that may enhance immune responses:
            - RGD motif (integrin binding)
            - Basic patches (potential adjuvant activity)
            - Amphipathic helices

        Args:
            sequence: Amino acid sequence (one-letter codes).

        Returns:
            Dictionary with keys:
                - positions: List of 1-indexed start positions.
                - scores: List of motif scores.
                - motifs: List of matched sequences.
                - annotations: List of motif descriptions.
                - matches: List of MotifMatch objects.
        """
        sequence = sequence.upper().replace(" ", "").replace("\n", "")
        matches: List[MotifMatch] = []

        # RGD motif (integrin binding)
        rgd_pattern = re.compile(r"RGD[SFWV]?")
        for match in rgd_pattern.finditer(sequence):
            score = 0.8 if len(match.group()) > 3 else 0.7
            matches.append(MotifMatch(
                position=match.start() + 1,
                end=match.end(),
                sequence=match.group(),
                motif_name="rgd",
                score=score,
                annotation="RGD motif - integrin binding, cell adhesion"
            ))

        # Extended RGD variants
        rgd_variants = re.compile(r"[KR]GD[SFWV]|RG[DN][SFWV]")
        for match in rgd_variants.finditer(sequence):
            # Skip if overlaps with canonical RGD
            pos = match.start() + 1
            if not any(m.position == pos for m in matches if m.motif_name == "rgd"):
                matches.append(MotifMatch(
                    position=pos,
                    end=match.end(),
                    sequence=match.group(),
                    motif_name="rgd_variant",
                    score=0.5,
                    annotation="RGD variant - potential integrin interaction"
                ))

        # Basic patches (potential adjuvant activity)
        basic_patch = re.compile(r"[KR]{4,}")
        for match in basic_patch.finditer(sequence):
            patch_len = len(match.group())
            score = min(0.9, 0.5 + 0.1 * patch_len)
            matches.append(MotifMatch(
                position=match.start() + 1,
                end=match.end(),
                sequence=match.group(),
                motif_name="basic_patch",
                score=score,
                annotation="Basic patch - potential adjuvant activity, membrane disruption"
            ))

        # Amphipathic helix detection (simplified)
        # Look for alternating hydrophobic/hydrophilic patterns
        amphipathic_matches = self._detect_amphipathic_helices(sequence)
        matches.extend(amphipathic_matches)

        # KK, KR, RK, RR dipeptides (cathelin-like motifs)
        cationic_dipep = re.compile(r"[KR][KR]")
        dipep_positions = set()
        for match in cationic_dipep.finditer(sequence):
            pos = match.start() + 1
            if pos not in dipep_positions:
                dipep_positions.add(pos)
                matches.append(MotifMatch(
                    position=pos,
                    end=match.end(),
                    sequence=match.group(),
                    motif_name="cationic_dipeptide",
                    score=0.3,
                    annotation="Cationic dipeptide - antimicrobial peptide feature"
                ))

        # PAMP-mimicking sequences (poly-proline)
        polyproline = re.compile(r"P{3,}")
        for match in polyproline.finditer(sequence):
            matches.append(MotifMatch(
                position=match.start() + 1,
                end=match.end(),
                sequence=match.group(),
                motif_name="polyproline",
                score=0.4,
                annotation="Polyproline stretch - PPII helix, potential PRR interaction"
            ))

        return {
            "positions": [m.position for m in matches],
            "scores": [m.score for m in matches],
            "motifs": [m.sequence for m in matches],
            "annotations": [m.annotation for m in matches],
            "matches": [m.to_dict() for m in matches],
        }

    def _detect_amphipathic_helices(
        self,
        sequence: str,
        min_length: int = 11
    ) -> List[MotifMatch]:
        """
        Detect potential amphipathic helical regions.

        Amphipathic helices have one hydrophobic face and one hydrophilic face.
        In an alpha helix, residues i and i+3/i+4 are on the same face.

        Args:
            sequence: Amino acid sequence.
            min_length: Minimum helix length to consider.

        Returns:
            List of MotifMatch objects for detected amphipathic regions.
        """
        matches: List[MotifMatch] = []

        if len(sequence) < min_length:
            return matches

        # Sliding window for amphipathic helix detection
        window_size = min_length

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i:i + window_size]

            # Calculate hydrophobic moment (simplified)
            # Check if positions 0, 3, 4, 7, 8, 11 are hydrophobic (one face)
            # and positions 1, 2, 5, 6, 9, 10 are hydrophilic (other face)
            face_a = [0, 3, 4, 7, 8]  # One helical face (mod 11)
            face_b = [1, 2, 5, 6, 9, 10]  # Other helical face (mod 11)

            hydrophobic_face_a = sum(
                1 for j in face_a if j < len(window) and window[j] in HYDROPHOBIC
            )
            hydrophilic_face_b = sum(
                1 for j in face_b if j < len(window) and window[j] in (POLAR | BASIC | ACIDIC)
            )

            # Score amphipathicity
            max_a = min(len(face_a), len(window))
            max_b = min(len(face_b), len(window))

            amphipathic_score = 0.0
            if max_a > 0 and max_b > 0:
                amphipathic_score = (
                    (hydrophobic_face_a / max_a) * 0.5 +
                    (hydrophilic_face_b / max_b) * 0.5
                )

            if amphipathic_score >= 0.6:
                matches.append(MotifMatch(
                    position=i + 1,
                    end=i + window_size,
                    sequence=window,
                    motif_name="amphipathic_helix",
                    score=round(amphipathic_score, 3),
                    annotation="Potential amphipathic helix - membrane active"
                ))

        # Merge overlapping matches, keeping highest scoring
        if matches:
            matches = self._merge_overlapping_matches(matches)

        return matches

    def _merge_overlapping_matches(
        self,
        matches: List[MotifMatch]
    ) -> List[MotifMatch]:
        """Merge overlapping motif matches, keeping highest scoring."""
        if not matches:
            return []

        # Sort by position
        sorted_matches = sorted(matches, key=lambda m: m.position)
        merged: List[MotifMatch] = []

        current = sorted_matches[0]
        for next_match in sorted_matches[1:]:
            # Check for overlap
            if next_match.position <= current.end:
                # Keep the higher scoring one
                if next_match.score > current.score:
                    current = next_match
            else:
                merged.append(current)
                current = next_match

        merged.append(current)
        return merged


def predict_all_innate(
    sequence: str,
    nucleic_acid: Optional[str] = None
) -> Dict[str, Any]:
    """
    Run all innate immunity predictions on a sequence.

    Args:
        sequence: Amino acid sequence (one-letter codes).
        nucleic_acid: Optional nucleic acid sequence for CpG detection.

    Returns:
        Dictionary containing results from all prediction methods:
            - n_glycosylation: N-glycosylation site predictions.
            - o_glycosylation: O-glycosylation site predictions.
            - signal_peptide: Signal peptide prediction.
            - tlr_motifs: TLR ligand motif predictions.
            - immunostim: Immunostimulatory motif predictions.
            - summary: Summary statistics across all predictions.
    """
    predictor = InnatePredictor()

    n_glyco = predictor.predict_n_glyco(sequence)
    o_glyco = predictor.predict_o_glyco(sequence)
    signal_pep = predictor.predict_signal_peptide(sequence)
    tlr = predictor.predict_tlr_motifs(sequence, nucleic_acid)
    immunostim = predictor.predict_immunostim(sequence)

    # Generate summary
    summary = {
        "sequence_length": len(sequence.replace(" ", "").replace("\n", "")),
        "n_glyco_sites": len(n_glyco["positions"]),
        "canonical_n_glyco": sum(1 for a in n_glyco["annotations"] if a == "canonical"),
        "o_glyco_sites": len(o_glyco["positions"]),
        "high_confidence_o_glyco": sum(1 for a in o_glyco["annotations"] if a == "high"),
        "has_signal_peptide": signal_pep["is_signal_peptide"],
        "signal_peptide_confidence": signal_pep["confidence"],
        "tlr_motif_count": len(tlr["positions"]),
        "immunostim_motif_count": len(immunostim["positions"]),
        "has_rgd": any("rgd" in m.lower() for m in immunostim.get("motifs", [])),
    }

    return {
        "n_glycosylation": n_glyco,
        "o_glycosylation": o_glyco,
        "signal_peptide": signal_pep,
        "tlr_motifs": tlr,
        "immunostim": immunostim,
        "summary": summary,
    }


# Also provide class-based method that follows EpitopeMethod pattern
class NLinkedGlycosylationMethod(EpitopeMethod):
    """N-linked glycosylation site predictor as EpitopeMethod."""

    name = "N-Glycosylation Sites"
    key = "n_glyco"
    category = MethodCategory.INNATE
    description = "Predicts N-linked glycosylation sites (N-X-S/T sequon)."
    citation = "Marshall RD. Annu Rev Biochem. 1972;41:673-702."
    default_window = 3
    default_threshold = 0.5

    def predict(self, sequence: str, **kwargs) -> MethodResult:
        """Predict N-glycosylation sites."""
        sequence = self.validate_sequence(sequence)
        predictor = InnatePredictor()
        result = predictor.predict_n_glyco(sequence)

        # Convert to MethodResult format
        positions = np.array(result["positions"]) if result["positions"] else np.array([])
        scores = np.array(result["scores"]) if result["scores"] else np.array([])

        return MethodResult(
            scores=scores,
            positions=positions,
            epitopes=[],  # Glyco sites aren't epitopes per se
            metadata={
                "sites": result["sites"],
                "motifs": result["motifs"],
                "annotations": result["annotations"],
            }
        )


class OLinkedGlycosylationMethod(EpitopeMethod):
    """O-linked glycosylation site predictor as EpitopeMethod."""

    name = "O-Glycosylation Sites"
    key = "o_glyco"
    category = MethodCategory.INNATE
    description = "Predicts O-linked glycosylation sites on Ser/Thr."
    citation = "Julenius K et al. Glycobiology. 2005;15(2):153-164."
    default_window = 5
    default_threshold = 0.5

    def predict(self, sequence: str, **kwargs) -> MethodResult:
        """Predict O-glycosylation sites."""
        sequence = self.validate_sequence(sequence)
        predictor = InnatePredictor()
        result = predictor.predict_o_glyco(sequence)

        positions = np.array(result["positions"]) if result["positions"] else np.array([])
        scores = np.array(result["scores"]) if result["scores"] else np.array([])

        return MethodResult(
            scores=scores,
            positions=positions,
            epitopes=[],
            metadata={
                "sites": result["sites"],
                "motifs": result["motifs"],
                "annotations": result["annotations"],
            }
        )


class SignalPeptideMethod(EpitopeMethod):
    """Signal peptide predictor as EpitopeMethod."""

    name = "Signal Peptide"
    key = "signal_peptide"
    category = MethodCategory.INNATE
    description = "Predicts signal peptides and cleavage sites."
    citation = "von Heijne G. Nucleic Acids Res. 1986;14(11):4683-4690."
    default_window = 15
    default_threshold = 0.5

    def predict(self, sequence: str, **kwargs) -> MethodResult:
        """Predict signal peptide."""
        sequence = self.validate_sequence(sequence)
        predictor = InnatePredictor()
        result = predictor.predict_signal_peptide(sequence)

        # For signal peptide, we return scores for the N-terminal region
        cleavage_pos = result.get("cleavage_position")
        if cleavage_pos and result.get("is_signal_peptide"):
            positions = np.arange(1, cleavage_pos + 1)
            scores = np.ones(cleavage_pos) * result["confidence"]
        else:
            positions = np.array([])
            scores = np.array([])

        return MethodResult(
            scores=scores,
            positions=positions,
            epitopes=[],
            metadata={
                "is_signal_peptide": result["is_signal_peptide"],
                "cleavage_position": result["cleavage_position"],
                "confidence": result["confidence"],
                "n_region_score": result["n_region_score"],
                "h_region_score": result["h_region_score"],
                "c_region_score": result["c_region_score"],
            }
        )


# Convenience exports
__all__ = [
    "InnatePredictor",
    "GlycoSite",
    "SignalPeptideResult",
    "MotifMatch",
    "predict_all_innate",
    "NLinkedGlycosylationMethod",
    "OLinkedGlycosylationMethod",
    "SignalPeptideMethod",
    "O_GLYCO_POSITIVE",
    "O_GLYCO_NEGATIVE",
    "HYDROPHOBICITY",
]
