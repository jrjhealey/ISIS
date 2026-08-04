"""
Sequence input parsing for IEDB tools.
Supports FASTA, one-sequence, and whitespace-separated formats.
"""

import re


class SequenceInput:
    """Base class for sequence input parsing."""

    def __init__(self, sequence_text):
        self.sequence_text = sequence_text.strip()
        self.sequences = []
        self.sequence_names = []
        self._parse()

    def _parse(self):
        raise NotImplementedError("Subclasses must implement _parse()")

    def as_amino_acid_text(self):
        """Yield each sequence as amino acid text."""
        for seq in self.sequences:
            yield seq

    @classmethod
    def valid_sequence_text(cls, text):
        raise NotImplementedError("Subclasses must implement valid_sequence_text()")


class FASTASequenceInput(SequenceInput):
    """Parse FASTA format sequences."""

    VALID_AA = set('ACDEFGHIKLMNPQRSTVWYX*')

    def _parse(self):
        lines = self.sequence_text.split('\n')
        current_name = None
        current_seq = []

        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name is not None:
                    self.sequence_names.append(current_name)
                    self.sequences.append(''.join(current_seq).upper())
                current_name = line[1:].strip()
                if not current_name:
                    current_name = f"Sequence_{len(self.sequences) + 1}"
                current_seq = []
            else:
                current_seq.append(line.replace(' ', ''))

        if current_name is not None:
            self.sequence_names.append(current_name)
            self.sequences.append(''.join(current_seq).upper())

    @classmethod
    def valid_sequence_text(cls, text):
        """Check if text is valid FASTA format."""
        text = text.strip()
        if not text:
            return False
        lines = text.split('\n')
        has_header = False
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                has_header = True
            elif has_header:
                cleaned = line.replace(' ', '').upper()
                if not all(c in cls.VALID_AA for c in cleaned):
                    return False
        return has_header


class OneSequenceInput(SequenceInput):
    """Parse a single sequence without header."""

    VALID_AA = set('ACDEFGHIKLMNPQRSTVWYX*')

    def _parse(self):
        seq = ''.join(self.sequence_text.split()).upper()
        self.sequences.append(seq)
        self.sequence_names.append("Sequence_1")

    @classmethod
    def valid_sequence_text(cls, text):
        """Check if text is a valid single amino acid sequence."""
        text = text.strip()
        if not text:
            return False
        if text.startswith('>'):
            return False
        cleaned = ''.join(text.split()).upper()
        if len(cleaned) < 6:
            return False
        return all(c in cls.VALID_AA for c in cleaned)


class WhitespaceSeparatedSequenceInput(SequenceInput):
    """Parse whitespace-separated sequences with optional names."""

    VALID_AA = set('ACDEFGHIKLMNPQRSTVWYX*')

    def _parse(self):
        lines = self.sequence_text.split('\n')
        for i, line in enumerate(lines):
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                name = parts[0]
                seq = ''.join(parts[1:]).upper()
            else:
                name = f"Sequence_{len(self.sequences) + 1}"
                seq = parts[0].upper()

            if seq and all(c in self.VALID_AA for c in seq):
                self.sequence_names.append(name)
                self.sequences.append(seq)

    @classmethod
    def valid_sequence_text(cls, text):
        """Check if text is valid whitespace-separated format."""
        text = text.strip()
        if not text:
            return False
        if text.startswith('>'):
            return False
        lines = text.split('\n')
        valid_lines = 0
        for line in lines:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 2:
                seq = ''.join(parts[1:]).upper()
                if all(c in cls.VALID_AA for c in seq) and len(seq) >= 6:
                    valid_lines += 1
        return valid_lines >= 1
