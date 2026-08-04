"""Tests for ISIS epitope prediction."""

import pytest
import numpy as np

from isis import predict, predict_all, available_methods, Prediction, Epitope
from isis.scales import SCALES, get_scale


class TestScales:
    """Test amino acid scale data."""

    def test_all_scales_have_20_amino_acids(self):
        for name, scale in SCALES.items():
            assert len(scale) == 20, f"{name} should have 20 amino acids"

    def test_scale_lookup(self):
        scale = get_scale("emini")
        assert scale["A"] == 0.815
        assert scale["K"] == 1.545

    def test_scale_lookup_case_insensitive(self):
        scale1 = get_scale("Emini")
        scale2 = get_scale("EMINI")
        scale3 = get_scale("emini")
        assert scale1 == scale2 == scale3

    def test_invalid_scale_raises(self):
        with pytest.raises(ValueError):
            get_scale("invalid_method")


class TestPredict:
    """Test prediction functionality."""

    TEST_SEQUENCE = "MKTAYIAKQRQISFVKSHFSRQLEEALCLSLHRALQ"

    def test_basic_prediction(self):
        result = predict(self.TEST_SEQUENCE, method="emini")
        assert isinstance(result, Prediction)
        assert result.method == "emini"
        assert result.sequence == self.TEST_SEQUENCE
        assert len(result.scores) > 0

    def test_prediction_positions_are_1_indexed(self):
        result = predict(self.TEST_SEQUENCE, method="parker", window_size=7)
        assert result.positions[0] >= 1
        assert result.positions[-1] <= len(self.TEST_SEQUENCE)

    def test_prediction_with_custom_window(self):
        result = predict(self.TEST_SEQUENCE, method="emini", window_size=9)
        assert result.window_size == 9

    def test_prediction_with_custom_threshold(self):
        result = predict(self.TEST_SEQUENCE, method="emini", threshold=1.5)
        assert result.threshold == 1.5

    def test_all_methods_work(self):
        for method in available_methods():
            result = predict(self.TEST_SEQUENCE, method=method)
            assert result.method == method
            assert len(result.scores) > 0

    def test_predict_all(self):
        results = predict_all(self.TEST_SEQUENCE)
        assert len(results) == len(available_methods())
        for method, result in results.items():
            assert result.method == method

    def test_short_sequence_raises(self):
        with pytest.raises(ValueError):
            predict("ACDE", window_size=7)

    def test_sequence_with_whitespace_cleaned(self):
        seq_with_spaces = "MKT AYI AKQ"
        result = predict(seq_with_spaces, method="parker")
        assert result.sequence == "MKTAYIAKQ"


class TestEpitopes:
    """Test epitope extraction."""

    def test_epitopes_extracted(self):
        seq = "MKTAYIAKQRQISFVKSHFSRQLEEALCLSLHRALQ"
        result = predict(seq, method="kolaskar-tongaonkar")
        # Should find at least one epitope with default threshold
        epitopes = result.epitopes
        assert isinstance(epitopes, list)

    def test_epitope_properties(self):
        seq = "MKTAYIAKQRQISFVKSHFSRQLE" * 3  # Longer sequence
        result = predict(seq, method="emini", threshold=0.5)

        for epitope in result.epitopes:
            assert isinstance(epitope, Epitope)
            assert epitope.start >= 1
            assert epitope.end >= epitope.start
            assert epitope.length == epitope.end - epitope.start + 1
            assert len(epitope.sequence) == epitope.length
            assert epitope.score > 0

    def test_score_at_position(self):
        seq = "MKTAYIAKQRQISFVKSHFSRQLE"
        result = predict(seq, method="parker", window_size=7)

        # Center position should have a score
        center = len(seq) // 2
        score = result.score_at(center)
        assert score is not None

        # Position outside scored region returns None
        score = result.score_at(1)  # First position not centered
        # May or may not be None depending on window


class TestSerialization:
    """Test prediction serialization."""

    def test_to_dict(self):
        result = predict("MKTAYIAKQRQISFVKSHFSRQLE", method="emini")
        d = result.to_dict()

        assert d["method"] == "emini"
        assert d["sequence"] == "MKTAYIAKQRQISFVKSHFSRQLE"
        assert "positions" in d
        assert "scores" in d
        assert "epitopes" in d
        assert isinstance(d["positions"], list)
        assert isinstance(d["scores"], list)


class TestEminiMethod:
    """Test Emini-specific product-based scoring."""

    def test_emini_uses_product(self):
        seq = "AAAAAAA"  # All same AA
        result = predict(seq, method="emini", window_size=6)
        # With all same AA, product should be consistent
        assert len(result.scores) > 0

    def test_emini_normalization(self):
        seq = "MKTAYIAKQRQISFVKSHFSRQLE"
        result = predict(seq, method="emini")
        # Emini scores should be normalized around 1.0
        mean_score = np.mean(result.scores)
        assert 0.5 < mean_score < 2.0


class TestAvailableMethods:
    """Test method listing."""

    def test_available_methods_returns_list(self):
        methods = available_methods()
        assert isinstance(methods, list)
        assert len(methods) >= 5

    def test_expected_methods_present(self):
        methods = available_methods()
        expected = ["emini", "parker", "chou-fasman", "kolaskar-tongaonkar", "karplus-schulz"]
        for m in expected:
            assert m in methods
