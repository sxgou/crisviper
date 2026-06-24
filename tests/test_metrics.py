"""Tests for the metrics module (crisviper/metrics.py)."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
from crisviper import (
    effective_alleles, diversity_index,
    alleles_per_cell, singletons_per_cell, carlin_potential,
)


class TestEffectiveAlleles:
    def test_empty(self):
        assert effective_alleles([]) == 0.0

    def test_single_allele(self):
        # Single non-template allele → 2^0 - 1 = 0
        assert effective_alleles([10], is_template=[False]) == 0.0

    def test_two_equal(self):
        # Two equally frequent non-template alleles → H = 1, 2^1 - 1 = 1
        result = effective_alleles([10, 10], is_template=[False, False])
        assert abs(result - 1.0) < 1e-6

    def test_with_template_excluded(self):
        # Template + 2 non-template equally frequent
        result = effective_alleles([10, 10, 10], is_template=[True, False, False])
        # Non-template: [10, 10] → H = 1, 2^1 - 1 = 1
        assert abs(result - 1.0) < 1e-6

    def test_all_template(self):
        result = effective_alleles([10], is_template=[True])
        # No non-template alleles → 0
        assert result == 0.0

    def test_three_unequal(self):
        # Three alleles: 5, 3, 2 — H ≈ 1.4855, 2^H - 1 ≈ 1.8001
        result = effective_alleles([5, 3, 2], is_template=[False, False, False])
        assert abs(result - 1.8001) < 0.001, \
            f"Expected ~1.8001 for frequencies [5,3,2], got {result:.4f}"


class TestDiversityIndex:
    def test_empty(self):
        assert diversity_index([]) == 0.0

    def test_basic(self):
        result = diversity_index([10, 10], is_template=[False, False])
        assert abs(result - 0.05) < 0.001, \
            f"Expected diversity_index=0.05, got {result:.4f}"

    def test_all_template(self):
        assert diversity_index([10], is_template=[True]) == 0.0

    def test_normalize_by_edited(self):
        result = diversity_index([10, 5], is_template=[True, False], normalize_by_edited=True)
        assert result == 0.0, \
            f"Single non-template allele should yield 0.0, got {result}"


class TestAllelesPerCell:
    def test_empty(self):
        assert alleles_per_cell(0, 0) == 0.0

    def test_basic(self):
        assert alleles_per_cell(10, 100) == 0.1


class TestSingletonsPerCell:
    def test_empty(self):
        assert singletons_per_cell(0, 100) == 0.0

    def test_basic(self):
        assert singletons_per_cell(5, 100) == 0.05


class TestCarlinPotential:
    def test_basic(self):
        assert carlin_potential(10, 3) == 7.0

    def test_all_modified(self):
        assert carlin_potential(10, 10) == 0.0

    def test_no_negative(self):
        assert carlin_potential(5, 10) == 0.0
