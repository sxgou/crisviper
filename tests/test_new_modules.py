"""Tests for denoiser, threshold, and caller modules."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pytest

from crisviper.denoiser import directional_adjacency_top_down_denoiser
from crisviper.threshold import compute_threshold
from crisviper.caller import (
    call_alleles_coarse_grain, call_alleles_exact,
    CalledAllele, _event_structure,
)
from crisviper.models import (
    AlignmentResult, AlignmentStats, QueryRecord,
    MutationEvent, MutationType,
)


# ═══════════════════════════════════════════════════════════════
# denoiser tests
# ═══════════════════════════════════════════════════════════════

class TestDirectionalAdjacencyDenoiser:
    def test_single_tag(self):
        tags = ["AAAA"]
        weights = np.array([10])
        result = directional_adjacency_top_down_denoiser(tags, weights)
        assert result == [0]

    def test_two_tags_hamming_1(self):
        tags = ["AAAA", "AAAT"]
        weights = np.array([10, 5])
        result = directional_adjacency_top_down_denoiser(tags, weights)
        assert result[1] == result[0], "HD=1 child should map to parent"

    def test_two_tags_hamming_2(self):
        tags = ["AAAA", "AATT"]
        weights = np.array([10, 5])
        result = directional_adjacency_top_down_denoiser(tags, weights)
        assert result[1] != result[0], "HD=2 should NOT cluster"
        assert result[1] == 1, "HD=2 should be its own parent"

    def test_weight_threshold_no_cluster(self):
        """Child too close in weight to parent → no cluster."""
        # parent=10 can only absorb child with weight ≤ (10+1)/2 = 5
        tags = ["AAAA", "AAAT"]
        weights = np.array([10, 6])  # 6 > 5 → no cluster
        result = directional_adjacency_top_down_denoiser(tags, weights)
        assert result[1] != result[0], "Child too close in weight should not cluster"

    def test_different_lengths_no_cluster(self):
        tags = ["AAAA", "AAAAT"]
        weights = np.array([10, 5])
        result = directional_adjacency_top_down_denoiser(tags, weights)
        assert result[1] == 1, "Different lengths should not cluster"

    def test_exclude_mask(self):
        tags = ["AAAA", "AAAT", "AAAC"]
        weights = np.array([10, 8, 6])
        exclude = np.array([False, True, False])
        result = directional_adjacency_top_down_denoiser(tags, weights, exclude)
        assert result[1] == 1, "Excluded tag should be its own parent"

    def test_three_tag_chain(self):
        """AAAA → AAAT (HD=1) → AATT (HD=1 from AAAT, HD=2 from AAAA).

        AAAA absorbs AAAT. AAAT (already absorbed) also absorbs AATT
        via transitive clustering (weight 3 ≤ (5+1)/2 = 3).
        """
        tags = ["AAAA", "AAAT", "AATT"]
        weights = np.array([10, 5, 3])
        result = directional_adjacency_top_down_denoiser(tags, weights)
        assert result[1] == result[0]
        assert result[2] == result[0], "Transitive clustering via AAAT"

    def test_descending_order_preserves_high_weight(self):
        tags = ["AAAA", "AAAT", "AAAG", "AACA"]
        weights = np.array([10, 5, 4, 3])
        result = directional_adjacency_top_down_denoiser(tags, weights)
        # AAAA (10) can absorb <=5: AAAT (5, HD=1) and AAAG (4, HD=1)
        # AACA (3, HD=1 from AAAA) also absorbed
        assert result[1] == result[0]

    def test_weight_edge_case(self):
        """Boundary: weight 10 can absorb weight 5 (10 >= 2*5-1 = 9)."""
        tags = ["AAAA", "AAAT"]
        weights = np.array([10, 5])
        result = directional_adjacency_top_down_denoiser(tags, weights)
        assert result[1] == result[0]


# ═══════════════════════════════════════════════════════════════
# threshold tests
# ═══════════════════════════════════════════════════════════════

class TestComputeThreshold:
    def test_basic_umi_threshold(self):
        freqs = np.array([100, 80, 60, 40, 20, 10, 5, 3, 2, 1])
        result = compute_threshold(freqs, max_elem=5, n_reads=1000, threshold_type="UMI")
        assert "max_molecules" in result
        assert result["chosen"] >= result["one_tenth_99_pctl"]
        assert result["chosen"] >= result["read_floor"]

    def test_override(self):
        freqs = np.array([100, 80, 60])
        result = compute_threshold(freqs, max_elem=3, n_reads=500,
                                    read_override=50, threshold_type="UMI")
        assert result["chosen"] == 50

    def test_cb_type(self):
        freqs = np.array([100, 80, 60])
        result = compute_threshold(freqs, max_elem=3, n_reads=500,
                                    threshold_type="CB")
        assert "max_cells" in result

    def test_empty_edge(self):
        freqs = np.array([10])
        result = compute_threshold(freqs, max_elem=1, n_reads=10, threshold_type="UMI")
        assert result["chosen"] >= 1


# ═══════════════════════════════════════════════════════════════
# caller tests
# ═══════════════════════════════════════════════════════════════

def make_result(
    mutations: list,
    read_count: int = 10,
    aligned_ref: str = "ACGT",
    aligned_query: str = "ACGT",
    success: bool = True,
) -> AlignmentResult:
    stats = AlignmentStats(
        matches=sum(1 for a, b in zip(aligned_ref, aligned_query) if a == b and a != '-'),
        mismatches=0,
        gaps_in_ref=aligned_ref.count("-"),
        gaps_in_query=aligned_query.count("-"),
    )
    return AlignmentResult(
        query=QueryRecord(readName="test", readCount=read_count, seq=aligned_query.replace("-", "")),
        success=success,
        score=10.0,
        aligned_ref=aligned_ref,
        aligned_query=aligned_query,
        stats=stats,
        mutations=mutations,
    )


class TestEventStructure:
    def test_wt(self):
        assert _event_structure([]) == "WT"

    def test_single_substitution(self):
        m = MutationEvent(type=MutationType.SUBSTITUTION, ref_pos=5, ref_base="A", query_base="T")
        assert _event_structure([m]) == "S@5"

    def test_single_deletion(self):
        m = MutationEvent(type=MutationType.DELETION, ref_pos=3, length=5)
        assert _event_structure([m]) == "D5@3"

    def test_single_insertion(self):
        m = MutationEvent(type=MutationType.INSERTION, ref_pos=2, length=3)
        assert _event_structure([m]) == "I3@2"

    def test_multiple_sorted_by_position(self):
        m1 = MutationEvent(type=MutationType.DELETION, ref_pos=10, length=2)
        m2 = MutationEvent(type=MutationType.SUBSTITUTION, ref_pos=3, ref_base="G", query_base="T")
        result = _event_structure([m1, m2])
        assert result == "S@3+D2@10"


def _make_wt_result(rc=10):
    return make_result([], read_count=rc, aligned_ref="ACGTACGT", aligned_query="ACGTACGT")


def _make_del3_result(rc=10):
    m = MutationEvent(type=MutationType.DELETION, ref_pos=3, length=3)
    return make_result(
        [m], read_count=rc,
        aligned_ref="ACGT---CGT",
        aligned_query="ACG---CGT",
    )


def _make_sub_result(rc=5):
    m = MutationEvent(type=MutationType.SUBSTITUTION, ref_pos=4, ref_base="G", query_base="T")
    return make_result(
        [m], read_count=rc,
        aligned_ref="ACGTCGT",
        aligned_query="ACTTCGT",
    )


class TestCallAllelesCoarseGrain:
    def test_all_wt(self):
        results = [_make_wt_result(10), _make_wt_result(5)]
        alleles = call_alleles_coarse_grain(results)
        assert len(alleles) == 1
        assert alleles[0].event_structure == "WT"
        assert alleles[0].weight == 15

    def test_mixed_dominant_wt(self):
        results = [_make_wt_result(10), _make_del3_result(4)]
        alleles = call_alleles_coarse_grain(results)
        assert len(alleles) == 1  # WT dominates 10/14 > 0.5

    def test_no_dominant(self):
        results = [_make_wt_result(5), _make_del3_result(5)]
        alleles = call_alleles_coarse_grain(results)
        assert len(alleles) == 0  # neither > 50%

    def test_empty_input(self):
        assert call_alleles_coarse_grain([]) == []

    def test_failed_results_skipped(self):
        failed = make_result([], read_count=10, success=False)
        ok = _make_wt_result(10)
        alleles = call_alleles_coarse_grain([failed, ok])
        assert len(alleles) == 1


class TestCallAllelesExact:
    def test_identical_sequences(self):
        results = [_make_wt_result(10), _make_wt_result(5)]
        alleles = call_alleles_exact(results)
        assert len(alleles) == 1

    def test_mixed_dominant(self):
        results = [_make_wt_result(10), _make_sub_result(4)]
        alleles = call_alleles_exact(results)
        assert len(alleles) == 1
        assert alleles[0].event_structure == "WT"

    def test_no_dominant(self):
        results = [_make_wt_result(5), _make_sub_result(6)]
        alleles = call_alleles_exact(results)
        # 6/11 ≈ 0.545 > 0.5, so sub is dominant
        assert len(alleles) == 1
        assert alleles[0].event_structure == "S@4"

    def test_empty_input(self):
        assert call_alleles_exact([]) == []
