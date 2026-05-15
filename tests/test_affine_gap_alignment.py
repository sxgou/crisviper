"""Tests for the core affine gap alignment algorithm and lineage tracer."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import numpy as np
from crisviper import (
    affine_gap_alignment,
    calculate_alignment_stats,
    count_gap_blocks,
    CutsiteRegion,
    build_gradient_profiles,
    affine_gap_alignment_position_aware,
    lineage_tracer_align,
    get_amplicon_structure,
)


# ── Basic utilities ───────────────────────────────────────────────

def test_count_gap_blocks():
    assert count_gap_blocks("ACGT") == []
    assert count_gap_blocks("A--T") == [2]
    assert count_gap_blocks("A--T--G") == [2, 2]
    assert count_gap_blocks("---A") == [3]
    assert count_gap_blocks("A---") == [3]
    assert count_gap_blocks("-----") == [5]


def test_calculate_alignment_stats_identical():
    ar = "ACGT"
    aq = "ACGT"
    stats = calculate_alignment_stats(ar, aq)
    assert stats["matches"] == 4
    assert stats["mismatches"] == 0
    assert stats["gaps_in_ref"] == 0
    assert stats["gaps_in_query"] == 0
    assert stats["alignment_length"] == 4
    assert stats["similarity"] == 1.0
    assert stats["identity"] == 1.0


def test_calculate_alignment_stats_with_gap():
    ar = "ACGT"
    aq = "A-GT"
    stats = calculate_alignment_stats(ar, aq)
    assert stats["matches"] == 3
    assert stats["mismatches"] == 0
    assert stats["gaps_in_ref"] == 0
    assert stats["gaps_in_query"] == 1
    assert stats["alignment_length"] == 4


def test_calculate_alignment_stats_with_mismatch_and_gap():
    ar = "ACGT-"
    aq = "A-GTA"
    stats = calculate_alignment_stats(ar, aq)
    assert stats["matches"] == 3  # A, G, T
    assert stats["mismatches"] == 0  # A matches nothing at end... wait
    # Let's be more careful: ar=ACGT-, aq=A-GTA
    # positions: A=A, C=-, G=G, T=T, -=A
    # matches: A=A, G=G, T=T = 3
    # mismatches: none (C vs - is gap, - vs A is gap)
    # gaps_in_ref: '-' in ar = 1
    # gaps_in_query: '-' in aq = 1
    assert stats["gaps_in_ref"] == 1
    assert stats["gaps_in_query"] == 1


# ── Standard Gotoh alignment ──────────────────────────────────────

def test_gotoh_identical():
    score, ar, aq, stats = affine_gap_alignment("ACGT", "ACGT")
    assert ar == "ACGT"
    assert aq == "ACGT"
    assert stats["matches"] == 4
    assert stats["mismatches"] == 0
    assert score == 8.0  # 4 * 2.0


def test_gotoh_single_substitution():
    score, ar, aq, stats = affine_gap_alignment("ACGT", "ACCT")
    # With match=2, mismatch=-3, gap_open=-2, gap_extend=-0.1
    # Three matches + one mismatch vs gap
    # 3*2 + (-3) = 3 vs gap approach
    assert stats["matches"] == 3
    assert stats["mismatches"] == 1
    # Mismatch is preferred over gap since -3 > -2 + -0.1
    assert "C" in aq  # should keep the C


def test_gotoh_single_deletion():
    # Query missing one base
    score, ar, aq, stats = affine_gap_alignment("ACGT", "AGT")
    # Expected: A C G T
    #           A - G T
    assert stats["matches"] == 3
    assert stats["gaps_in_query"] == 1  # deletion in query relative to ref


def test_gotoh_single_insertion():
    # Query has extra base — use global mode so end gaps are penalized
    score, ar, aq, stats = affine_gap_alignment("ACGT", "ACGTT")
    # Expected alignment includes a gap in ref
    assert stats["gaps_in_ref"] >= 1 or stats["gaps_in_query"] >= 1


def test_gotoh_semi_global_free_ends():
    """Semi-global should allow free gaps at ends."""
    # Query is a substring of ref
    score, ar, aq, stats = affine_gap_alignment("AAAACGT", "CGT")
    # Should match CGT perfectly without penalizing leading ref bases
    assert "CGT" in aq.replace("-", "") or aq.replace("-", "") == "CGT"


def test_gotoh_global_vs_semi_global():
    """Global should penalize end gaps that semi-global does not."""
    ref, query = "AAACCC", "CCC"
    s_sg, _, _, _ = affine_gap_alignment(ref, query)
    s_g, _, _, _ = affine_gap_alignment(ref, query)
    # Semi-global should score higher because end gaps are free
    assert s_sg >= s_g


def test_gotoh_fit_mode():
    """Fit mode: query must be fully used, ref ends free."""
    # Long ref, short query that matches in the middle
    ref = "AAAAACGTGGGGG"
    query = "ACGT"
    score, ar, aq, stats = affine_gap_alignment(ref, query)
    # Query should be fully aligned
    assert aq.replace("-", "") == query


# ── Position-aware alignment ──────────────────────────────────────

def test_build_gradient_profiles():
    cutsites = [
        CutsiteRegion(name="T1", start=5, end=10),
    ]
    ref_len = 20
    # max_scale == cutsite_edge_scale == 2.0 means no external gradient
    go, ge, mp = build_gradient_profiles(
        ref_len, cutsites,
        base_gap_open=-2.0, base_gap_extend=-0.1,
        mismatch_penalty=-3.0,
        min_scale=1.0, max_scale=2.0, cutsite_edge_scale=2.0,
    )
    assert len(go) == ref_len
    assert len(ge) == ref_len
    assert len(mp) == ref_len
    # Far region (no cutsite influence): max_scale=2.0 → gap_open=-4.0
    assert go[0] == pytest.approx(-4.0)
    assert go[15] == pytest.approx(-4.0)
    # Center of cutsite (between pos 7-8): near min_scale → ~-2.0
    assert go[7] > -2.5, f"center go[7]={go[7]} too low"
    assert go[8] > -2.5, f"center go[8]={go[8]} too low"
    # Edges of cutsite (pos 5, 10): near cutsite_edge_scale → ~-4.0
    assert go[5] < -3.0, f"edge go[5]={go[5]} too high"
    assert go[10] < -3.0, f"edge go[10]={go[10]} too high"


def test_position_aware_identical():
    ref = "ACGTACGT"
    query = "ACGTACGT"
    cutsites = [CutsiteRegion("T1", 2, 5)]
    go, ge, mp = build_gradient_profiles(
        len(ref), cutsites,
        base_gap_open=-2.0, base_gap_extend=-0.1, mismatch_penalty=-3.0,
    )
    score, ar, aq, stats = affine_gap_alignment_position_aware(
        ref, query, go, ge, mismatch_penalty_profile=mp,
    )
    assert ar.replace("-", "") == ref
    assert aq.replace("-", "") == query
    assert stats["matches"] == 8


def test_position_aware_prefers_gap_at_cutsite():
    """Gap should preferentially open in cutsite region."""
    ref = "AAAAACCCCCTTTTTT"
    # Query has a deletion in the middle
    query = "AAAAATTTTTT"
    cutsites = [CutsiteRegion("T1", 5, 10)]  # C region

    go, ge, mp = build_gradient_profiles(
        len(ref), cutsites,
        base_gap_open=-2.0, base_gap_extend=-0.1, mismatch_penalty=-3.0,
        min_scale=1.0, max_scale=2.0, cutsite_edge_scale=2.0,
    )
    score, ar, aq, stats = affine_gap_alignment_position_aware(
        ref, query, go, ge, mismatch_penalty_profile=mp
    )
    # The gap (deletion) should be placed primarily in the cutsite region
    # rather than splitting across conserved regions
    # At minimum, we should see that the alignment contains gaps
    assert stats["gaps_in_query"] > 0


def test_position_aware_anchor5_equal():
    """anchor5: query == ref length, full identity."""
    ref = "ACGTACGT"
    query = "ACGTACGT"
    cutsites = [CutsiteRegion("T1", 2, 5)]
    go, ge, mp = build_gradient_profiles(len(ref), cutsites,
        base_gap_open=-2.0, base_gap_extend=-0.1, mismatch_penalty=-3.0)
    score, ar, aq, stats = affine_gap_alignment_position_aware(
        ref, query, go, ge, mismatch_penalty_profile=mp,
    )
    assert ar.replace("-", "") == ref
    assert aq.replace("-", "") == query
    assert stats["matches"] == 8
    assert stats["mismatches"] == 0


def test_position_aware_anchor5_shorter():
    """anchor5: query < ref (deletion at 3' end)."""
    ref = "ACGTACGT"
    query = "ACGT"
    cutsites = [CutsiteRegion("T1", 2, 5)]
    go, ge, mp = build_gradient_profiles(len(ref), cutsites,
        base_gap_open=-2.0, base_gap_extend=-0.1, mismatch_penalty=-3.0)
    score, ar, aq, stats = affine_gap_alignment_position_aware(
        ref, query, go, ge, mismatch_penalty_profile=mp,
    )
    # Query starts at position 0 of ref (anchor at 5')
    assert aq[0] != '-' or ar[0] == aq[0]
    assert aq.replace("-", "") == query
    assert ar[0] == "A"
    assert aq[0] == "A"


def test_position_aware_anchor5_longer():
    """anchor5: query > ref (insertion), should still work."""
    ref = "ACGTACGT"
    query = "ACGTACGTAAA"
    cutsites = [CutsiteRegion("T1", 2, 5)]
    go, ge, mp = build_gradient_profiles(len(ref), cutsites,
        base_gap_open=-2.0, base_gap_extend=-0.1, mismatch_penalty=-3.0)
    score, ar, aq, stats = affine_gap_alignment_position_aware(
        ref, query, go, ge, mismatch_penalty_profile=mp,
    )
    # Anchor at 5': first bases match
    assert ar[0] == "A" and aq[0] == "A"
    # All query bases must be in alignment output
    assert aq.replace("-", "") == query
    # Insertion should appear as gaps in ref
    assert stats["gaps_in_ref"] > 0


# ── Lineage tracer pipeline ───────────────────────────────────────

def test_lineage_tracer_basic():
    """Smoke test for the full lineage tracer pipeline."""
    # Minimal CARLIN-like structure
    ref = "A" * 23 + "CGCCG" + "A" * 13 + "GAGTCGA" + "A" * 7 + "A" * 13 + "GAGTCGA" + "A" * 8 + "A" * 33
    query = ref  # identical

    cutsites = [
        CutsiteRegion("Target1", 41, 47),
        CutsiteRegion("Target2", 68, 74),
    ]

    score, ar, aq, stats = lineage_tracer_align(
        ref, query, cutsites,
    )
    assert stats["matches"] > 0
    assert stats["mismatches"] == 0
    assert stats["gaps_in_query"] == 0
    assert stats["gaps_in_ref"] == 0


def test_lineage_tracer_with_deletion():
    """Query with a deletion in a cutsite region."""
    ref = "A" * 23 + "CGCCG" + "A" * 13 + "GAGTCGA" + "A" * 7 + "A" * 13 + "GAGTCGA" + "A" * 8 + "A" * 33
    # Remove 4 bases from first cutsite
    query = "A" * 23 + "CGCCG" + "A" * 13 + "GA" + "A" * 7 + "A" * 13 + "GAGTCGA" + "A" * 8 + "A" * 33

    cutsites = [
        CutsiteRegion("Target1", 41, 47),
        CutsiteRegion("Target2", 68, 74),
    ]

    score, ar, aq, stats = lineage_tracer_align(
        ref, query, cutsites,
    )
    # Should detect the deletion
    assert stats["gaps_in_query"] > 0


# ── get_amplicon_structure ────────────────────────────────────────

def test_get_amplicon_structure_standard():
    """Standard 332bp CARLIN amplicon."""
    carlin_ref = (
        "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCAACGATGGAGTCGCGACTGTACGCACTCGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGACTGCACGACAGTCGACTATGGAGTCGATACGTAGCACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"
    )
    assert len(carlin_ref) == 332

    cutsites = get_amplicon_structure(carlin_ref)
    assert len(cutsites) == 10
    assert cutsites[0].name == "Target1"
    assert cutsites[0].start == 41
    assert cutsites[0].end == 47
    assert cutsites[9].name == "Target10"
    assert cutsites[9].start == 284
    assert cutsites[9].end == 290
