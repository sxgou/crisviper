"""Tests for the mutation detection module (crisviper/mutation.py)."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
from crisviper import (
    extract_mutations, classify_mutation_type,
    MutationEvent, MutationType, AlignmentStats,
    CutsiteRegion,
)


class TestExtractMutations:
    """extract_mutations — Extract mutation events from alignment results"""

    def test_identical_sequences(self):
        """Identical sequences → no mutations"""
        ar = "ACGTACGT"
        aq = "ACGTACGT"
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 0

    def test_empty_input(self):
        """Empty alignment → empty list"""
        assert extract_mutations("", "") == []

    def test_mismatched_lengths(self):
        """Mismatched lengths → empty list"""
        assert extract_mutations("ACGT", "AC") == []

    def test_single_substitution(self):
        """Single point mutation"""
        ar = "ACGT"
        aq = "ACCT"  # G→C at position 2
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 1
        m = mutations[0]
        assert m.type == MutationType.SUBSTITUTION
        assert m.ref_pos == 2
        assert m.ref_base == "G"
        assert m.query_base == "C"
        assert m.length == 1

    def test_two_substitutions(self):
        """Two independent point mutations"""
        ar = "ACGTACGT"
        aq = "ACCTACAT"  # G→C at pos2, G→A at pos6
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 2
        assert mutations[0].type == MutationType.SUBSTITUTION
        assert mutations[0].ref_pos == 2
        assert mutations[1].type == MutationType.SUBSTITUTION
        assert mutations[1].ref_pos == 6

    def test_single_deletion(self):
        """Single deletion (1bp)"""
        ar = "ACGT"
        aq = "A-GT"  # C deleted
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 1
        m = mutations[0]
        assert m.type == MutationType.DELETION
        assert m.length == 1
        assert m.ref_base == "C"

    def test_multi_base_deletion(self):
        """Multi-base deletion"""
        ar = "ACGTACG"
        aq = "A-----G"  # CGTAC deleted (5bp)
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 1
        m = mutations[0]
        assert m.type == MutationType.DELETION
        assert m.length == 5

    def test_single_insertion(self):
        """Single insertion (1bp)"""
        ar = "ACGT"
        aq = "ACGGT"  # G inserted between C and G
        # alignment: ref=AC-T, query=ACGT → wait, this is wrong
        # For insertion: ref has gap, query has base
        # Let me construct a clearer case
        ar = "AC--T"
        aq = "ACGGT"  # GG inserted between C and T
        mutations = extract_mutations(ar, aq)
        assert len(mutations) >= 1
        # At least one insertion event
        has_ins = any(m.type == MutationType.INSERTION for m in mutations)
        assert has_ins

    def test_missing_primer5_region_no_mutation(self):
        """Primer5 region matches should not be detected as mutations"""
        # Simulating full-length alignment where primer regions are identical
        # and mutation is in internal region
        ar = "TATGTGTGGGAGGGCTAAGAGG" + "ACGTACGT" + "TAGTTGCCAGCCATCTGTTGT"
        aq = "TATGTGTGGGAGGGCTAAGAGG" + "ACGTACAT" + "TAGTTGCCAGCCATCTGTTGT"
        # Mutation at internal region (ref: G→A at internal pos 5)
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 1
        assert mutations[0].type == MutationType.SUBSTITUTION

    def test_cutsite_window_detection(self):
        """Mutation inside cutsite window → in_cutsite_window=True"""
        ar = "ACGTACGTACGT"
        aq = "ACGTACGAACGT"  # T→A at pos 7
        cutsites = [CutsiteRegion("T1", 5, 10)]
        mutations = extract_mutations(ar, aq, cutsites=cutsites, sub_window=3)
        # pos 7 is inside cutsite 5-10 → in_window should be True
        for m in mutations:
            if m.ref_pos == 7:
                assert m.in_cutsite_window

    def test_outside_cutsite_window(self):
        """Mutation outside cutsite window → in_cutsite_window=False"""
        ar = "ACGTACGTACGT"
        aq = "ACGTACGAACGT"  # T→A at pos 7
        cutsites = [CutsiteRegion("T1", 1, 3)]
        mutations = extract_mutations(ar, aq, cutsites=cutsites, sub_window=1)
        # pos 7 is outside cutsite+1 → in_window should be False
        for m in mutations:
            if m.ref_pos == 7:
                assert not m.in_cutsite_window

    def test_complex_indel_adjacent(self):
        """Adjacent insertion + deletion merged into a complex event"""
        # This is an edge case that may or may not merge
        # Test that complex events can form
        ar = "ACGT--AC"
        aq = "AC--GGAC"  # DEL at pos 3-4 (GT), INS at pos 5-6 (GG)
        mutations = extract_mutations(ar, aq)
        # May be merged or separate depending on adjacency calculation
        assert len(mutations) >= 1
        # Should at least detect both indels somehow
        types = {m.type for m in mutations}
        assert MutationType.DELETION in types or MutationType.INDEL in types

    def test_ref_pos_mapping(self):
        """ref_pos correctly points to reference sequence coordinates"""
        # With gaps in the alignment:
        # Ref:  A C G T A C G T
        # Query:A C - T A C G T
        # Position:0 1 2 3 4 5 6 7
        # Query has deletion at pos2 (C deleted)
        ar = "ACGTACGT"
        aq = "AC-TACGT"
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 1
        m = mutations[0]
        assert m.type == MutationType.DELETION
        assert m.ref_pos == 2  # C at ref position 2


class TestGreedyIndelMerge:
    """Tests for the greedy INDEL merge algorithm in _merge_adjacent_indels"""

    def test_adjacent_del_plus_ins_merged(self):
        """Adjacent deletion and insertion merge into single INDEL"""
        # Ref:  A C G T - - A C
        # Query: A C - - G G A C
        # DEL(2,2) + INS(4,2) adjacent -> merged INDEL
        ar = "ACGT--AC"
        aq = "AC--GGAC"
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 1
        assert mutations[0].type == MutationType.INDEL

    def test_separate_del_and_ins_not_merged(self):
        """Non-adjacent deletion and insertion remain as separate INDEL events"""
        # Ref:  A C G T A C - - T
        # Query: A C - - A C G G T
        # DEL(2,2) at ref_pos=2, INS(6,2) at ref_pos=6
        # >1bp apart, so not merged — but each converted to INDEL individually
        ar = "ACGTAC--T"
        aq = "AC--ACGGT"
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 2
        assert all(m.type == MutationType.INDEL for m in mutations)

    def test_multi_event_group_merged(self):
        """DEL + INS + SUB in one adjacent group merge into single INDEL"""
        # Ref:  A C G T - - A C
        # Query: A C - - G G T C
        ar = "ACGT--AC"
        aq = "AC--GGTC"
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 1
        assert mutations[0].type == MutationType.INDEL

    def test_pure_substitutions_not_merged(self):
        """Adjacent substitutions remain as individual events"""
        # Ref:  A C G T A C
        # Query: A C C T T C
        ar = "ACGTAC"
        aq = "ACCTTC"
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 2
        assert all(m.type == MutationType.SUBSTITUTION for m in mutations)

    def test_single_event_unchanged(self):
        """Single deletion event remains as-is"""
        ar = "ACGTAC"
        aq = "AC--AC"
        mutations = extract_mutations(ar, aq)
        assert len(mutations) == 1
        assert mutations[0].type == MutationType.DELETION


class TestClassifyMutationType:
    """classify_mutation_type — Classify mutation type based on statistics"""

    def test_unmutated(self):
        stats = AlignmentStats(mismatches=0, gaps_in_ref=0, gaps_in_query=0)
        assert classify_mutation_type(stats) == "unmutated"

    def test_only_substitution(self):
        stats = AlignmentStats(mismatches=3, gaps_in_ref=0, gaps_in_query=0)
        assert classify_mutation_type(stats) == "only_substitution"

    def test_only_deletion(self):
        stats = AlignmentStats(mismatches=0, gaps_in_ref=0, gaps_in_query=5)
        assert classify_mutation_type(stats) == "only_deletion"

    def test_only_insertion(self):
        stats = AlignmentStats(mismatches=0, gaps_in_ref=3, gaps_in_query=0)
        assert classify_mutation_type(stats) == "only_insertion"

    def test_insertion_and_deletion(self):
        stats = AlignmentStats(mismatches=0, gaps_in_ref=2, gaps_in_query=3)
        assert classify_mutation_type(stats) == "insertion_and_deletion"

    def test_all_three(self):
        stats = AlignmentStats(mismatches=1, gaps_in_ref=2, gaps_in_query=3)
        assert classify_mutation_type(stats) == "insertion_deletion_substitution"


class TestBuildMutationSummary:
    """build_mutation_summary — Aggregate statistics"""

    def test_empty(self):
        from crisviper import build_mutation_summary
        summary = build_mutation_summary([])
        assert summary["type_counts"] == {}
        assert summary["ins_lengths"] == []
        assert summary["total_mismatches"] == 0

    def test_single_mutated(self):
        from crisviper import build_mutation_summary, AlignmentResult, QueryRecord, AlignmentStats
        result = AlignmentResult(
            query=QueryRecord(readName="test", readCount=5, seq="ACGT"),
            success=True,
            stats=AlignmentStats(mismatches=1, gaps_in_ref=0, gaps_in_query=0, gap_blocks_ref=[], gap_blocks_query=[]),
        )
        summary = build_mutation_summary([result])
        assert summary["type_counts"]["only_substitution"]["sequences"] == 1
        assert summary["type_counts"]["only_substitution"]["reads"] == 5
        assert summary["total_mismatches"] == 1

    def test_skips_failed_results(self):
        from crisviper import build_mutation_summary, AlignmentResult, QueryRecord
        result = AlignmentResult(
            query=QueryRecord(readName="test", seq="ACGT"),
            success=False,
            error="alignment failed",
        )
        summary = build_mutation_summary([result])
        assert summary["type_counts"] == {}


# ═══════════════════════════════════════════════════════════════
# HGVS annotation tests
# ═══════════════════════════════════════════════════════════════

class TestAnnotateMutation:
    def test_substitution(self):
        m = MutationEvent(type=MutationType.SUBSTITUTION, ref_pos=40, ref_base="G", query_base="T")
        from crisviper.mutation import annotate_mutation
        assert annotate_mutation(m) == "41G>T"

    def test_deletion_single(self):
        m = MutationEvent(type=MutationType.DELETION, ref_pos=40, length=1)
        from crisviper.mutation import annotate_mutation
        assert annotate_mutation(m) == "41del"

    def test_deletion_multi(self):
        m = MutationEvent(type=MutationType.DELETION, ref_pos=40, length=3)
        from crisviper.mutation import annotate_mutation
        assert annotate_mutation(m) == "41_43del"

    def test_insertion_short(self):
        m = MutationEvent(type=MutationType.INSERTION, ref_pos=40, length=3, query_base="---XYZ")
        from crisviper.mutation import annotate_mutation
        assert annotate_mutation(m) == "41_42ins3"
        assert annotate_mutation(m, full=True) == "41_42insXYZ"

    def test_complex_short(self):
        m = MutationEvent(type=MutationType.INDEL, ref_pos=40, length=5, query_base="---XYZ")
        from crisviper.mutation import annotate_mutation
        assert annotate_mutation(m) == "41_45delins3"
        assert annotate_mutation(m, full=True) == "41_45delinsXYZ"

    def test_multiple_mutations(self):
        from crisviper.mutation import annotate_mutations
        m1 = MutationEvent(type=MutationType.SUBSTITUTION, ref_pos=40, ref_base="G", query_base="T")
        m2 = MutationEvent(type=MutationType.DELETION, ref_pos=100, length=5)
        result = annotate_mutations([m1, m2])
        assert result == "41G>T;101_105del"

    def test_empty_list(self):
        from crisviper.mutation import annotate_mutations
        assert annotate_mutations([]) == ""



class TestClassifyBpEvent:
    """classify_bp_event — per-base event classification"""

    def test_no_events(self):
        from crisviper.mutation import classify_bp_event
        bp, de, ins = classify_bp_event([], 10)
        assert bp == ['N'] * 10
        assert de == [False] * 10
        assert ins == [False] * 10

    def test_substitution(self):
        from crisviper.mutation import classify_bp_event, MutationEvent, MutationType
        ev = MutationEvent(type=MutationType.SUBSTITUTION, ref_pos=5, ref_base='G', query_base='T', length=1)
        bp, de, ins = classify_bp_event([ev], 10)
        assert bp[5] == 'I'
        assert de[5] is True
        assert ins[5] is True

    def test_deletion(self):
        from crisviper.mutation import classify_bp_event, MutationEvent, MutationType
        ev = MutationEvent(type=MutationType.DELETION, ref_pos=3, length=4)
        bp, de, ins = classify_bp_event([ev], 10)
        for i in range(3, 7):
            assert bp[i] == 'D'
            assert de[i] is True
        assert not any(ins)

    def test_insertion(self):
        from crisviper.mutation import classify_bp_event, MutationEvent, MutationType
        ev = MutationEvent(type=MutationType.INSERTION, ref_pos=2, query_base='XYZ', length=3)
        bp, de, ins = classify_bp_event([ev], 10)
        assert ins[2] is True
        # insertion extends to min(2+3-1, ...)
        assert bp[2] == 'I'

    def test_complex(self):
        from crisviper.mutation import classify_bp_event, MutationEvent, MutationType
        ev = MutationEvent(type=MutationType.INDEL, ref_pos=3, query_base='XYZ', length=5)
        bp, de, ins = classify_bp_event([ev], 10)
        # deletion span
        assert any(de[3:8])
        # insertion at start
        assert ins[3] is True

