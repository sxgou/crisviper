"""Tests for summary.py helper functions (_allele_label, _split_indel, etc.)."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from crisviper.summary import (
    _allele_label,
    _split_indel,
    _mutation_type_label,
    _mutation_overlaps_target,
)


# ═══════════════════════════════════════════════════════════════
# _allele_label tests
# ═══════════════════════════════════════════════════════════════

class TestAlleleLabel:
    def test_empty_mutations_returns_wt(self):
        """Empty mutation list should return 'wt'."""
        assert _allele_label([]) == "wt"

    def test_single_substitution(self):
        """Single substitution should produce 'sub:posA>T' format."""
        mutations = [{"type": "substitution", "ref_pos": 5, "ref_base": "A", "query_base": "T"}]
        assert _allele_label(mutations) == "sub:5A>T"

    def test_single_deletion(self):
        """Single deletion should produce 'del:pos-3bp' format."""
        mutations = [{"type": "deletion", "ref_pos": 10, "length": 3}]
        assert _allele_label(mutations) == "del:10-3bp"

    def test_single_insertion_with_seq(self):
        """Single insertion with sequence should produce 'ins:pos+SEQ' format."""
        mutations = [{"type": "insertion", "ref_pos": 7, "query_base": "GGG"}]
        assert _allele_label(mutations) == "ins:7+GGG"

    def test_single_insertion_no_seq(self):
        """Single insertion without query_base should fall back to length."""
        mutations = [{"type": "insertion", "ref_pos": 7, "length": 3}]
        assert _allele_label(mutations) == "ins:7+3bp"

    def test_multiple_mutations_sorted_by_position(self):
        """Multiple mutations should be sorted by ref_pos and joined by ';'."""
        mutations = [
            {"type": "substitution", "ref_pos": 15, "ref_base": "C", "query_base": "G"},
            {"type": "deletion", "ref_pos": 3, "length": 2},
        ]
        # Deletion at 3 should sort before substitution at 15
        assert _allele_label(mutations) == "del:3-2bp;sub:15C>G"

    def test_indel_calls_split_indel(self):
        """INDEL type should pass through _split_indel."""
        mutations = [{"type": "indel", "ref_pos": 5, "length": 3,
                      "ref_base": "ACG", "query_base": "A--"}]
        result = _allele_label(mutations)
        assert "del:" in result


# ═══════════════════════════════════════════════════════════════
# _split_indel tests
# ═══════════════════════════════════════════════════════════════

class TestSplitIndel:
    def test_mixed_del_ins(self):
        """INDEL with both deletion and insertion parts should split."""
        m = {"ref_base": "ACG--", "query_base": "A--TT", "ref_pos": 5, "length": 4}
        parts = _split_indel(m)
        assert len(parts) >= 1
        assert any("del:" in p for p in parts)
        assert any("ins:" in p for p in parts)

    def test_indel_without_bases_falls_back(self):
        """INDEL without ref_base/query_base should use positional fallback."""
        m = {"ref_pos": 10, "length": 2}
        parts = _split_indel(m)
        assert parts == ["indel:10-2bp"]

    def test_deletion_only_indel(self):
        """INDEL containing only deletion should produce a 'del:' part.

        With ref='ACGTT', query='A----', the 'A' at pos=1 matches, then
        'CGTT' at pos=2-5 are deleted. So ref_pos advances past the match
        before the deletion starts."""
        m = {"ref_base": "ACGTT", "query_base": "A----", "ref_pos": 1, "length": 4}
        parts = _split_indel(m)
        del_parts = [p for p in parts if p.startswith("del:")]
        assert len(del_parts) >= 1
        # ref_base="ACGTT", query_base="A----":
        # 'A' at pos=1 matches; 'CGTT' at pos=2..5 deleted
        # ref_pos advances past the match before the deletion starts
        assert del_parts[0] == "del:2-4bp"


# ═══════════════════════════════════════════════════════════════
# _mutation_type_label tests
# ═══════════════════════════════════════════════════════════════

class TestMutationTypeLabel:
    def test_empty_returns_wt(self):
        assert _mutation_type_label([]) == "wt"

    def test_single_type_del(self):
        mutations = [{"type": "deletion"}]
        assert _mutation_type_label(mutations) == "del"

    def test_single_type_sub(self):
        mutations = [{"type": "substitution"}]
        assert _mutation_type_label(mutations) == "sub"

    def test_combined_types(self):
        mutations = [{"type": "deletion"}, {"type": "insertion"}]
        assert _mutation_type_label(mutations) == "del+ins"

    def test_all_types_sorted(self):
        mutations = [{"type": "deletion"}, {"type": "insertion"},
                     {"type": "indel"}, {"type": "substitution"}]
        assert _mutation_type_label(mutations) == "del+ins+indel+sub"


# ═══════════════════════════════════════════════════════════════
# _mutation_overlaps_target tests
# ═══════════════════════════════════════════════════════════════

class TestMutationOverlapsTarget:
    def test_deletion_overlaps(self):
        """Deletion spanning [10, 14] should overlap target [12, 18]."""
        m = {"type": "deletion", "ref_pos": 10, "length": 5}
        assert _mutation_overlaps_target(m, 12, 18) is True

    def test_deletion_before_target(self):
        """Deletion ending at 9 should NOT overlap target starting at 10."""
        m = {"type": "deletion", "ref_pos": 5, "length": 5}
        assert _mutation_overlaps_target(m, 10, 20) is False

    def test_substitution_inside_target(self):
        """Substitution at pos 15 should overlap target [10, 20]."""
        m = {"type": "substitution", "ref_pos": 15}
        assert _mutation_overlaps_target(m, 10, 20) is True

    def test_substitution_outside_target(self):
        """Substitution at pos 25 should NOT overlap target [10, 20]."""
        m = {"type": "substitution", "ref_pos": 25}
        assert _mutation_overlaps_target(m, 10, 20) is False

    def test_negative_ref_pos(self):
        """Mutation with negative ref_pos should return False."""
        m = {"type": "substitution", "ref_pos": -1}
        assert _mutation_overlaps_target(m, 0, 100) is False

    def test_indel_spans_target_edge(self):
        """INDEL starting before target but ending inside should overlap."""
        m = {"type": "indel", "ref_pos": 8, "length": 5}
        assert _mutation_overlaps_target(m, 10, 20) is True

    def test_deletion_exactly_at_target_start(self):
        """Deletion starting exactly at target_start should overlap."""
        m = {"type": "deletion", "ref_pos": 10, "length": 1}
        assert _mutation_overlaps_target(m, 10, 20) is True
