"""Validation tests — Verify pipeline accuracy using synthetic sequences with known ground truth

Existing tests only verify "the function runs without errors"; this file validates accuracy:

1. Mutation detection accuracy: whether mutations at known positions/types/lengths are correctly identified
2. Editing efficiency verification: whether editing efficiency computed from synthetic data with known efficiency is consistent
3. Position-aware gap penalty: whether gaps preferentially open in cutsite regions
"""

import sys
import os
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from crisviper import (
    Pipeline, PipelineConfig, QueryRecord,
    align_single,
    get_amplicon_structure, CutsiteRegion,
    affine_gap_alignment,
    affine_gap_alignment_position_aware,
    build_gradient_profiles,
    lineage_tracer_align,
    extract_mutations,
    MutationEvent, MutationType,
)

# ═══════════════════════════════════════════════════════════════
# Standard CARLIN reference sequence (332bp)
# ═══════════════════════════════════════════════════════════════
CARLIN_REF = (
    "TATGTGTGGGAGGGCTAAGAGG"   # Primer5 (23bp)
    "CCGCC"                     # Prefix (5bp)
    "GACTGCACGACAGTCGA"         # Target1: conserved region 13bp + cutsite 7bp
    "CGATGGAG"                  # Linker (7bp)
    "TCGACACGACTCGCGCA"         # Target2
    "TACGATGG"                  # Linker
    "AGTCGACTACAGTCGCTA"        # Target3
    "CGACGATG"                  # Linker
    "GAGTCGCGAGCGCTATG"         # Target4
    "AGCGACTA"                  # Linker
    "TGGAGTCGATACGATACG"        # Target5
    "CGCACGCT"                  # Linker
    "ATGGAGTCGAGAGCGCGC"        # Target6
    "TCGTCAAC"                  # Linker
    "GATGGAGTCGCGACTGTA"        # Target7
    "CGCACTCG"                  # Linker
    "CGATGGAGTCGATAGTAT"        # Target8
    "GCGTACAC"                  # Linker
    "GCGATGGAGTCGACTGCA"        # Target9
    "CGACAGTC"                  # Linker
    "GACTATGGAGTCGATACGTAGC"    # Target10
    "ACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"  # Postfix(8) + Primer3(33)
)

# Cutsite coordinates on the full-length reference sequence (0-indexed, inclusive)
CUTSITES = [
    CutsiteRegion("Target1",  41,  47),
    CutsiteRegion("Target2",  68,  74),
    CutsiteRegion("Target3",  95, 101),
    CutsiteRegion("Target4", 122, 128),
    CutsiteRegion("Target5", 149, 155),
    CutsiteRegion("Target6", 176, 182),
    CutsiteRegion("Target7", 203, 209),
    CutsiteRegion("Target8", 230, 236),
    CutsiteRegion("Target9", 257, 263),
    CutsiteRegion("Target10", 284, 290),
]


# ═══════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════

def _make_deletion(seq: str, pos: int, length: int) -> str:
    """Delete `length` bases starting at position `pos`"""
    return seq[:pos] + seq[pos+length:]


def _make_insertion(seq: str, pos: int, insert: str) -> str:
    """Insert `insert` at position `pos`"""
    return seq[:pos] + insert + seq[pos:]


def _make_substitution(seq: str, pos: int, new_base: str) -> str:
    """Replace the base at position `pos` with `new_base`"""
    lst = list(seq)
    lst[pos] = new_base
    return ''.join(lst)


def _run_align_single(
    query_seq: str,
    read_name: str = "test",
    read_count: int = 50,
    lineage_mode: bool = True,
    min_reads_sub: int = 0,
    min_reads_indel: int = 0,
) -> 'AlignmentResult':
    """Align a single sequence via align_single, returns AlignmentResult"""
    config = PipelineConfig(
        lineage_mode=lineage_mode,
        primer5_len=23,
        primer3_len=33,
        primer5_threshold=19,
        primer3_threshold=29,
        min_reads_sub=min_reads_sub,
        min_reads_indel=min_reads_indel,
    )
    query = QueryRecord(
        readName=read_name, cellBC="test", UMI="UMI",
        readCount=read_count, seq=query_seq,
    )
    cutsites = get_amplicon_structure(CARLIN_REF) if lineage_mode else None
    return align_single(query, CARLIN_REF, config, cutsites)


def _run_pipeline(
    queries,
    lineage_mode: bool = True,
    min_reads_sub: int = 0,
    min_reads_indel: int = 0,
):
    """Run the full pipeline, returns PipelineResult"""
    config = PipelineConfig(
        lineage_mode=lineage_mode,
        primer5_len=23,
        primer3_len=33,
        primer5_threshold=19,
        primer3_threshold=29,
        min_reads_sub=min_reads_sub,
        min_reads_indel=min_reads_indel,
        threads=1,
    )
    pipeline = Pipeline(config=config, ref_seq=CARLIN_REF)
    if lineage_mode:
        pipeline.load_cutsites()
    return pipeline.run(queries)


# ═══════════════════════════════════════════════════════════════
# Test 1: Mutation detection accuracy
# ═══════════════════════════════════════════════════════════════

class TestMutationDetectionAccuracy:
    """Verify mutation identification accuracy using known ground truth"""

    def test_wildtype_no_mutations(self):
        """Wildtype sequence → no mutations"""
        result = _run_align_single(CARLIN_REF, "wt")
        assert result.success, f"Wildtype alignment failed: {result.error}"
        assert len(result.mutations) == 0, \
            f"Wildtype should have no mutations but detected {len(result.mutations)}: {[(m.type, m.ref_pos) for m in result.mutations]}"
        assert result.stats.mismatches == 0
        assert not result.stats.has_indel

    def test_3bp_deletion_at_cutsite(self):
        """3bp deletion inside Target1 cutsite → identified as DELETION, length=3, inside window"""
        query = _make_deletion(CARLIN_REF, 42, 3)
        result = _run_align_single(query, "del_3bp_cutsite")
        assert result.success, f"3bp deletion alignment failed: {result.error}"

        deletions = [m for m in result.mutations if m.type == MutationType.DELETION]
        assert len(deletions) >= 1, \
            f"No deletion detected, found: {[(m.type, m.length) for m in result.mutations]}"

        del_event = deletions[0]
        assert del_event.length == 3, f"Deletion length should be 3, got {del_event.length}"
        assert del_event.in_cutsite_window, \
            f"Deletion (ref_pos={del_event.ref_pos}) should be inside cutsite window"

    def test_5bp_insertion_at_cutsite(self):
        """5bp insertion inside Target1 cutsite → identified as INSERTION, length=5, inside window"""
        query = _make_insertion(CARLIN_REF, 42, "ACGTA")
        result = _run_align_single(query, "ins_5bp_cutsite")
        assert result.success, f"5bp insertion alignment failed: {result.error}"

        insertions = [m for m in result.mutations if m.type == MutationType.INSERTION]
        assert len(insertions) >= 1, \
            f"No insertion detected, found: {[(m.type, m.length) for m in result.mutations]}"

        ins_event = insertions[0]
        assert ins_event.length == 5, f"Insertion length should be 5, got {ins_event.length}"
        assert ins_event.in_cutsite_window, \
            f"Insertion (ref_pos={ins_event.ref_pos}) should be inside cutsite window"

    def test_point_mutation_at_cutsite(self):
        """Point mutation inside Target1 cutsite → identified as SUBSTITUTION, length=1, inside window"""
        query = _make_substitution(CARLIN_REF, 44, "T")
        result = _run_align_single(query, "sub_cutsite")
        assert result.success, f"Point mutation alignment failed: {result.error}"

        subs = [m for m in result.mutations if m.type == MutationType.SUBSTITUTION]
        assert len(subs) >= 1, \
            f"No substitution detected, found: {[(m.type, m.length) for m in result.mutations]}"

        sub_event = subs[0]
        assert sub_event.length == 1, f"Substitution length should be 1, got {sub_event.length}"
        assert sub_event.in_cutsite_window, \
            f"Substitution (ref_pos={sub_event.ref_pos}) should be inside cutsite window"

    def test_7bp_deletion_full_cutsite(self):
        """Complete deletion of a cutsite (7bp) → identified as DELETION, length=7"""
        query = _make_deletion(CARLIN_REF, 41, 7)  # Delete entire Target1 cutsite
        result = _run_align_single(query, "del_7bp_full_cutsite")
        assert result.success, f"7bp deletion alignment failed: {result.error}"

        deletions = [m for m in result.mutations if m.type == MutationType.DELETION]
        if not deletions:
            deletions = [m for m in result.mutations if m.type == MutationType.INDEL]
        assert len(deletions) >= 1, \
            f"No deletion detected: {[(m.type, m.length) for m in result.mutations]}"

        # Total deletion length (may be distributed across multiple events)
        total_del_len = sum(m.length for m in result.mutations
                            if m.type in (MutationType.DELETION, MutationType.INDEL))
        assert total_del_len >= 7, \
            f"Total deletion length should be >= 7, got {total_del_len}"

    def test_mutations_in_multiple_cutsites(self):
        """Mutations on multiple cutsites → all correctly identified"""
        # Target1 cutsite 3bp deletion + Target2 cutsite point mutation
        query = _make_deletion(CARLIN_REF, 42, 3)
        # Note: coordinates shift after deletion, Target2 goes from 68 to 68-3=65
        # CARLIN_REF[68]=='A', use 'T' to create a truly different base
        query = _make_substitution(query, 65, "T")  # Original Target2 position
        result = _run_align_single(query, "multi_cutsites")
        assert result.success, f"Multi-cutsite mutation alignment failed: {result.error}"

        types = {m.type for m in result.mutations}
        has_deletion = MutationType.DELETION in types or MutationType.INDEL in types
        has_substitution = MutationType.SUBSTITUTION in types
        assert has_deletion, f"No deletion detected, types found: {types}"
        assert has_substitution, f"No substitution detected, types found: {types}"

    def test_large_deletion_across_target(self):
        """Large deletion spanning a target → correctly identified"""
        query = _make_deletion(CARLIN_REF, 41, 20)  # Delete most of Target1 region
        result = _run_align_single(query, "del_large_20bp")
        assert result.success, f"Large deletion alignment failed: {result.error}"

        # Should have a deletion event
        deletions = [m for m in result.mutations if m.type in
                     (MutationType.DELETION, MutationType.INDEL)]
        total_del_len = sum(m.length for m in deletions)
        assert total_del_len >= 15, \
            f"Large deletion length should be >= 15, detected {total_del_len}"

    def test_standard_mode_also_detects_mutations(self):
        """Standard mode (non-lineage) also correctly detects known mutations"""
        query = _make_deletion(CARLIN_REF, 42, 3)
        result = _run_align_single(query, "del_std", lineage_mode=False)
        assert result.success, f"Standard mode 3bp deletion alignment failed: {result.error}"

        deletions = [m for m in result.mutations if m.type == MutationType.DELETION]
        assert len(deletions) >= 1, \
            f"Standard mode did not detect deletion: {[(m.type, m.length) for m in result.mutations]}"
        assert deletions[0].length == 3


# ═══════════════════════════════════════════════════════════════
# Test 2: Editing efficiency verification
# ═══════════════════════════════════════════════════════════════

class TestEditingEfficiency:
    """Verify editing efficiency calculation using synthetic data with known efficiency"""

    def test_known_editing_efficiency(self):
        """Known editing efficiency = 60% (12/20 sequences with mutations)"""
        queries = []

        # 12 mutated sequences (different types)
        for i in range(6):
            queries.append(QueryRecord(
                readName=f"del_{i}", cellBC="test", UMI=f"D{i}",
                readCount=10,
                seq=_make_deletion(CARLIN_REF, 42, 3),
            ))
        for i in range(6):
            queries.append(QueryRecord(
                readName=f"ins_{i}", cellBC="test", UMI=f"I{i}",
                readCount=10,
                seq=_make_insertion(CARLIN_REF, 42, "ACG"),
            ))

        # 8 wildtype sequences
        for i in range(8):
            queries.append(QueryRecord(
                readName=f"wt_{i}", cellBC="test", UMI=f"W{i}",
                readCount=10, seq=CARLIN_REF,
            ))

        result = _run_pipeline(queries)

        assert result.stats.total_queries == 20
        assert result.stats.successful == 20, \
            f"Expected 20 successful, got {result.stats.successful}, failed {result.stats.failed}"
        assert result.stats.mutated_sequences == 12, \
            f"Expected 12 mutated sequences, detected {result.stats.mutated_sequences}"

        expected_eff = 12.0 / 20.0 * 100.0
        actual_eff = result.stats.editing_efficiency_pct
        assert actual_eff == pytest.approx(expected_eff, abs=1.0), \
            f"Editing efficiency: expected {expected_eff:.1f}%, actual {actual_eff:.1f}%"

    def test_zero_editing_efficiency(self):
        """All wildtype → editing efficiency = 0%"""
        queries = [
            QueryRecord(readName=f"wt_{i}", cellBC="test", UMI=f"W{i}",
                        readCount=10, seq=CARLIN_REF)
            for i in range(10)
        ]
        result = _run_pipeline(queries)
        assert result.stats.editing_efficiency_pct == pytest.approx(0.0, abs=0.1)

    def test_full_editing_efficiency(self):
        """All mutated → editing efficiency = 100%"""
        queries = [
            QueryRecord(readName=f"mut_{i}", cellBC="test", UMI=f"M{i}",
                        readCount=10, seq=_make_deletion(CARLIN_REF, 42, 3))
            for i in range(10)
        ]
        result = _run_pipeline(queries)
        assert result.stats.editing_efficiency_pct == pytest.approx(100.0, abs=0.1)

    def test_editing_efficiency_with_low_readcount_filtering(self):
        """Low readCount sequences should be filtered without affecting editing efficiency"""
        queries = []
        # High readCount mutated sequences (should be kept)
        for i in range(5):
            queries.append(QueryRecord(
                readName=f"good_del_{i}", cellBC="test", UMI=f"G{i}",
                readCount=50,
                seq=_make_deletion(CARLIN_REF, 42, 3),
            ))
        # Low readCount mutated sequences (should be filtered)
        for i in range(5):
            queries.append(QueryRecord(
                readName=f"bad_sub_{i}", cellBC="test", UMI=f"B{i}",
                readCount=2,
                seq=_make_substitution(CARLIN_REF, 44, "T"),
            ))
        # High readCount wildtype sequences
        for i in range(5):
            queries.append(QueryRecord(
                readName=f"good_wt_{i}", cellBC="test", UMI=f"H{i}",
                readCount=50, seq=CARLIN_REF,
            ))

        config = PipelineConfig(
            lineage_mode=True,
            primer5_len=23, primer3_len=33,
            primer5_threshold=19, primer3_threshold=29,
            min_reads_sub=10,
            min_reads_indel=3,
            threads=1,
        )
        pipeline = Pipeline(config=config, ref_seq=CARLIN_REF)
        pipeline.load_cutsites()
        result = pipeline.run(queries)

        # 5 good_del + 5 wt = 10 successful (5 bad_sub filtered as false positives)
        assert result.stats.successful == 10, \
            f"Expected 10 successful, got {result.stats.successful} (failed: {result.stats.failed})"
        # 5 good_del (bad_sub filtered as false positives, only deletions kept)
        assert result.stats.mutated_sequences == 5
        # Editing efficiency = 5/10 * 100 = 50%
        assert result.stats.editing_efficiency_pct == pytest.approx(50.0, abs=0.1)

    def test_stats_consistency(self):
        """Statistical data consistency check"""
        queries = []
        for i in range(8):
            queries.append(QueryRecord(
                readName=f"mut_{i}", cellBC="test", UMI=f"M{i}",
                readCount=10, seq=_make_deletion(CARLIN_REF, 42, 3),
            ))
        for i in range(5):
            queries.append(QueryRecord(
                readName=f"wt_{i}", cellBC="test", UMI=f"W{i}",
                readCount=10, seq=CARLIN_REF,
            ))

        result = _run_pipeline(queries)
        stats = result.stats

        # Consistency: successful = mutated + unmutated
        assert stats.mutated_sequences + stats.unmutated_sequences == stats.successful
        # Consistency: total_queries = successful + failed
        assert stats.successful + stats.failed == stats.total_queries
        # Editing efficiency bounds
        assert 0.0 <= stats.editing_efficiency_pct <= 100.0


# ═══════════════════════════════════════════════════════════════
# Test 3: Position-aware gap penalty
# ═══════════════════════════════════════════════════════════════

class TestPositionAwareGapPenalty:
    """Verify that gaps preferentially open in cutsite regions rather than conserved regions"""

    def test_build_gradient_profiles_correctness(self):
        """Numerical correctness of gradient penalty profiles"""
        ref_len = 20
        cutsites = [CutsiteRegion("T1", 5, 10)]
        go, ge, mp, _ = build_gradient_profiles(
            ref_len, cutsites,
            base_gap_open=-2.0, base_gap_extend=-0.1,
            mismatch_penalty=-3.0,
            min_scale=1.0, max_scale=2.0, cutsite_edge_scale=2.0,
        )

        # far region (0-1, 14-19): max_scale=2.0 → open=-4.0, extend=-0.2
        for i in [0, 1, 15, 16, 17, 18, 19]:
            assert go[i] == pytest.approx(-4.0), f"far go[{i}]"
            assert ge[i] == pytest.approx(-0.2), f"far ge[{i}]"

        # cutsite center (7-8): close to min_scale=1.0 → open ≈ -2.0
        assert go[7] > -2.5, f"center go[7]={go[7]} should be near -2.0"
        assert go[8] > -2.5, f"center go[8]={go[8]} should be near -2.0"
        # cutsite edge (5, 10): close to edge_scale=2.0 → open ≈ -4.0
        assert go[5] < -3.0, f"edge go[5]={go[5]} should be near -4.0"
        assert go[10] < -3.0, f"edge go[10]={go[10]} should be near -4.0"

    def test_position_aware_alignment_gap_in_cutsite(self):
        """Position-aware alignment places gaps in cutsite regions (low-cost regions)"""
        # Design a sequence with conserved ends and a cutsite in the middle
        # Conserved regions (0-3, 8-11) + cutsite(4-7)
        ref = "AAATTTCCCGGG"    # 12bp
        query = "AAAGGG"         # 6bp, missing the middle TTTCCC (6bp)
        cutsites = [CutsiteRegion("T1", 4, 7)]

        # Position-aware alignment (low cutsite gap penalty)
        go, ge, mp, _ = build_gradient_profiles(
            len(ref), cutsites,
            base_gap_open=-2.0, base_gap_extend=-0.1,
            mismatch_penalty=-3.0,
            min_scale=1.0, max_scale=2.0, cutsite_edge_scale=2.0,
        )
        _, ar, aq, stats = affine_gap_alignment_position_aware(
            ref, query, go, ge,
            mismatch_penalty_profile=mp,
        )

        # Should have a gap (6bp deletion)
        assert stats["gaps_in_query"] >= 6, \
            f"Expected gap >= 6, got {stats['gaps_in_query']}"

        # Gap should be located in the cutsite region (4-7)
        # Ideal case: AAA-----GGG → gap from position 3 to 8
        gap_indices = [i for i, c in enumerate(aq) if c == '-']
        assert len(gap_indices) >= 6, f"Gap column count: {len(gap_indices)}"

        # Gap should include the cutsite region (4-7)
        gap_set = set(gap_indices)
        cutsite_set = set(range(4, 8))
        overlap = gap_set & cutsite_set
        assert len(overlap) >= 3, \
            f"Insufficient gap-cutsite overlap: gap={gap_indices}, cutsite=4-7, overlap={overlap}"

    def test_position_aware_vs_standard_gap_placement(self):
        """Gap placement difference between position-aware and standard alignment: position-aware should prefer cutsite"""
        ref = "AAATTTCCCGGGAAA"  # 15bp
        query = "AAACCCAAA"       # 9bp, missing TTTGGG (6bp)
        cutsites = [CutsiteRegion("T1", 3, 5)]  # TTT region

        # Standard alignment (uniform penalty)
        _, ar_std, aq_std, _ = affine_gap_alignment(
            ref, query,
        )

        # Position-aware alignment (low cutsite gap penalty, higher flanking penalty)
        go, ge, mp, _ = build_gradient_profiles(
            len(ref), cutsites,
            base_gap_open=-2.0, base_gap_extend=-0.1,
            mismatch_penalty=-3.0,
            min_scale=1.0, max_scale=5.0, cutsite_edge_scale=5.0,
        )
        _, ar_pa, aq_pa, stats_pa = affine_gap_alignment_position_aware(
            ref, query, go, ge,
            mismatch_penalty_profile=mp,
        )

        # Position-aware alignment should have a gap
        assert stats_pa["gaps_in_query"] >= 6

        # Gaps should be mostly in the cutsite region (3-5)
        pa_gaps = {i for i, c in enumerate(aq_pa) if c == '-'}
        pa_in_cutsite = sum(1 for i in pa_gaps if 3 <= i <= 5)
        pa_outside = len(pa_gaps) - pa_in_cutsite

        # Most gaps should be in cutsite
        assert pa_in_cutsite >= pa_outside, \
            f"Position-aware: gaps in cutsite={pa_in_cutsite}, outside={pa_outside}"

    def test_lineage_tracer_align_detects_cutsite_indel(self):
        """lineage_tracer_align correctly detects indels within the cutsite region"""
        ref = "A" * 15 + "GAGTCG" + "A" * 15  # GAGTCG is the cutsite motif
        query = "A" * 15 + "A" * 15              # cutsite completely deleted
        cutsites = [CutsiteRegion("T1", 15, 20)]  # GAGTCG position

        score, ar, aq, stats = lineage_tracer_align(
            ref, query, cutsites,
            min_scale=1.0,
            max_scale=2.0,
        )

        # Should detect a 6bp deletion
        assert stats["gaps_in_query"] >= 5, \
            f"lineage_tracer_align did not detect a large enough deletion: gaps_in_query={stats['gaps_in_query']}"

        # Gap should be near the cutsite
        gap_start = -1
        gap_end = -1
        for i, c in enumerate(aq):
            if c == '-':
                if gap_start == -1:
                    gap_start = i
                gap_end = i

        assert gap_start >= 14 or gap_start == -1, \
            f"Gap start position {gap_start} should be near cutsite(15-20)"


# ═══════════════════════════════════════════════════════════════
# Test 4: extract_mutations direct verification
# ═══════════════════════════════════════════════════════════════

class TestExtractMutationsAccuracy:
    """Directly verify extract_mutations function accuracy on known mutations"""

    def test_extract_single_substitution(self):
        """Known point mutation position and base are correct"""
        ar = "ACGTACGT"
        aq = "ACCTACGT"  # G→C at ref pos 2
        mutations = extract_mutations(ar, aq)

        assert len(mutations) == 1
        m = mutations[0]
        assert m.type == MutationType.SUBSTITUTION
        assert m.ref_pos == 2
        assert m.ref_base == "G"
        assert m.query_base == "C"
        assert m.length == 1

    def test_extract_deletion_length(self):
        """Known deletion length is correct"""
        ar = "ACGTACGT"
        aq = "A-----GT"  # CGTAC deleted (5bp)
        mutations = extract_mutations(ar, aq)

        deletions = [m for m in mutations if m.type == MutationType.DELETION]
        assert len(deletions) == 1
        assert deletions[0].length == 5

    def test_extract_deletion_with_cutsite_window(self):
        """Deletion inside cutsite window, in_cutsite_window=True"""
        # 3bp deletion at ref positions 3-5 (TAC removed from ACGTACGTACGT)
        ar = "ACGTACGTACGT"
        aq = "ACG---GTACGT"  # dashes at positions 3-5 for TAC deletion
        cutsites = [CutsiteRegion("T1", 3, 5)]
        mutations = extract_mutations(ar, aq, cutsites=cutsites, sub_window=3)

        deletions = [m for m in mutations if m.type == MutationType.DELETION]
        assert len(deletions) >= 1
        assert deletions[0].in_cutsite_window, "Deletion inside cutsite window should be marked in_window"

    def test_extract_deletion_outside_cutsite_window(self):
        """Deletion outside cutsite window, in_cutsite_window=False"""
        # 3bp deletion at ref positions 3-5, cutsite at 8-10 (far away)
        ar = "ACGTACGTACGT"
        aq = "ACG---GTACGT"  # dashes at positions 3-5 for TAC deletion
        cutsites = [CutsiteRegion("T1", 8, 10)]  # cutsite far from deletion
        mutations = extract_mutations(ar, aq, cutsites=cutsites, sub_window=1)

        deletions = [m for m in mutations if m.type == MutationType.DELETION]
        assert len(deletions) >= 1
        assert not deletions[0].in_cutsite_window, \
            "Deletion outside cutsite window should be marked not_in_window (false)"

    def test_extract_insertion_length(self):
        """Known insertion length is correct"""
        ar = "AC--GT"
        aq = "ACGGGT"  # GG inserted (2bp)
        mutations = extract_mutations(ar, aq)

        insertions = [m for m in mutations if m.type == MutationType.INSERTION]
        assert len(insertions) == 1
        assert insertions[0].length == 2

    def test_no_mutations_for_identical(self):
        """Identical sequences have no mutations"""
        assert extract_mutations("ACGT", "ACGT") == []

    def test_adjacent_indel_merged_to_complex(self):
        """Adjacent insertion + deletion merged into a complex event"""
        ar = "ACGT--AC"
        aq = "AC--GGAC"  # DEL at pos 3-4 (GT), INS at pos 5-6 (GG)
        mutations = extract_mutations(ar, aq)

        types = {m.type for m in mutations}
        assert MutationType.INDEL in types, \
            f"Adjacent insertion + deletion should merge to INDEL: {types}"

    def test_empty_inputs(self):
        """Empty input returns empty list"""
        assert extract_mutations("", "") == []

    def test_mismatched_lengths(self):
        """Mismatched lengths return empty list"""
        assert extract_mutations("ACGT", "AC") == []


# ═══════════════════════════════════════════════════════════════
# Standalone test functions (for import by run_tests.py)
# ═══════════════════════════════════════════════════════════════

def run_accuracy_checks(check_func):
    """Run core accuracy assertions using check_func(name, condition)

    This function is designed for run_tests.py, using a consistent check API.
    The check_func parameter should be a callable of (str, bool) -> None.
    """
    print("\n── Accuracy Validation Tests ──\n")

    # ── 1a: Mutation detection ──
    result = _run_align_single(CARLIN_REF, "wt", read_count=50)
    check_func("Wildtype has no mutations", result.success and len(result.mutations) == 0)

    result = _run_align_single(_make_deletion(CARLIN_REF, 42, 3), "del_3bp")
    deletions = [m for m in result.mutations if m.type == MutationType.DELETION]
    check_func("3bp deletion detection", result.success
               and len(deletions) >= 1
               and deletions[0].length == 3
               and deletions[0].in_cutsite_window)

    result = _run_align_single(_make_insertion(CARLIN_REF, 42, "ACGTA"), "ins_5bp")
    insertions = [m for m in result.mutations if m.type == MutationType.INSERTION]
    check_func("5bp insertion detection", result.success
               and len(insertions) >= 1
               and insertions[0].length == 5
               and insertions[0].in_cutsite_window)

    result = _run_align_single(_make_substitution(CARLIN_REF, 44, "T"), "sub_1bp")
    subs = [m for m in result.mutations if m.type == MutationType.SUBSTITUTION]
    check_func("Point mutation detection", result.success
               and len(subs) >= 1
               and subs[0].length == 1
               and subs[0].in_cutsite_window)

    # ── Editing efficiency verification ──
    queries = []
    for i in range(6):
        queries.append(QueryRecord(
            readName=f"mut_{i}", cellBC="test", UMI=f"M{i}",
            readCount=10,
            seq=_make_deletion(CARLIN_REF, 42, 3),
        ))
    for i in range(4):
        queries.append(QueryRecord(
            readName=f"wt_{i}", cellBC="test", UMI=f"W{i}",
            readCount=10, seq=CARLIN_REF,
        ))
    result = _run_pipeline(queries)
    expected_eff = 60.0
    actual_eff = result.stats.editing_efficiency_pct
    check_func(f"Editing efficiency {expected_eff}% (expected={expected_eff}, actual={actual_eff:.1f})",
               abs(actual_eff - expected_eff) < 1.0)

    # ── Position-aware gap ──
    ref = "AAATTTCCCGGG"
    query = "AAAGGG"
    cutsites = [CutsiteRegion("T1", 4, 7)]
    go, ge, mp = build_gradient_profiles(
        len(ref), cutsites,
        base_gap_open=-2.0, base_gap_extend=-0.1,
        mismatch_penalty=-3.0,
        min_scale=1.0, max_scale=2.0, cutsite_edge_scale=2.0,
    )
    _, ar, aq, stats = affine_gap_alignment_position_aware(
        ref, query, go, ge,
        mismatch_penalty_profile=mp,
    )
    gap_indices = {i for i, c in enumerate(aq) if c == '-'}
    cutsite_range = set(range(4, 8))
    overlap = len(gap_indices & cutsite_range)
    check_func("Position-aware gap in cutsite region",
               stats["gaps_in_query"] >= 6 and overlap >= 3)

    # ── extract_mutations direct verification ──
    m = extract_mutations("ACGT", "ACCT")[0]
    check_func("extract_mutations point mutation ref_pos",
               m.type == MutationType.SUBSTITUTION and m.ref_pos == 2)
    m = extract_mutations("ACGTACGT", "A-----GT")[0]
    check_func("extract_mutations 5bp deletion",
               m.type == MutationType.DELETION and m.length == 5)
    check_func("extract_mutations identical sequences no mutation",
               extract_mutations("ACGT", "ACGT") == [])

    print()  # blank line separator
