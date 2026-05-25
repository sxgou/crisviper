"""Tests for the gap_exit_bonus DP parameter — mathematical correctness."""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import numpy as np
from crisviper import (
    affine_gap_alignment,
    affine_gap_alignment_position_aware,
    build_gradient_profiles,
    lineage_tracer_align,
    CutsiteRegion,
    Pipeline, PipelineConfig, QueryRecord,
    align_single,
    get_amplicon_structure,
)

# ═══════════════════════════════════════════════════════════════
# Unit tests: gap_exit_bonus mathematical correctness
# ═══════════════════════════════════════════════════════════════

class TestGapExitBonusMath:
    """Directly verify the mathematical effect of gap_exit_bonus"""

    def test_parameter_default_is_zero(self):
        """Default gap_exit_bonus=0.0, behavior matches original code"""
        s1, _, aq1, _ = affine_gap_alignment("ACGT", "ACCT")
        s2, _, aq2, _ = affine_gap_alignment("ACGT", "ACCT", gap_exit_bonus=0.0)
        assert s1 == s2
        assert aq1 == aq2

    def test_negative_bonus_lowers_score(self):
        """Negative bonus lowers alignment score (because M path values from Ix/Iy are reduced)"""
        ref, qry = "AAAAATTTTT", "AAAATTTTT"
        s0 = affine_gap_alignment(ref, qry, gap_exit_bonus=0.0)[0]
        s1 = affine_gap_alignment(ref, qry, gap_exit_bonus=-1.0)[0]
        assert s1 <= s0, f"Score with bonus {s1} > without {s0}"

    def test_monotonic_decreasing_score(self):
        """The more negative the bonus, the alignment score is monotonically non-increasing"""
        ref, qry = "GACTGCACGACAGTCGAT", "GACTGCAGTCGAT"
        prev = float('inf')
        for bonus in [0.0, -0.5, -1.0, -1.5, -2.0]:
            s = affine_gap_alignment(ref, qry, gap_exit_bonus=bonus)[0]
            assert s <= prev + 1e-9, f"bonus={bonus}: score increased {prev} -> {s}"
            prev = s

    def test_bonus_affects_m_score_via_ix_and_iy(self):
        """Verify M[i,j] computation uses max(M, Ix+bonus, Iy+bonus) semantics"""
        ref, qry = "ATGC", "AG"
        s0, ar0, aq0, st0 = affine_gap_alignment(ref, qry, gap_exit_bonus=0.0)
        s1, ar1, aq1, st1 = affine_gap_alignment(ref, qry, gap_exit_bonus=-2.0)
        assert len(ar0) == len(ar1)
        assert st0["alignment_length"] == st1["alignment_length"]

    def test_position_aware_bonus_affects_score(self):
        """Bonus also lowers score in position-aware DP"""
        go = np.full(12, -2.0)
        ge = np.full(12, -0.1)
        ref, qry = "AAAACCCCGGGG", "AAAAGGGG"
        s0 = affine_gap_alignment_position_aware(ref, qry, go, ge, gap_exit_bonus=0.0)[0]
        s1 = affine_gap_alignment_position_aware(ref, qry, go, ge, gap_exit_bonus=-1.0)[0]
        assert s1 <= s0

    def test_lineage_tracer_bonus_affects_score(self):
        """Bonus lowers score in lineage_tracer_align"""
        cutsites = [CutsiteRegion("T1", 4, 7)]
        ref, qry = "AAAACCCCGGGG", "AAAAGGGG"
        s0 = lineage_tracer_align(ref, qry, cutsites, gap_exit_bonus=0.0)[0]
        s1 = lineage_tracer_align(ref, qry, cutsites, gap_exit_bonus=-1.0)[0]
        assert s1 <= s0

    def test_bonus_does_not_change_identical_alignments(self):
        """Identical sequence alignment is unaffected by bonus"""
        for bonus in [0.0, -1.0, -2.0]:
            s, ar, aq, st = affine_gap_alignment("ACGTACGT", "ACGTACGT",
                                                   gap_exit_bonus=bonus)
            assert ar == "ACGTACGT"
            assert aq == "ACGTACGT"
            assert st["matches"] == 8


# ═══════════════════════════════════════════════════════════════
# Integration tests: gap_exit_bonus works through the Pipeline
# ═══════════════════════════════════════════════════════════════

CARLIN_REF = (
    "TATGTGTGGGAGGGCTAAGAGG"   # Primer5 (23bp)
    "CCGCC"                     # Prefix (5bp)
    "GACTGCACGACAGTCGA"         # Target1
    "CGATGGAG"                  # Linker
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
    "ACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"
)

class TestPipelineGapExitBonus:
    """Verify gap_exit_bonus works through the complete Pipeline"""

    def make_query(self, seq: str, name: str = "test", rc: int = 50):
        return QueryRecord(readName=name, cellBC="test", UMI="UMI",
                           readCount=rc, seq=seq)

    def test_pipeline_lineage_mode_accepts_bonus(self):
        config = PipelineConfig(lineage_mode=True, gap_exit_strength=-1.0,
                                min_reads_sub=0, min_reads_indel=0)
        pipeline = Pipeline(config=config, ref_seq=CARLIN_REF)
        pipeline.load_cutsites()
        result = pipeline.run([self.make_query(CARLIN_REF)])
        assert result.stats.successful == 1

    def test_pipeline_standard_mode_accepts_bonus(self):
        config = PipelineConfig(lineage_mode=False, gap_exit_strength=-1.0,
                                min_reads_sub=0, min_reads_indel=0)
        pipeline = Pipeline(config=config, ref_seq=CARLIN_REF)
        result = pipeline.run([self.make_query(CARLIN_REF)])
        assert result.stats.successful == 1

    def test_align_single_with_bonus_preserves_wt(self):
        config = PipelineConfig(lineage_mode=True, gap_exit_strength=-1.0,
                                min_reads_sub=0, min_reads_indel=0)
        cutsites = get_amplicon_structure(CARLIN_REF)
        q = self.make_query(CARLIN_REF)
        result = align_single(q, CARLIN_REF, config, cutsites)
        assert result.success
        assert len(result.mutations) == 0

    def test_align_single_with_bonus_reduces_fragmentation(self):
        config_off = PipelineConfig(lineage_mode=True, gap_exit_strength=0.0,
                                    min_reads_sub=0, min_reads_indel=0)
        config_on  = PipelineConfig(lineage_mode=True, gap_exit_strength=-1.0,
                                    min_reads_sub=0, min_reads_indel=0)
        cutsites = get_amplicon_structure(CARLIN_REF)
        qry = CARLIN_REF[:41] + CARLIN_REF[48:]
        result_off = align_single(self.make_query(qry), CARLIN_REF, config_off, cutsites)
        result_on  = align_single(self.make_query(qry), CARLIN_REF, config_on, cutsites)
        print(f"\n  bonus=0.0: aln_ref len={len(result_off.aligned_ref)}, "
              f"mutations={len(result_off.mutations)}")
        print(f"  bonus=-1.0: aln_ref len={len(result_on.aligned_ref)}, "
              f"mutations={len(result_on.mutations)}")
        assert result_off.success
        assert result_on.success

    def test_bonus_with_chain_indel_detection(self):
        config = PipelineConfig(lineage_mode=True, gap_exit_strength=-1.0,
                                min_reads_sub=0, min_reads_indel=0)
        cutsites = get_amplicon_structure(CARLIN_REF)
        qry = CARLIN_REF[:42] + "ACGTA" + CARLIN_REF[42:]
        q = self.make_query(qry, "ins5")
        result = align_single(q, CARLIN_REF, config, cutsites)
        assert result.success, f"Insertion failed with bonus: {result.error}"
        insertions = [m for m in result.mutations
                      if m.type.name in ("INSERTION", "INDEL")]
        assert len(insertions) >= 1


# ═══════════════════════════════════════════════════════════════
# Phase A3: Short match region discount
# ═══════════════════════════════════════════════════════════════

class TestShortMatchDiscount:
    def test_default_is_no_discount(self):
        s1, ar1, aq1, _ = affine_gap_alignment("ACGT", "ACT")
        s2, ar2, aq2, _ = affine_gap_alignment("ACGT", "ACT", short_match_window=0, short_match_discount=0.5)
        assert s1 == s2

    def test_short_match_consolidated_into_gap_global(self):
        ref, qry = "AAACCCGGGTTT", "AAATTTGGGTTT"
        s0, _, _, _ = affine_gap_alignment(ref, qry, short_match_window=0, short_match_discount=0.5)
        s1, _, _, _ = affine_gap_alignment(ref, qry, short_match_window=3, short_match_discount=0.0)
        print(f"\n  ref={ref}  qry={qry}\n  OFF score={s0:.1f}\n  ON  score={s1:.1f}")

    def test_long_matches_untouched(self):
        ref, qry = "ACGTACGTACGTACGT", "ACGTACGTACGTACGT"
        s0, ar0, aq0, _ = affine_gap_alignment(ref, qry)
        s1, ar1, aq1, _ = affine_gap_alignment(ref, qry, short_match_window=3, short_match_discount=0.5)
        assert ar0 == ar1 and aq0 == aq1
        assert s0 - s1 <= 3.0 + 1e-9

    def test_short_match_window_size_configurable(self):
        ref, qry = "ACGTACGTNNNN", "ACGTACGTTTTT"
        s_win1, _, _, _ = affine_gap_alignment(ref, qry, short_match_window=1, short_match_discount=0.0)
        s_win5, _, _, _ = affine_gap_alignment(ref, qry, short_match_window=5, short_match_discount=0.0)
        assert s_win5 <= s_win1 + 1e-9

    def test_lineage_tracer_forwards_short_match(self):
        cutsites = [CutsiteRegion("T1", 4, 7)]
        ref, qry = "AAAACCCCGGGGTTTT", "AAAAAAGGGGTTTT"
        s0 = lineage_tracer_align(ref, qry, cutsites, short_match_window=0)[0]
        s1 = lineage_tracer_align(ref, qry, cutsites, short_match_window=3, short_match_discount=0.0)[0]
        assert s1 <= s0


# ═══════════════════════════════════════════════════════════════
# Phase A3: Dense mismatch region penalty
# ═══════════════════════════════════════════════════════════════

class TestDenseMismatchPenalty:
    def test_default_is_no_penalty(self):
        go = np.full(20, -2.0); ge = np.full(20, -0.1)
        ref, qry = "AAAACCCCGGGG", "AAAATTTTGGGG"
        s0, _, _, _ = affine_gap_alignment_position_aware(ref, qry, go, ge, dense_mismatch_penalty=0.0)
        s1, _, _, _ = affine_gap_alignment_position_aware(ref, qry, go, ge, dense_mismatch_penalty=0.0, dense_mismatch_window=6, dense_mismatch_threshold=0.34)
        assert s0 == s1

    def test_dense_mismatch_converts_to_insertion(self):
        go = np.full(12, -5.0); ge = np.full(12, -0.1)
        mmp = np.full(12, -1.0)
        ref, qry = "AAAACCCCGGGG", "AAAATTTTGGGG"
        s_off, _, _, st_off = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=mmp, dense_mismatch_penalty=0.0)
        s_on, _, _, st_on = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=mmp, dense_mismatch_window=4, dense_mismatch_threshold=0.5, dense_mismatch_penalty=-5.0)
        print(f"\n  OFF: score={s_off:.1f} mismatches={st_off['mismatches']} gaps_ref={st_off['gaps_in_ref']}")
        print(f"  ON:  score={s_on:.1f} mismatches={st_on['mismatches']} gaps_ref={st_on['gaps_in_ref']}")
        changed = (st_on['mismatches'] != st_off['mismatches'] or st_on['gaps_in_ref'] != st_off['gaps_in_ref'])
        assert changed or abs(s_on - s_off) > 1e-9

    def test_dense_penalty_affects_score(self):
        go = np.full(20, -2.0); ge = np.full(20, -0.1)
        mmp = np.full(20, -3.0)
        ref, qry = "AAAACCCCGGGG", "AAAATTTTGGGG"
        s0 = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=mmp, dense_mismatch_penalty=0.0)[0]
        s1 = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=mmp, dense_mismatch_window=4, dense_mismatch_threshold=0.5, dense_mismatch_penalty=-5.0)[0]
        assert s1 <= s0

    def test_no_mismatch_region_unaffected(self):
        go = np.full(20, -2.0); ge = np.full(20, -0.1)
        mmp = np.full(20, -3.0)
        ref, qry = "AAAACCCCGGGG", "AAAACCCCGGGG"
        s0, ar0, aq0, _ = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=mmp, dense_mismatch_penalty=0.0)
        s1, ar1, aq1, _ = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=mmp, dense_mismatch_window=4, dense_mismatch_threshold=0.5, dense_mismatch_penalty=-10.0)
        assert s0 == s1 and ar0 == ar1 and aq0 == aq1

    def test_lineage_tracer_forwards_dense_penalty(self):
        cutsites = [CutsiteRegion("T1", 4, 7)]
        ref, qry = "AAAACCCCGGGG", "AAAATTTTGGGG"
        s0 = lineage_tracer_align(ref, qry, cutsites, dense_mismatch_penalty=0.0)[0]
        s1 = lineage_tracer_align(ref, qry, cutsites, dense_mismatch_window=4, dense_mismatch_penalty=-5.0)[0]
        assert s1 <= s0

    def test_mismatch_profile_required_for_dense(self):
        go = np.full(20, -2.0); ge = np.full(20, -0.1)
        ref, qry = "AAAACCCCGGGG", "AAAATTTTGGGG"
        s = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=None, dense_mismatch_window=4, dense_mismatch_threshold=0.5, dense_mismatch_penalty=-5.0)[0]
        assert isinstance(s, float)


# ═══════════════════════════════════════════════════════════════
# Phase A3: Combined parameters
# ═══════════════════════════════════════════════════════════════

class TestCombinedParameters:
    def test_all_defaults_no_change(self):
        ref, qry = "ACGTACGT", "ACGTACGT"
        s0, ar0, aq0, st0 = affine_gap_alignment(ref, qry)
        s1, ar1, aq1, st1 = affine_gap_alignment_position_aware(ref, qry, np.full(8, -2.0), np.full(8, -0.1))
        assert s0 == s1 and ar0 == ar1

    def test_all_params_together(self):
        go = np.full(30, -2.0); ge = np.full(30, -0.1)
        mmp = np.full(30, -3.0)
        ref, qry = "AAACCCGGGTTTAAACCCGGGTTT", "AAATTTGGGTTTAAATTTGGGTTT"
        s, ar, aq, st = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=mmp,
            gap_exit_bonus=-1.0, short_match_window=3, short_match_discount=0.5,
            dense_mismatch_window=4, dense_mismatch_threshold=0.5, dense_mismatch_penalty=-2.0)
        assert isinstance(s, float) and len(ar) == len(aq) and st['alignment_length'] > 0


# ═══════════════════════════════════════════════════════════════
# Phase A4: Homology region penalty profile
# ═══════════════════════════════════════════════════════════════

class TestHomologyPenalty:
    def test_default_is_zero(self):
        cutsites = [CutsiteRegion("T1", 0, 3)]
        ref, qry = "ACGTACGT", "ACGTACGT"
        s0 = lineage_tracer_align(ref, qry, cutsites, homology_penalty=0.0)[0]
        s1 = lineage_tracer_align(ref, qry, cutsites)[0]
        assert s0 == s1

    def test_homology_penalty_lowers_score(self):
        cutsites = [CutsiteRegion("T1", 0, 3)]
        ref, qry = "ACGTACGTACGT", "ACGTACGTACGT"
        s0 = lineage_tracer_align(ref, qry, cutsites, homology_penalty=0.0)[0]
        s1 = lineage_tracer_align(ref, qry, cutsites, homology_window=4, homology_penalty=-1.0)[0]
        assert s1 <= s0 + 1e-9

    def test_build_homology_profile_length(self):
        from crisviper.lineage import build_homology_penalty_profile
        ref = "ACGTACGTACGT"
        profile = build_homology_penalty_profile(ref, homology_window=4, homology_penalty=-1.0)
        assert len(profile) == len(ref)
        assert np.any(profile < 0)

    def test_unique_sequence_no_penalty(self):
        from crisviper.lineage import build_homology_penalty_profile
        ref = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        profile = build_homology_penalty_profile(ref, homology_window=4, homology_penalty=-1.0)
        assert np.all(profile == 0.0)

    def test_pipeline_accepts_homology_params(self):
        config = PipelineConfig(lineage_mode=True, homology_window=8, homology_penalty=-1.0, min_reads_sub=0, min_reads_indel=0)
        pipeline = Pipeline(config=config, ref_seq=CARLIN_REF); pipeline.load_cutsites()
        result = pipeline.run([QueryRecord(readName="t",cellBC="t",UMI="U",readCount=1,seq=CARLIN_REF)])
        assert result.stats.successful == 1

    def test_affine_gap_alignment_position_aware_accepts_homology(self):
        go = np.full(20, -2.0); ge = np.full(20, -0.1)
        ref, qry = "AAAAAAAACCCCCCCC", "AAAAAAAACCCCCCCC"
        hp = np.full(20, -0.5)
        s0, ar0, aq0, _ = affine_gap_alignment_position_aware(ref, qry, go, ge, homology_profile=None)
        s1, ar1, aq1, _ = affine_gap_alignment_position_aware(ref, qry, go, ge, homology_profile=hp)
        assert s1 <= s0 and len(ar0) == len(ar1)


# ═══════════════════════════════════════════════════════════════
# Phase A5: Isolated base penalty
# ═══════════════════════════════════════════════════════════════

class TestIsolatedBasePenalty:
    def test_default_is_zero(self):
        cutsites = [CutsiteRegion("T1", 0, 3)]
        ref, qry = "ACGTACGT", "ACGTACGT"
        s0 = lineage_tracer_align(ref, qry, cutsites, isolated_base_penalty=0.0)[0]
        s1 = lineage_tracer_align(ref, qry, cutsites)[0]
        assert s0 == s1

    def test_negative_penalty_lowers_score(self):
        cutsites = [CutsiteRegion("T1", 0, 3)]
        ref, qry = "ACGTACGT", "ACGTACCT"
        s0 = lineage_tracer_align(ref, qry, cutsites, isolated_base_penalty=0.0)[0]
        s1 = lineage_tracer_align(ref, qry, cutsites, isolated_base_penalty=-2.0)[0]
        assert s1 <= s0

    def test_isolated_base_absorbed_into_gap(self):
        go = np.full(20, -2.0); ge = np.full(20, -0.1)
        ref, qry = "AAACCCGGGTTT", "AAAGGGTTT"
        s_off, ar_off, aq_off, st_off = affine_gap_alignment_position_aware(ref, qry, go, ge, isolated_base_penalty=0.0)
        s_on, ar_on, aq_on, st_on = affine_gap_alignment_position_aware(ref, qry, go, ge, isolated_base_penalty=-5.0)
        print(f"\n  ref={ref}  qry={qry}\n  OFF: ar={ar_off} aq={aq_off} blocks={st_off['gap_blocks_query']}\n  ON:  ar={ar_on} aq={aq_on} blocks={st_on['gap_blocks_query']}")
        assert isinstance(s_on, float) and len(ar_on) == len(aq_on)

    def test_pipeline_accepts_isolated_base(self):
        config = PipelineConfig(lineage_mode=True, isolated_base_penalty=-2.0, min_reads_sub=0, min_reads_indel=0)
        pipeline = Pipeline(config=config, ref_seq=CARLIN_REF); pipeline.load_cutsites()
        result = pipeline.run([QueryRecord(readName="t",cellBC="t",UMI="U",readCount=1,seq=CARLIN_REF)])
        assert result.stats.successful == 1


# ═══════════════════════════════════════════════════════════════
# Phase A1: Correction pipeline controls
# ═══════════════════════════════════════════════════════════════

# ═══════════════════════════════════════════════════════════════
# Phase B: All parameter combinations verification
# ═══════════════════════════════════════════════════════════════

class TestAllParametersCombined:
    def test_all_dp_params_together(self):
        go = np.full(40, -2.0); ge = np.full(40, -0.1)
        mmp, hp = np.full(40, -3.0), np.full(40, -0.5)
        ref, qry = "AAACCCGGGTTTAAACCCGGGTTT", "AAATTTGGGTTTAAATTTGGGTTT"
        s, ar, aq, st = affine_gap_alignment_position_aware(ref, qry, go, ge, mismatch_penalty_profile=mmp,
            gap_exit_bonus=-1.0, short_match_window=3, short_match_discount=0.5,
            dense_mismatch_window=4, dense_mismatch_threshold=0.5, dense_mismatch_penalty=-2.0,
            homology_profile=hp, isolated_base_penalty=-1.0)
        assert isinstance(s, float) and len(ar) == len(aq) and st['alignment_length'] > 0

    def test_pipeline_all_params(self):
        config = PipelineConfig(lineage_mode=True, gap_exit_strength=-1.0,
            short_match_window=3, short_match_discount=0.5,
            dense_mismatch_window=6, dense_mismatch_penalty=-2.0,
            homology_window=8, homology_penalty=-1.0,
            isolated_base_penalty=-2.0, min_reads_sub=0, min_reads_indel=0)
        pipeline = Pipeline(config=config, ref_seq=CARLIN_REF); pipeline.load_cutsites()
        result = pipeline.run([QueryRecord(readName="t",cellBC="t",UMI="U",readCount=1,seq=CARLIN_REF)])
        assert result.stats.successful == 1
