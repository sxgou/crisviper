"""crisviper/pipeline.py — Full analysis pipeline orchestration.

Organizes the complete lineage tracing analysis workflow into a
configurable Pipeline class with the following steps:
  1. Full-length global alignment (standard Gotoh or lineage-tracer)
  2. Primer region quality check
  3. Internal region extraction (trim primers)
  4. Background substitution correction (before allele merging)
  5. WT primer assembly for full-length output
  6. Internal region re-scoring
  7. Mutation extraction
  8. Allele confidence filtering

Each step is independently testable, replaceable, and skippable.
"""

import os
# Prevent fork + NumPy thread conflicts (must set before importing numpy)
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from functools import partial
from typing import List, Dict, Optional, Tuple

from crisviper.logging_config import get_logger
from crisviper.models import (
    QueryRecord, AlignmentResult, AlignmentStats,
    PipelineConfig, PipelineResult, PipelineStats,
)
from crisviper.config import CutsiteRegion
from crisviper.mutation import extract_mutations, build_mutation_summary, _build_ref_pos_map_full
from crisviper.alignment import (
    affine_gap_alignment,
    affine_gap_alignment_position_aware,
    calculate_alignment_stats,
)
from crisviper.lineage import (
    lineage_tracer_align,
    build_gradient_profiles,
    build_homology_penalty_profile,
    get_amplicon_structure,
)
from crisviper.caller import call_alleles_coarse_grain, call_alleles_exact

log = get_logger(__name__)


# ═══════════════════════════════════════════════════════════════
# Step 1: Primer quality check (post-alignment)
# ═══════════════════════════════════════════════════════════════

def _check_primer_quality(
    ar: str,
    aq: str,
    ref_seq: str,
    p5: int,
    p3: int,
    p5_threshold: int,
    p3_threshold: int,
) -> Tuple[bool, int, bool, int]:
    """Check primer region alignment quality from the full-length global alignment.

    Called after global alignment to evaluate whether the query sequence
    matches the primer regions well enough. Replaces the older
    check_primer_anchoring pre-filter approach.

    Args:
        ar: Aligned reference sequence (with gaps).
        aq: Aligned query sequence (with gaps).
        ref_seq: Original (ungapped) reference sequence.
        p5: 5' primer length in bp.
        p3: 3' primer length in bp.
        p5_threshold: Minimum required matches for the 5' primer.
        p3_threshold: Minimum required matches for the 3' primer.

    Returns:
        Tuple of (p5_passed, p5_match_count, p3_passed, p3_match_count).
    """
    pos_map, total = _build_ref_pos_map_full(ar)
    if total < p5 + p3:
        return False, 0, False, 0

    # 5' primer region: ref positions 0..p5-1
    p5_end = next((i for i, p in enumerate(pos_map) if p == p5 - 1), None)
    if p5_end is None:
        return False, 0, False, 0
    p5_match = sum(1 for i in range(p5_end + 1)
                   if ar[i] != '-' and aq[i] != '-' and ar[i] == aq[i])

    # 3' primer region: ref positions total-p3..total-1
    # The trailing columns beyond total-1 are excluded by the ar[i] != '-' and
    # aq[i] != '-' filters (Iy termination → aq gap, Ix termination → ar gap).
    p3_start = next((i for i, p in enumerate(pos_map) if p == total - p3), None)
    if p3_start is None:
        return False, p5_match, False, 0
    p3_match = sum(1 for i in range(p3_start, len(ar))
                   if ar[i] != '-' and aq[i] != '-' and ar[i] == aq[i])

    return p5_match >= p5_threshold, p5_match, p3_match >= p3_threshold, p3_match


def _extract_internal_region(
    ar: str,
    aq: str,
    r_seq: str,
    p5: int,
    p3: int,
) -> Tuple[Optional[str], Optional[str], str]:
    """Extract the primer-trimmed internal region from a full-length alignment.

    Removes the 5' and 3' primer columns from the aligned sequences,
    returning only the internal (target-containing) region for downstream
    mutation analysis.

    Args:
        ar: Aligned reference sequence (full-length, with gaps).
        aq: Aligned query sequence (full-length, with gaps).
        r_seq: Original (ungapped) reference sequence.
        p5: 5' primer length in bp.
        p3: 3' primer length in bp.

    Returns:
        Tuple of (aligned_internal_ref, aligned_internal_query, ungapped_internal_ref).
        Returns (None, None, "") if the internal region cannot be extracted.
    """
    pos_map, total = _build_ref_pos_map_full(ar)

    int_start = next((i for i, p in enumerate(pos_map) if p == p5), None)
    int_end = next((i for i, p in enumerate(pos_map) if p == total - p3 - 1), None)
    if int_start is None or int_end is None:
        return None, None, ""

    ar_int = ar[int_start:int_end + 1]
    aq_int = aq[int_start:int_end + 1]
    int_r = r_seq[p5:total - p3] if p3 > 0 else r_seq[p5:]

    return ar_int, aq_int, int_r


# ═══════════════════════════════════════════════════════════════
# Step 2: Allele confidence filtering
# ═══════════════════════════════════════════════════════════════

def check_allele_confidence(
    stats: AlignmentStats,
    read_count: int,
    min_reads_sub: int = 5,
    min_reads_indel: int = 0,
) -> Tuple[bool, str]:
    """Check whether an allele passes confidence filtering thresholds.

    Rules:
      - Substitution-only (no indel): requires readCount >= min_reads_sub.
      - Indel-containing: requires readCount >= min_reads_indel.
      - Wild-type (no mutations): passes automatically.

    Args:
        stats: Alignment statistics for this allele.
        read_count: Read count supporting this allele.
        min_reads_sub: Minimum read count for pure-substitution alleles (inclusive).
        min_reads_indel: Minimum read count for indel-containing alleles (inclusive).

    Returns:
        Tuple of (passed: bool, reason: str). Empty reason means passed.
    """
    if not stats.has_indel and stats.mismatches > 0:
        # Substitution-only: must reach threshold
        if read_count < min_reads_sub:
            return False, f"False positive: substitution-only with insufficient reads ({read_count}<{min_reads_sub})"
    elif stats.has_indel:
        # Indel-containing: must reach threshold
        if read_count < min_reads_indel:
            return False, f"False positive: indel with insufficient reads ({read_count}<{min_reads_indel})"
    return True, ""


# ═══════════════════════════════════════════════════════════════
# Step 3: Background substitution correction
# ═══════════════════════════════════════════════════════════════

def correct_background_substitutions(
    ar_int: str,
    aq_int: str,
    cutsites: Optional[List[CutsiteRegion]],
    sub_window: int = 3,
    keep_sub_indel_window: int = 3,
    lineage_mode: bool = True,
) -> str:
    """Correct background sequencing errors by reverting out-of-context substitutions.

    Called before allele merging to eliminate allele fragmentation caused
    by sequencing errors that appear as isolated point mutations outside
    expected cut sites.

    A substitution is kept (not reverted) if it falls within:
      1. Lineage mode: cutsite center ± sub_window
      2. Standard mode: the cutsite region itself
      3. Any mode: within keep_sub_indel_window of an indel endpoint

    Args:
        ar_int: Internal-region aligned reference (with gaps).
        aq_int: Internal-region aligned query (with gaps).
        cutsites: Cutsite regions (already adjusted to internal coordinates).
        sub_window: Cutsite-adjacent retention window (bp).
        keep_sub_indel_window: Indel-adjacent retention window (bp).
        lineage_mode: Whether lineage-tracer mode is active.

    Returns:
        Corrected aq_int string with background substitutions reverted to reference.
    """
    if cutsites is None or not ar_int or len(ar_int) != len(aq_int):
        return aq_int

    pos_map, total = _build_ref_pos_map_full(ar_int)

    # Step 1: Collect indel-affected ref coordinate intervals (with flank extension)
    indel_keep_regions = []
    i = 0
    alen = len(ar_int)
    while i < alen:
        # Deletion: ref has base, query is gap
        if ar_int[i] != '-' and aq_int[i] == '-':
            start = pos_map[i]
            while i < alen and ar_int[i] != '-' and aq_int[i] == '-':
                i += 1
            end = pos_map[i - 1]
            rs = max(0, start - keep_sub_indel_window)
            re = min(total - 1, end + keep_sub_indel_window)
            indel_keep_regions.append((rs, re))
        # Insertion: ref is gap, query has base
        elif ar_int[i] == '-' and aq_int[i] != '-':
            left_pos = pos_map[i - 1] if i > 0 and pos_map[i - 1] >= 0 else 0
            while i < alen and ar_int[i] == '-' and aq_int[i] != '-':
                i += 1
            rs = max(0, left_pos - keep_sub_indel_window)
            re = min(total - 1, left_pos + keep_sub_indel_window)
            indel_keep_regions.append((rs, re))
        else:
            i += 1

    # Step 2: Build retention region set (merged intervals)
    keep_intervals = []  # (start, end) inclusive

    # 2a: Cutsite regions
    for cs in cutsites:
        if lineage_mode:
            ks = max(0, cs.start - sub_window)
            ke = min(total - 1, cs.end + sub_window)
        else:
            ks = max(0, cs.start)
            ke = min(total - 1, cs.end)
        keep_intervals.append((ks, ke))

    # 2b: Indel flanking regions
    keep_intervals.extend(indel_keep_regions)

    if not keep_intervals:
        return aq_int

    # Merge overlapping intervals
    keep_intervals.sort()
    merged = [keep_intervals[0]]
    for s, e in keep_intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    # Step 3: Scan substitution columns and revert those outside retention intervals
    aq_list = list(aq_int)
    for i in range(alen):
        if ar_int[i] != '-' and aq_int[i] != '-' and ar_int[i] != aq_int[i]:
            ref_pos = pos_map[i]
            if ref_pos < 0:
                continue
            in_keep = False
            for ks, ke in merged:
                if ks <= ref_pos <= ke:
                    in_keep = True
                    break
            if not in_keep:
                aq_list[i] = ar_int[i]

    return ''.join(aq_list)


# ═══════════════════════════════════════════════════════════════
# Step 4: Single-sequence alignment pipeline
# ═══════════════════════════════════════════════════════════════


def align_single(
    query: QueryRecord,
    ref_seq: str,
    config: PipelineConfig,
    cutsites: Optional[List[CutsiteRegion]] = None,
) -> AlignmentResult:
    """Run the complete per-sequence alignment pipeline.

    Pipeline steps:
      1. Full-length global alignment (standard Gotoh or lineage-tracer mode)
      2. Primer region quality check (replaces older pre-filter approach)
      3. Internal region extraction (trim primer columns)
      4. Background substitution correction (eliminate sequencing noise)
      5. WT primer assembly for full-length output
      6. Internal region re-scoring
      7. Mutation extraction
      8. Allele confidence filtering

    Args:
        query: The query sequence record.
        ref_seq: Full reference sequence.
        config: Pipeline configuration.
        cutsites: Cutsite regions (optional, required for gradient mode).

    Returns:
        AlignmentResult with success/failure status and all extracted data.
    """
    query_seq = query.seq
    q_seq = query_seq
    r_seq = ref_seq
    p5 = config.primer5_len
    p3 = config.primer3_len

    # ── Step 1: Full-length global alignment ──
    if config.lineage_mode and cutsites:
        score, ar, aq, raw_stats = _align_full_lineage(
            r_seq, q_seq, cutsites, config,
        )
    else:
        score, ar, aq, raw_stats = _align_full_standard(
            r_seq, q_seq, config, cutsites=cutsites,
        )

    if not ar:
        return AlignmentResult.error_result(query, "Full-length alignment failed", category="alignment")

    # ── Step 2: Primer region alignment quality check (post-alignment) ──
    p5_ok, p5_match, p3_ok, p3_match = _check_primer_quality(
        ar, aq, r_seq, p5, p3,
        config.primer5_threshold, config.primer3_threshold,
    )
    if not (p5_ok and p3_ok):
        n5 = config.primer5_threshold
        n3 = config.primer3_threshold
        return AlignmentResult.error_result(
            query,
            f"Primer anchoring failed: Primer5({p5_match}/{p5}<{n5}) "
            f"Primer3({p3_match}/{p3}<{n3})",
            category="anchor",
        )

    # ── Step 3: Extract internal (primer-trimmed) region ──
    ar_int, aq_int, int_r = _extract_internal_region(ar, aq, r_seq, p5, p3)
    if ar_int is None:
        return AlignmentResult.error_result(query, "Failed to extract internal region from full-length alignment", category="extraction")

    # ── Step 4: Background substitution correction (pre-allele-merging) ──
    if config.correct_bg_sub and cutsites is not None:
        int_r_len = len(r_seq) - p5 - p3
        internal_cutsites = _adjust_cutsites_to_internal(cutsites, p5, int_r_len)
        aq_int = correct_background_substitutions(
            ar_int, aq_int,
            cutsites=internal_cutsites,
            sub_window=config.sub_window,
            keep_sub_indel_window=config.keep_sub_indel_window,
            lineage_mode=config.lineage_mode,
        )

    # ── Step 4b: Merge DEL→INS→DEL artifact patterns ──
    ar_int, aq_int = _merge_del_ins_del(ar_int, aq_int)

    # ── Step 5: WT primer assembly for full-length output ──
    aligned_ref, aligned_query = _assemble_full_length(
        ar_int, aq_int, r_seq, p5, p3,
    )

    # ── Step 6: Internal region re-scoring ──
    stats_dict = calculate_alignment_stats(ar_int, aq_int)
    score = (stats_dict['matches'] * 2 +
             stats_dict['mismatches'] * (-3) +
             stats_dict['gaps_in_ref'] * (-2) +
             stats_dict['gaps_in_query'] * (-2))
    stats_dict['score'] = score
    stats = AlignmentStats.from_dict(stats_dict)

    # ── Step 7: Mutation extraction (internal region only) ──
    if cutsites is not None:
        int_r_len = len(r_seq) - p5 - p3
        internal_cutsites = _adjust_cutsites_to_internal(cutsites, p5, int_r_len)
    else:
        internal_cutsites = None
    mutations = extract_mutations(
        ar_int, aq_int,
        cutsites=internal_cutsites,
        sub_window=config.sub_window,
    )

    # Convert ref_pos from internal coordinates back to full-length reference coordinates
    for m in mutations:
        m.ref_pos += config.primer5_len

    # ── Step 8: Allele confidence filtering ──
    passed, reason = check_allele_confidence(
        stats, query.readCount,
        min_reads_sub=config.min_reads_sub,
        min_reads_indel=config.min_reads_indel,
    )
    if not passed:
        return AlignmentResult.error_result(query, reason, category="noise")

    return AlignmentResult(
        query=query,
        success=True,
        score=score,
        aligned_ref=aligned_ref,
        aligned_query=aligned_query,
        stats=stats,
        mutations=mutations,
        mode="lineage" if config.lineage_mode else "standard",
    )


# ═══════════════════════════════════════════════════════════════
# Lineage tracer alignment mode
# ═══════════════════════════════════════════════════════════════

def _align_full_lineage(
    r_seq: str,
    q_seq: str,
    cutsites: List[CutsiteRegion],
    config: PipelineConfig,
) -> Tuple[float, str, str, Dict]:
    """Full-length lineage-tracer alignment of ref vs query (332bp global alignment).

    Delegates to lineage_tracer_align for the full-length alignment.
    Cutsite coordinates are in full-length reference coordinates.
    No special CGCCG prefix handling is needed — the DP naturally handles
    context in the 332bp alignment window.
    """
    if not q_seq:
        return 0.0, "", "", {"matches": 0, "mismatches": 0, "gaps_in_ref": 0,
                             "gaps_in_query": 0, "gap_blocks_ref": [],
                             "gap_blocks_query": [], "alignment_length": 0,
                             "similarity": 0.0, "identity": 0.0, "score": 0.0}

    score, ar, aq, raw_stats = lineage_tracer_align(
        r_seq, q_seq, cutsites,
        match_score=config.match_score,
        mismatch_penalty=config.mismatch_penalty,
        base_gap_open=config.gap_open,
        base_gap_extend=config.gap_extend,
        min_scale=config.min_scale,
        max_scale=config.max_scale,
        cutsite_edge_scale=config.cutsite_edge_scale,
        gradient_radius=config.gradient_radius,
        mismatch_density_threshold=config.mismatch_density_threshold,
        sub_window=config.sub_window,
        gap_exit_bonus=0.0,
        base_gap_exit=config.gap_exit_strength,
        short_match_window=config.short_match_window,
        short_match_discount=config.short_match_discount,
        dense_mismatch_window=config.dense_mismatch_window,
        dense_mismatch_penalty=config.dense_mismatch_penalty,
        homology_window=config.homology_window,
        homology_penalty=config.homology_penalty,
        isolated_base_penalty=config.isolated_base_penalty,
    )
    return score, ar, aq, raw_stats


# ═══════════════════════════════════════════════════════════════
# Standard Gotoh alignment mode
# ═══════════════════════════════════════════════════════════════

def _align_full_standard(
    r_seq: str,
    q_seq: str,
    config: PipelineConfig,
    cutsites: Optional[List[CutsiteRegion]] = None,
) -> Tuple[float, str, str, Dict]:
    """Full-length standard Gotoh alignment of ref vs query.

    When cutsites are available and config.gradient_mode is enabled,
    uses smoothstep position-aware gap/mismatch penalties; otherwise
    uses global fixed penalties.
    """
    if not q_seq:
        return 0.0, "", "", {"matches": 0, "mismatches": 0, "gaps_in_ref": 0,
                             "gaps_in_query": 0, "gap_blocks_ref": [],
                             "gap_blocks_query": [], "alignment_length": 0,
                             "similarity": 0.0, "identity": 0.0, "score": 0.0}

    if cutsites and config.gradient_mode:
        gap_open_p, gap_extend_p, mismatch_p, _ = build_gradient_profiles(
            ref_length=len(r_seq), cutsites=cutsites,
            base_gap_open=config.gap_open, base_gap_extend=config.gap_extend,
            mismatch_penalty=config.mismatch_penalty,
            min_scale=config.min_scale, max_scale=config.max_scale,
            cutsite_edge_scale=config.cutsite_edge_scale,
            gradient_radius=config.gradient_radius,
        )
        homology_p = build_homology_penalty_profile(
            r_seq, homology_window=config.homology_window,
            homology_penalty=config.homology_penalty,
        )
        score, ar, aq, raw_stats = affine_gap_alignment_position_aware(
            r_seq, q_seq, gap_open_p, gap_extend_p,
            match_score=config.match_score,
            mismatch_penalty=config.mismatch_penalty,
            mismatch_penalty_profile=mismatch_p,
            gap_exit_bonus=config.gap_exit_strength,
            short_match_window=config.short_match_window,
            short_match_discount=config.short_match_discount,
            dense_mismatch_window=config.dense_mismatch_window,
            dense_mismatch_threshold=config.mismatch_density_threshold,
            dense_mismatch_penalty=config.dense_mismatch_penalty,
            homology_profile=homology_p,
            isolated_base_penalty=config.isolated_base_penalty,
        )
    else:
        score, ar, aq, raw_stats = affine_gap_alignment(
            r_seq, q_seq,
            match_score=config.match_score,
            mismatch_penalty=config.mismatch_penalty,
            gap_open=config.gap_open,
            gap_extend=config.gap_extend,
            gap_exit_bonus=config.gap_exit_strength,
            short_match_window=config.short_match_window,
            short_match_discount=config.short_match_discount,
            isolated_base_penalty=config.isolated_base_penalty,
        )
    return score, ar, aq, raw_stats


# ═══════════════════════════════════════════════════════════════
# DEL→INS→DEL artifact correction
# ═══════════════════════════════════════════════════════════════

def _merge_del_ins_del(ar_int: str, aq_int: str) -> Tuple[str, str]:
    """Merge DEL→INS→DEL artifact patterns in aligned strings.

    The position-dependent gradient DP can split a single INDEL event into
    DEL→INS→DEL when gap costs differ sharply across the cutsite boundary.
    This detects such patterns and re-arranges to contiguous DEL followed
    by INS, matching the biological reality of NHEJ repair.

    Args:
        ar_int: Internal-region aligned reference (may contain gaps).
        aq_int: Internal-region aligned query (may contain gaps).

    Returns:
        Corrected (ar_int, aq_int) — unchanged if no pattern found.
    """
    alen = len(ar_int)
    if alen == 0 or len(ar_int) != len(aq_int):
        return ar_int, aq_int

    result_ar = list(ar_int)
    result_aq = list(aq_int)
    changed = False

    i = 0
    while i < alen:
        # Block 1: deletion (ref base, query gap)
        if ar_int[i] != '-' and aq_int[i] == '-':
            b1_start = i
            while i < alen and ar_int[i] != '-' and aq_int[i] == '-':
                i += 1
            b1_end = i
            b1_len = b1_end - b1_start

            # Block 2: insertion (ref gap, query base) — must immediately follow
            if i < alen and ar_int[i] == '-' and aq_int[i] != '-':
                b2_start = i
                while i < alen and ar_int[i] == '-' and aq_int[i] != '-':
                    i += 1
                b2_end = i
                b2_len = b2_end - b2_start

                # Block 3: deletion (ref base, query gap) — must immediately follow
                if i < alen and ar_int[i] != '-' and aq_int[i] == '-':
                    b3_start = i
                    while i < alen and ar_int[i] != '-' and aq_int[i] == '-':
                        i += 1
                    b3_end = i
                    b3_len = b3_end - b3_start
                    total_del = b1_len + b3_len

                    # Extract bases to re-arrange
                    ref_b1 = ar_int[b1_start:b1_end]
                    ref_b3 = ar_int[b3_start:b3_end]
                    q_b2 = aq_int[b2_start:b2_end]

                    # Rearrange: [combined ref] + [insertion gaps in ref]
                    #            [combined del gaps in query] + [insertion bases]
                    new_ar_seg = ref_b1 + ref_b3 + '-' * b2_len
                    new_aq_seg = '-' * total_del + q_b2

                    result_ar[b1_start:b3_end] = list(new_ar_seg)
                    result_aq[b1_start:b3_end] = list(new_aq_seg)
                    changed = True

                    # Continue scanning after the merged region
                    i = b3_end
                    continue  # skip the else clause
        else:
            i += 1

    if not changed:
        return ar_int, aq_int
    return ''.join(result_ar), ''.join(result_aq)


# ═══════════════════════════════════════════════════════════════
# Full-length alignment assembly
# ═══════════════════════════════════════════════════════════════

def _assemble_full_length(
    ar_int: str,
    aq_int: str,
    r_seq: str,
    p5: int,
    p3: int,
) -> Tuple[str, str]:
    """Assemble full-length alignment using wild-type primer sequences.

    Replaces the query's primer regions with the reference (WT) primer
    sequences to avoid allele fragmentation caused by primer-region
    mutations in the query.
    """
    aligned_ref = r_seq[:p5] + ar_int + (r_seq[-p3:] if p3 > 0 else "")
    aligned_query = r_seq[:p5] + aq_int + (r_seq[-p3:] if p3 > 0 else "")
    return aligned_ref, aligned_query


# ═══════════════════════════════════════════════════════════════
# Cutsite coordinate adjustment
# ═══════════════════════════════════════════════════════════════

def _adjust_cutsites_to_internal(
    cutsites: List[CutsiteRegion],
    primer5_len: int,
    int_r_len: int,
) -> List[CutsiteRegion]:
    """Convert cutsite coordinates from full-length to internal (primer-trimmed) coordinates."""
    result = []
    for cs in cutsites:
        ns = cs.start - primer5_len
        ne = cs.end - primer5_len
        if ns < 0 and ne < 0:
            continue
        if ne >= int_r_len:
            ne = int_r_len - 1
        if ns < 0:
            ns = 0
        if ns <= ne:
            result.append(CutsiteRegion(name=cs.name, start=ns, end=ne))
    return result


# ═══════════════════════════════════════════════════════════════
# Chunk processing (for parallel execution)
# ═══════════════════════════════════════════════════════════════

def _process_chunk(
    queries_chunk: List[QueryRecord],
    ref_seq: str,
    config: PipelineConfig,
    cutsites: Optional[List[CutsiteRegion]] = None,
) -> List[AlignmentResult]:
    """Process a chunk of query sequences sequentially (for parallel executor pool calls)."""
    return [
        align_single(q, ref_seq, config, cutsites)
        for q in queries_chunk
    ]


# ═══════════════════════════════════════════════════════════════
# Pipeline main class
# ═══════════════════════════════════════════════════════════════

class Pipeline:
    """Lineage tracing analysis pipeline orchestrator.

    Full workflow:
      1. Full-length global alignment        (affine_gap_alignment / lineage_tracer_align)
      2. Primer region alignment quality     (_check_primer_quality)
      3. Internal region extraction          (_extract_internal_region)
      4. Background substitution correction  (correct_background_substitutions)
      5. WT primer assembly                  (_assemble_full_length)
      6. Internal region re-scoring          (calculate_alignment_stats)
      7. Mutation extraction                 (extract_mutations)
      8. Allele confidence filtering         (check_allele_confidence)

    Each step is independently testable and replaceable.
    """

    def __init__(
        self,
        config: Optional[PipelineConfig] = None,
        ref_seq: str = "",
    ):
        self.config = config or PipelineConfig()
        self.ref_seq = ref_seq
        self.cutsites: Optional[List[CutsiteRegion]] = None

    def load_cutsites(self, ref_seq: Optional[str] = None) -> None:
        """Load or auto-detect cutsite positions.

        Priority (highest to lowest):
        1. explicit_cutsites — explicit cutsite list from YAML config
        2. cutsites_path — JSON file with cutsite definitions
        3. auto_detect — auto-detection using amplicon_config (if available)
        """
        seq = ref_seq or self.ref_seq
        if not seq:
            return

        # Priority 1: Explicit cutsites from YAML configuration
        if self.config.explicit_cutsites:
            self.cutsites = self.config.explicit_cutsites
            log.info("Loaded %d cutsite regions from YAML config", len(self.cutsites))
            return

        if self.config.cutsites_path:
            # Load from JSON file
            import json
            try:
                with open(self.config.cutsites_path) as f:
                    cs_data = json.load(f)
                raw = cs_data.get("cutsites", cs_data)
                self.cutsites = [
                    CutsiteRegion(name=c.get("name", f"Target{i+1}"),
                                  start=c["start"], end=c["end"])
                    for i, c in enumerate(raw)
                ]
                log.info("Loaded %d cutsite regions from config file", len(self.cutsites))
            except Exception as e:
                log.error("Failed to read cutsite config file - %s", e)
                sys.exit(1)
        elif self.config.auto_detect_cutsites:
            self.cutsites = get_amplicon_structure(seq, config=self.config.amplicon_config)
            if self.cutsites:
                log.info("Auto-detected %d cutsite regions", len(self.cutsites))
            else:
                log.warning("Could not auto-detect cutsite positions")

    def run(
        self,
        queries: List[QueryRecord],
        ref_seq: Optional[str] = None,
    ) -> PipelineResult:
        """Run the complete analysis pipeline on all query sequences.

        Args:
            queries: List of QueryRecord objects to align.
            ref_seq: Reference sequence (overrides constructor value if set).

        Returns:
            PipelineResult containing all alignment results and statistics.
        """
        if ref_seq:
            self.ref_seq = ref_seq
        if not self.ref_seq:
            raise ValueError("Reference sequence not set")

        # Auto-load cutsites when lineage mode or gradient mode is enabled
        # (standard mode with gradient_mode=True requires cutsites for
        # position-aware penalty profiles)
        if self.cutsites is None and (self.config.lineage_mode or self.config.gradient_mode):
            self.load_cutsites()

        total = len(queries)
        if total == 0:
            return PipelineResult(
                results=[], config=self.config,
                stats=PipelineStats(),
                ref_length=len(self.ref_seq),
            )

        # ── Parallel alignment ──
        mode_label = "lineage-tracer" if self.config.lineage_mode else "standard"

        threads = self.config.threads or 1
        threads = min(threads, total)

        log.info("  Processing %d sequences across %d parallel processes (%s mode)...",
                 total, threads, mode_label)

        # Compute chunk size for batching
        chunk_size = self.config.chunk_size
        chunks = [queries[i:i + chunk_size] for i in range(0, total, chunk_size)]
        log.info("  Split into %d chunks (max %d per chunk)", len(chunks), chunk_size)

        # Parallel execution with ProcessPoolExecutor
        results = []

        processed_indices = set()
        if threads > 1:
            try:
                with ProcessPoolExecutor(max_workers=threads) as executor:
                    chunk_func = partial(
                        _process_chunk,
                        ref_seq=self.ref_seq,
                        config=self.config,
                        cutsites=self.cutsites,
                    )
                    future_to_idx = {executor.submit(chunk_func, ch): i for i, ch in enumerate(chunks)}
                    for future in as_completed(future_to_idx):
                        idx = future_to_idx[future]
                        try:
                            results.extend(future.result())
                            processed_indices.add(idx)
                        except BrokenProcessPool:
                            log.error("  Chunk %d: process pool broken, will reprocess", idx)
                            break
                        except Exception as e:
                            log.error("  Chunk %d failed, will reprocess: %s", idx, e)
            except BrokenProcessPool:
                log.warning("ProcessPoolExecutor crashed (BrokenProcessPool), reprocessing remaining chunks")

        # Re-process remaining chunks after BrokenProcessPool (also handles threads <= 1)
        remaining_indices = set(range(len(chunks))) - processed_indices
        if remaining_indices:
            if threads > 1 and processed_indices:
                log.info("  Reprocessing %d remaining chunks...", len(remaining_indices))
            else:
                log.info("  Processing %d sequences in single-threaded mode...", total)
            for i in remaining_indices:
                results.extend(_process_chunk(chunks[i], self.ref_seq, self.config, self.cutsites))

        # ── Allele calling (optional) ──
        called_alleles = []
        successful_results = [r for r in results if r.success]
        if self.config.call_alleles_enabled and successful_results:
            method = call_alleles_coarse_grain if self.config.call_alleles_mode == "coarse" else call_alleles_exact
            called_alleles = method(successful_results, dominant_frac=self.config.dominant_frac)
            log.info("Allele calling (%s): %d alleles identified",
                     self.config.call_alleles_mode, len(called_alleles))

        # ── Build pipeline statistics ──
        return self._build_pipeline_result(results, called_alleles)

    def _build_pipeline_result(self, results: List[AlignmentResult],
                               called_alleles: list = None) -> PipelineResult:
        """Build pipeline statistics from alignment results."""
        stats = PipelineStats()
        stats.total_queries = len(results)

        # Categorize discard reasons using failure_category field
        n_anchor = sum(1 for r in results if not r.success and r.failure_category == "anchor")
        n_noise = sum(1 for r in results if not r.success and r.failure_category == "noise")

        # Tally successful and failed
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        stats.successful = len(successful)
        stats.failed = len(failed)
        stats.total_reads = sum(r.query.readCount for r in successful)
        stats.n_anchor_failed = n_anchor
        stats.n_noise_filtered = n_noise

        # Mutation counts
        mutated = [r for r in successful if r.stats and r.stats.has_mutation]
        unmutated = [r for r in successful if r.stats and not r.stats.has_mutation]
        stats.mutated_sequences = len(mutated)
        stats.unmutated_sequences = len(unmutated)
        stats.mutated_reads = sum(r.query.readCount for r in mutated)

        # Build mutation summary
        summary = build_mutation_summary(successful)

        # Log discard reasons
        if stats.failed > 0:
            parts = []
            if n_anchor:
                parts.append(f"primer anchoring failed {n_anchor}")
            if n_noise:
                parts.append(f"false positive allele {n_noise}")
            other = stats.failed - n_anchor - n_noise
            if other:
                parts.append(f"other {other}")
            log.warning("Discarded %d sequences: %s", stats.failed, ", ".join(parts))
        log.info("Batch alignment complete: %d successful results", stats.successful)

        return PipelineResult(
            results=results,
            config=self.config,
            stats=stats,
            ref_length=len(self.ref_seq),
            mutation_type_counts=summary["type_counts"],
            total_mismatches=summary["total_mismatches"],
            insertion_lengths=summary["ins_lengths"],
            deletion_lengths=summary["del_lengths"],
            called_alleles=called_alleles or [],
        )
