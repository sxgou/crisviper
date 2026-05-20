"""Lineage tracer alignment pipeline — structure-aware alignment for CARLIN amplicons."""

import numpy as np
import re
from typing import Tuple, Dict, List, Optional
from crisviper.config import CutsiteRegion, AmpliconConfig
from crisviper.alignment import affine_gap_alignment_position_aware, calculate_alignment_stats
from crisviper.logging_config import get_logger

log = get_logger(__name__)


def _smoothstep(t: np.ndarray) -> np.ndarray:
    """Smoothstep easing: 3t² - 2t³, derivative is 0 at both t=0 and t=1."""
    return t * t * (3.0 - 2.0 * t)


def _calc_gradient_radius(cutsites: List[CutsiteRegion]) -> float:
    """Auto-calculate gradient radius from adjacent cutsite spacing.

    Returns half the distance between adjacent cutsite centers minus the
    cutsite half-width, so that the gradient reaches max_scale exactly at
    the midpoint between neighboring cutsites.
    """
    if len(cutsites) < 2:
        return 30.0
    centers = [(cs.start + cs.end) / 2.0 for cs in cutsites]
    min_gap = min(centers[i + 1] - centers[i] for i in range(len(centers) - 1))
    half_cs = (cutsites[0].end - cutsites[0].start + 1) / 2.0
    radius = min_gap / 2.0 - half_cs
    return max(radius, 1.0)


def build_gradient_profiles(
    ref_length: int,
    cutsites: List[CutsiteRegion],
    base_gap_open: float = -2.0,
    base_gap_extend: float = -0.1,
    mismatch_penalty: float = -3.0,
    min_scale: float = 1.0,
    max_scale: float = 6.0,
    cutsite_edge_scale: float = 2.0,
    gradient_radius: Optional[float] = None,
    base_gap_exit: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build position-dependent gap and mismatch penalty profiles using smoothstep gradient.

    For each cutsite, the penalty scale follows:
      - Internal micro-gradient: min_scale at cut center → cutsite_edge_scale at edge
      - External gradient: cutsite_edge_scale at edge → max_scale at gradient_radius away

    Multiple cutsites are combined by min() — the nearest cutsite dominates.

    Returns: (gap_open_profile, gap_extend_profile, mismatch_profile, geb_profile)
      All four are 1-D arrays of length ref_length.
    """
    if gradient_radius is None:
        gradient_radius = _calc_gradient_radius(cutsites)

    effective_scale = np.full(ref_length, max_scale, dtype=float)

    for cs in cutsites:
        center = (cs.start + cs.end) / 2.0
        cs_half = (cs.end - cs.start + 1) / 2.0

        positions = np.arange(ref_length, dtype=float)
        d = np.abs(positions - center)

        # Internal micro-gradient (within cutsite)
        internal = d <= cs_half
        if internal.any():
            t = d[internal] / cs_half
            scale = min_scale + (cutsite_edge_scale - min_scale) * _smoothstep(t)
            effective_scale[internal] = np.minimum(effective_scale[internal], scale)

        # External gradient (from cutsite edge to gradient_radius)
        outer = (d > cs_half) & (d <= cs_half + gradient_radius)
        if outer.any():
            t = (d[outer] - cs_half) / gradient_radius
            scale = cutsite_edge_scale + (max_scale - cutsite_edge_scale) * _smoothstep(t)
            effective_scale[outer] = np.minimum(effective_scale[outer], scale)

    gap_open_profile = base_gap_open * effective_scale
    gap_extend_profile = base_gap_extend * effective_scale
    mismatch_profile = mismatch_penalty * effective_scale
    # geb_profile: inverted gradient scaling — strongest penalty at cutsite center
    if base_gap_exit != 0.0:
        geb_profile = base_gap_exit * (max_scale + min_scale - effective_scale)
    else:
        geb_profile = np.zeros(ref_length, dtype=float)
    return gap_open_profile, gap_extend_profile, mismatch_profile, geb_profile


def build_homology_penalty_profile(
    ref_seq: str,
    homology_window: int = 8,
    homology_penalty: float = 0.0,
) -> np.ndarray:
    """构建同源区域惩罚profile — 对参考序列中存在多拷贝的区域施加match_score打折。

    homology_window: 检测同源性的滑动窗口大小（bp）。
    homology_penalty: ≤0，同源位置match_score的额外惩罚值，0=关闭。

    工作原理：
      对参考序列的每个位置i，提取以其为中心的homology_window长度的子序列。
      如果该子序列在参考序列的其他位置也存在，则位置i被标记为同源区域，
      在此处进行match时额外施加homology_penalty惩罚，使DP更不倾向在此处匹配。

    返回: 长度为len(ref_seq)的数组，每个元素≤0，homology_penalty=0时全0。
    """
    m = len(ref_seq)
    profile = np.zeros(m, dtype=float)
    if homology_penalty >= 0:
        return profile
    half_w = homology_window // 2
    for i in range(m):
        start = max(0, i - half_w)
        end = min(m, i + half_w + 1)
        fragment = ref_seq[start:end]
        # Search for this fragment elsewhere in the ref
        count = 0
        p = 0
        while True:
            p = ref_seq.find(fragment, p)
            if p < 0:
                break
            # Allow partial overlap at ends but not exact same position
            if abs(p - start) >= min(3, homology_window // 2):
                count += 1
                if count > 0:
                    break  # one non-self occurrence is enough
            p += 1
        if count > 0:
            profile[i] = homology_penalty
    return profile


def lineage_tracer_align(
    ref_seq: str,
    query_seq: str,
    cutsites: List[CutsiteRegion],
    match_score: float = 2.0,
    mismatch_penalty: float = -3.0,
    base_gap_open: float = -2.0,
    base_gap_extend: float = -0.1,
    min_scale: float = 1.0,
    max_scale: float = 6.0,
    cutsite_edge_scale: float = 2.0,
    gradient_radius: Optional[float] = None,
    mismatch_density_threshold: float = 0.34,
    sub_window: int = 3,
    no_gap_prefix: int = 0,
    gap_exit_bonus: float = 0.0,
    base_gap_exit: float = 0.0,
    short_match_window: int = 0,
    short_match_discount: float = 1.0,
    dense_mismatch_window: int = 6,
    dense_mismatch_penalty: float = 0.0,
    homology_window: int = 8,
    homology_penalty: float = 0.0,
    isolated_base_penalty: float = 0.0,
) -> Tuple[float, str, str, Dict]:
    gap_open_profile, gap_extend_profile, mismatch_profile, geb_profile = build_gradient_profiles(
        ref_length=len(ref_seq), cutsites=cutsites,
        base_gap_open=base_gap_open, base_gap_extend=base_gap_extend,
        mismatch_penalty=mismatch_penalty,
        min_scale=min_scale, max_scale=max_scale,
        cutsite_edge_scale=cutsite_edge_scale,
        gradient_radius=gradient_radius,
        base_gap_exit=base_gap_exit)
    if no_gap_prefix > 0:
        gap_open_profile[:no_gap_prefix] = -1e6
        gap_extend_profile[:no_gap_prefix] = -1e6
    homology_profile = build_homology_penalty_profile(
        ref_seq, homology_window=homology_window, homology_penalty=homology_penalty)
    # Use geb_profile when base_gap_exit is set; otherwise fall back to scalar
    geb_param = {} if base_gap_exit == 0.0 else {"gap_exit_bonus_profile": geb_profile}
    score, aligned_ref, aligned_query, stats = affine_gap_alignment_position_aware(
        ref_seq, query_seq, gap_open_profile, gap_extend_profile,
        match_score=match_score, mismatch_penalty=mismatch_penalty,
        mismatch_penalty_profile=mismatch_profile,
        gap_exit_bonus=gap_exit_bonus,
        **geb_param,
        short_match_window=short_match_window,
        short_match_discount=short_match_discount,
        dense_mismatch_window=dense_mismatch_window,
        dense_mismatch_threshold=mismatch_density_threshold,
        dense_mismatch_penalty=dense_mismatch_penalty,
        homology_profile=homology_profile,
        isolated_base_penalty=isolated_base_penalty)
    final_ref, final_query = aligned_ref, aligned_query
    final_stats = calculate_alignment_stats(final_ref, final_query)
    corrected_score = (final_stats['matches'] * 2 +
                       final_stats['mismatches'] * (-3) +
                       final_stats['gaps_in_ref'] * (-2) +
                       final_stats['gaps_in_query'] * (-2))
    final_stats['score'] = corrected_score
    return corrected_score, final_ref, final_query, final_stats


def get_amplicon_structure(ref_seq: str, config: AmpliconConfig = None) -> List[CutsiteRegion]:
    if config is None:
        config = AmpliconConfig.carlin_standard()
    if abs(len(ref_seq) - config.expected_full_length) <= 1:
        target_start = config.primer5_len + len(config.prefix)
        cutsites = []
        for i in range(config.n_targets):
            t_start = target_start + i * config.period
            cs_start = t_start + config.cutsite_offset
            cs_end = cs_start + config.cutsite_len - 1
            cutsites.append(CutsiteRegion(name=f"Target{i+1}", start=cs_start, end=cs_end))
        return cutsites
    motif = 'GAGTCG'
    positions = [m.start() for m in re.finditer(motif, ref_seq)]
    if len(positions) < 3:
        log.warning("无法在参考序列中检测到足够的GAGTCG motif（仅发现%d处）", len(positions))
        log.warning("请手动提供cutsite位置")
        return []
    diffs = [positions[i+1] - positions[i] for i in range(len(positions)-1)]
    first_pos = positions[0]
    period = config.period
    start = first_pos % period
    while start + period <= first_pos:
        start += period
    while start > first_pos:
        start -= period
    if start < 0:
        start = 0
    cutsites = []
    for i in range(config.n_targets):
        cs_start = start + i * period
        cs_end = cs_start + config.cutsite_len - 1
        if cs_start >= len(ref_seq):
            break
        if cs_end >= len(ref_seq):
            cs_end = len(ref_seq) - 1
        cutsites.append(CutsiteRegion(name=f"Target{i+1}", start=cs_start, end=cs_end))
    if cutsites:
        log.info("auto检测: 发现 %d 个cutsite区域 (周期=%dbp, motif起始=%d)", len(cutsites), period, start)
    return cutsites
