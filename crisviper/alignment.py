"""Core sequence alignment algorithms — standard Gotoh and position-aware DP."""

import numpy as np
from typing import Tuple, Dict, List, Optional
from numba import jit


def count_gap_blocks(seq: str) -> List[int]:
    blocks = []
    in_gap = False
    current_length = 0
    for c in seq:
        if c == '-':
            if not in_gap:
                in_gap = True
                current_length = 1
            else:
                current_length += 1
        else:
            if in_gap:
                blocks.append(current_length)
                in_gap = False
    if in_gap:
        blocks.append(current_length)
    return blocks


def calculate_alignment_stats(aligned_ref: str, aligned_query: str) -> Dict:
    matches = sum(1 for a, b in zip(aligned_ref, aligned_query) if a == b and a != '-')
    mismatches = sum(1 for a, b in zip(aligned_ref, aligned_query)
                     if a != b and a != '-' and b != '-')
    gaps_in_ref = aligned_ref.count('-')
    gaps_in_query = aligned_query.count('-')
    gap_blocks_ref = count_gap_blocks(aligned_ref)
    gap_blocks_query = count_gap_blocks(aligned_query)
    avg_gap_len_ref = np.mean(gap_blocks_ref) if gap_blocks_ref else 0
    avg_gap_len_query = np.mean(gap_blocks_query) if gap_blocks_query else 0
    alignment_length = len(aligned_ref)
    similarity = matches / alignment_length if alignment_length > 0 else 0
    return {
        'matches': matches,
        'mismatches': mismatches,
        'gaps_in_ref': gaps_in_ref,
        'gaps_in_query': gaps_in_query,
        'gap_blocks_ref': gap_blocks_ref,
        'gap_blocks_query': gap_blocks_query,
        'avg_gap_len_ref': avg_gap_len_ref,
        'avg_gap_len_query': avg_gap_len_query,
        'alignment_length': alignment_length,
        'similarity': similarity,
        'identity': matches / (matches + mismatches) if (matches + mismatches) > 0 else 0,
    }


def _build_score_matrix(ref_seq: str, query_seq: str,
                        match_score: float, mismatch_penalty: float,
                        mismatch_penalty_profile: Optional[np.ndarray],
                        short_match_window: int, short_match_discount: float,
                        homology_profile: Optional[np.ndarray],
                        gap_exit_bonus: float,
                        gap_exit_bonus_profile: Optional[np.ndarray],
                        dense_mismatch_active: bool,
                        dense_mismatch_window: int,
                        isolated_base_penalty: float,
                        ) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray]]:
    """Pre-compute the m×n substitution score matrix and auxiliary arrays.

    Returns: (sub_score, run_len, density, is_match)
    """
    m, n = len(ref_seq), len(query_seq)

    # Vectorized match detection
    ref_codes = np.frombuffer(ref_seq.encode(), dtype=np.int8).copy()
    qry_codes = np.frombuffer(query_seq.encode(), dtype=np.int8).copy()
    is_match = ref_codes[:, None] == qry_codes[None, :]  # m×n

    # Base substitution score
    if mismatch_penalty_profile is not None:
        mpp = mismatch_penalty_profile[:m] if len(mismatch_penalty_profile) >= m else mismatch_penalty_profile
        sub_score = np.where(is_match, match_score, mpp[:, None])
    else:
        sub_score = np.where(is_match, float(match_score), float(mismatch_penalty))

    # Homology penalty: subtract from match positions
    if homology_profile is not None:
        hp = homology_profile[:m] if len(homology_profile) >= m else homology_profile
        if (hp < 0).any():
            sub_score = sub_score + hp[:, None] * is_match

    # Short-match discount on run lengths
    short_match_active = (short_match_window > 0 and short_match_discount < 1.0)
    needs_run_len = short_match_active or dense_mismatch_active or gap_exit_bonus != 0.0 or isolated_base_penalty < 0.0
    run_len = None
    if needs_run_len:
        run_len = np.zeros((m + 1, n + 1), dtype=np.int32)
        _compute_run_len(is_match, run_len, m, n)

        if short_match_active:
            mask = (run_len[:m, :n] > 0) & (run_len[:m, :n] <= short_match_window)
            sub_score[mask] = match_score * short_match_discount

        # Isolated base penalty: run_len == 1 at this cell means it's an isolated match
        if isolated_base_penalty < 0.0:
            iso_mask = run_len[:m, :n] == 1
            sub_score[iso_mask] += isolated_base_penalty

    # Dense mismatch density
    density = None
    if dense_mismatch_active:
        density = _compute_dense_mismatch_density(is_match, m, n, dense_mismatch_window)

    return sub_score, run_len, density, is_match


@jit(nopython=True, cache=True)
def _compute_run_len_numba(is_match, run_len, m, n):
    """Numba-JIT: forward-looking match run length."""
    for i in range(m - 1, -1, -1):
        for j in range(n - 1, -1, -1):
            if is_match[i, j]:
                run_len[i, j] = 1 + run_len[i + 1, j + 1]


def _compute_run_len(is_match: np.ndarray, run_len: np.ndarray, m: int, n: int) -> None:
    _compute_run_len_numba(is_match, run_len, m, n)


@jit(nopython=True, cache=True)
def _dense_mm_density_numba(is_match, m, n, window):
    """Numba-JIT: dense mismatch density via sliding window on each diagonal."""
    half_w = window // 2
    density = np.zeros((m + 1, n + 1), dtype=np.float64)

    for d in range(-(n - 1), m):
        if d >= 0:
            i0, j0 = d, 0
            length = min(m - d, n)
        else:
            i0, j0 = 0, -d
            length = min(m, n + d)
        if length < 1:
            continue

        # Initial window sum for k=0: positions [0, half_w] on diagonal
        win_sum = 0
        for t in range(min(half_w + 1, length)):
            if is_match[i0 + t, j0 + t]:
                win_sum += 0
            else:
                win_sum += 1

        for k in range(length):
            denom = min(k + half_w + 1, length - k + half_w, window)
            density[i0 + k, j0 + k] = win_sum / max(denom, 1)

            # Slide window: remove k-half_w, add k+half_w+1
            out_pos = k - half_w
            if out_pos >= 0 and not is_match[i0 + out_pos, j0 + out_pos]:
                win_sum -= 1
            in_pos = k + half_w + 1
            if in_pos < length and not is_match[i0 + in_pos, j0 + in_pos]:
                win_sum += 1

    return density


def _compute_dense_mismatch_density(is_match: np.ndarray, m: int, n: int,
                                    window: int) -> np.ndarray:
    return _dense_mm_density_numba(is_match, m, n, window)


@jit(nopython=True, cache=True)
def _dp_fill_numba(
    M, Ix, Iy,
    sub_score,
    gap_open_profile,
    gap_extend_profile,
    geb_profile,
    use_profile_geb,
    geb_scalar,
):
    """Numba-JIT compiled DP fill for affine gap alignment (Gotoh algorithm).

    Fills M, Ix, Iy matrices row by row. All arrays are modified in-place.
    """
    m = M.shape[0] - 1
    n = M.shape[1] - 1
    for i in range(1, m + 1):
        go_i = gap_open_profile[i - 1]
        ge_i = gap_extend_profile[i - 1]
        go_ge_i = go_i + ge_i

        # Gap exit bonus for this row
        geb_i = geb_profile[i - 1] if use_profile_geb else geb_scalar

        for j in range(1, n + 1):
            # Iy[i, j] = max(M[i-1,j] + go+ge, Ix[i-1,j] + go+ge, Iy[i-1,j] + ge)
            v1 = M[i - 1, j] + go_ge_i
            v2 = Ix[i - 1, j] + go_ge_i
            v3 = Iy[i - 1, j] + ge_i
            Iy[i, j] = v1 if v1 >= v2 and v1 >= v3 else (v2 if v2 >= v3 else v3)

            # M[i, j] = s(i-1,j-1) + max(M[i-1,j-1], Ix[i-1,j-1] + geb, Iy[i-1,j-1] + geb)
            v1 = M[i - 1, j - 1]
            v2 = Ix[i - 1, j - 1] + geb_i
            v3 = Iy[i - 1, j - 1] + geb_i
            best = v1 if v1 >= v2 and v1 >= v3 else (v2 if v2 >= v3 else v3)
            M[i, j] = sub_score[i - 1, j - 1] + best

            # Ix[i, j] = max(M[i,j-1] + go+ge, Ix[i,j-1] + ge, Iy[i,j-1] + go+ge)
            v1 = M[i, j - 1] + go_ge_i
            v2 = Ix[i, j - 1] + ge_i
            v3 = Iy[i, j - 1] + go_ge_i
            Ix[i, j] = v1 if v1 >= v2 and v1 >= v3 else (v2 if v2 >= v3 else v3)


@jit(nopython=True, cache=True)
def _dp_fill_numba_standard(M, Ix, Iy, ref_seq, query_seq,
                             match_score, mismatch_penalty,
                             gap_open, gap_extend, go_ge,
                             gap_exit_bonus, use_gap_exit,
                             run_len, use_short_match,
                             short_match_window, short_match_discount):
    """Numba-JIT: standard Gotoh DP fill (no profiles, no homology)."""
    m = len(ref_seq)
    n = len(query_seq)
    for i in range(1, m + 1):
        r_char = ref_seq[i - 1]
        for j in range(1, n + 1):
            q_char = query_seq[j - 1]

            # Substitution score
            if use_short_match and r_char == q_char:
                run = run_len[i - 1, j - 1]
                if 0 < run <= short_match_window:
                    s = match_score * short_match_discount
                else:
                    s = match_score
            else:
                s = match_score if r_char == q_char else mismatch_penalty

            # M[i, j]
            if use_gap_exit:
                best_prev = max(M[i - 1, j - 1],
                                Ix[i - 1, j - 1] + gap_exit_bonus,
                                Iy[i - 1, j - 1] + gap_exit_bonus)
            else:
                best_prev = max(M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1])
            M[i, j] = s + best_prev

            # Ix[i, j] = max(M[i,j-1] + go+ge, Ix[i,j-1] + ge, Iy[i,j-1] + go+ge)
            v1 = M[i, j - 1] + go_ge
            v2 = Ix[i, j - 1] + gap_extend
            v3 = Iy[i, j - 1] + go_ge
            Ix[i, j] = v1 if v1 >= v2 and v1 >= v3 else (v2 if v2 >= v3 else v3)

            # Iy[i, j] = max(M[i-1,j] + go+ge, Ix[i-1,j] + go+ge, Iy[i-1,j] + ge)
            v1 = M[i - 1, j] + go_ge
            v2 = Ix[i - 1, j] + go_ge
            v3 = Iy[i - 1, j] + gap_extend
            Iy[i, j] = v1 if v1 >= v2 and v1 >= v3 else (v2 if v2 >= v3 else v3)


def affine_gap_alignment(ref_seq: str, query_seq: str,
                         match_score: float = 2.0,
                         mismatch_penalty: float = -3.0,
                         gap_open: float = -2.0,
                         gap_extend: float = -0.1,
                         gap_exit_bonus: float = 0.0,
                         short_match_window: int = 0,
                         short_match_discount: float = 1.0) -> Tuple[float, str, str, Dict]:
    m, n = len(ref_seq), len(query_seq)
    M = np.zeros((m + 1, n + 1), dtype=float)
    Ix = np.zeros((m + 1, n + 1), dtype=float)
    Iy = np.zeros((m + 1, n + 1), dtype=float)
    M[0, 0] = 0.0
    Ix[0, 0] = gap_open
    Iy[0, 0] = gap_open
    for i in range(1, m + 1):
        M[i, 0] = -np.inf
        Ix[i, 0] = -np.inf
        Iy[i, 0] = gap_open + gap_extend * (i - 1)
    for j in range(1, n + 1):
        M[0, j] = -np.inf
        Ix[0, j] = gap_open + gap_extend * (j - 1)
        Iy[0, j] = -np.inf

    go_ge = gap_open + gap_extend
    use_gap_exit = (gap_exit_bonus != 0.0)
    use_short_match = (short_match_window > 0 and short_match_discount < 1.0)

    if use_short_match:
        is_match_arr = np.array([[ref_seq[i] == query_seq[j]
                                   for j in range(n)] for i in range(m)])
        run_len = np.zeros((m + 1, n + 1), dtype=np.int32)
        _compute_run_len_numba(is_match_arr, run_len, m, n)
    else:
        run_len = np.zeros((m + 1, n + 1), dtype=np.int32)

    _dp_fill_numba_standard(
        M, Ix, Iy, ref_seq, query_seq,
        match_score, mismatch_penalty,
        gap_open, gap_extend, go_ge,
        gap_exit_bonus, use_gap_exit,
        run_len, use_short_match,
        short_match_window, short_match_discount,
    )

    # Global alignment: only M[m, n] is valid at termination
    max_score = M[m, n]
    state = 'M'
    i, j = m, n
    aligned_ref, aligned_query = [], []
    while i > 0 and j > 0:
        if state == 'M':
            aligned_ref.append(ref_seq[i - 1])
            aligned_query.append(query_seq[j - 1])
            scores = [M[i - 1, j - 1],
                      Ix[i - 1, j - 1] + (gap_exit_bonus if use_gap_exit else 0),
                      Iy[i - 1, j - 1] + (gap_exit_bonus if use_gap_exit else 0)]
            state = ['M', 'Ix', 'Iy'][np.argmax(scores)]
            i -= 1; j -= 1
        elif state == 'Ix':
            aligned_ref.append('-')
            aligned_query.append(query_seq[j - 1])
            scores = [M[i, j - 1] + gap_open + gap_extend,
                      Ix[i, j - 1] + gap_extend,
                      Iy[i, j - 1] + gap_open + gap_extend]
            state = ['M', 'Ix', 'Iy'][np.argmax(scores)]
            j -= 1
        else:
            aligned_ref.append(ref_seq[i - 1])
            aligned_query.append('-')
            scores = [M[i - 1, j] + gap_open + gap_extend,
                      Ix[i - 1, j] + gap_open + gap_extend,
                      Iy[i - 1, j] + gap_extend]
            state = ['M', 'Ix', 'Iy'][np.argmax(scores)]
            i -= 1
    while i > 0:
        aligned_ref.append(ref_seq[i - 1]); aligned_query.append('-'); i -= 1
    while j > 0:
        aligned_ref.append('-'); aligned_query.append(query_seq[j - 1]); j -= 1
    aligned_ref = ''.join(reversed(aligned_ref))
    aligned_query = ''.join(reversed(aligned_query))
    stats = calculate_alignment_stats(aligned_ref, aligned_query)
    stats['score'] = max_score
    return max_score, aligned_ref, aligned_query, stats


def affine_gap_alignment_position_aware(
    ref_seq: str,
    query_seq: str,
    gap_open_profile: np.ndarray,
    gap_extend_profile: np.ndarray,
    match_score: float = 2.0,
    mismatch_penalty: float = -3.0,
    mismatch_penalty_profile: Optional[np.ndarray] = None,
    gap_exit_bonus: float = 0.0,
    gap_exit_bonus_profile: Optional[np.ndarray] = None,
    short_match_window: int = 0,
    short_match_discount: float = 1.0,
    dense_mismatch_window: int = 6,
    dense_mismatch_threshold: float = 0.34,
    dense_mismatch_penalty: float = 0.0,
    homology_profile: Optional[np.ndarray] = None,
    isolated_base_penalty: float = 0.0,
) -> Tuple[float, str, str, Dict]:
    m, n = len(ref_seq), len(query_seq)

    # ── Pre-compute score and auxiliary matrices (vectorized) ──
    dense_mismatch_active = (dense_mismatch_penalty < 0.0 and mismatch_penalty_profile is not None)
    sub_score, run_len, density, is_match = _build_score_matrix(
        ref_seq, query_seq, match_score, mismatch_penalty,
        mismatch_penalty_profile,
        short_match_window, short_match_discount,
        homology_profile, gap_exit_bonus, gap_exit_bonus_profile,
        dense_mismatch_active, dense_mismatch_window,
        isolated_base_penalty)

    # Apply dense mismatch penalty to substitution scores
    if dense_mismatch_active:
        dense_mask = density[:m, :n] >= dense_mismatch_threshold
        sub_score[dense_mask] += dense_mismatch_penalty
        # don't let it go below mismatch penalty
        if mismatch_penalty_profile is not None:
            mpp = mismatch_penalty_profile[:m] if len(mismatch_penalty_profile) >= m else mismatch_penalty_profile
            min_allowed = mpp[:, None]
            # but only clamp mismatch positions
            mismatch_positions = ~is_match
            sub_score = np.where(mismatch_positions & dense_mask,
                                 np.maximum(sub_score, min_allowed),
                                 sub_score)
        else:
            mismatch_positions = ~is_match
            sub_score = np.where(mismatch_positions & dense_mask,
                                 np.maximum(sub_score, mismatch_penalty),
                                 sub_score)

    # ── DP matrix initialization ──
    M = np.full((m + 1, n + 1), -np.inf, dtype=float)
    Ix = np.full((m + 1, n + 1), -np.inf, dtype=float)
    Iy = np.full((m + 1, n + 1), -np.inf, dtype=float)
    M[0, 0] = 0.0
    Ix[0, 0] = gap_open_profile[0]
    Iy[0, 0] = gap_open_profile[0]
    for i in range(1, m + 1):
        Iy[i, 0] = gap_open_profile[0] + gap_extend_profile[0] * (i - 1)
    for j in range(1, n + 1):
        Ix[0, j] = gap_open_profile[0] + gap_extend_profile[0] * (j - 1)

    # ── DP fill via numba JIT ──
    gap_exit_active = (gap_exit_bonus != 0.0 or gap_exit_bonus_profile is not None)
    use_profile_geb = gap_exit_bonus_profile is not None
    geb_profile = gap_exit_bonus_profile if gap_exit_bonus_profile is not None else np.zeros(m, dtype=np.float64)
    _dp_fill_numba(M, Ix, Iy, sub_score,
                   gap_open_profile.astype(np.float64),
                   gap_extend_profile.astype(np.float64),
                   geb_profile, use_profile_geb, gap_exit_bonus)

    # ── Backtrace ──
    # Global alignment: only M[m, n] is valid at termination
    max_score = M[m, n]
    state = 'M'
    i, j = m, n
    aligned_ref, aligned_query = [], []
    while i > 0 and j > 0:
        if state == 'M':
            aligned_ref.append(ref_seq[i - 1])
            aligned_query.append(query_seq[j - 1])
            if gap_exit_active:
                if use_profile_geb:
                    geb_i = gap_exit_bonus_profile[i - 1]
                else:
                    geb_i = gap_exit_bonus
                scores = [M[i - 1, j - 1], Ix[i - 1, j - 1] + geb_i, Iy[i - 1, j - 1] + geb_i]
            else:
                scores = [M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1]]
            state = ['M', 'Ix', 'Iy'][np.argmax(scores)]
            i -= 1; j -= 1
        elif state == 'Ix':
            aligned_ref.append('-')
            aligned_query.append(query_seq[j - 1])
            ri = max(0, i - 1)
            scores = [M[i, j - 1] + gap_open_profile[ri] + gap_extend_profile[ri],
                      Ix[i, j - 1] + gap_extend_profile[ri],
                      Iy[i, j - 1] + gap_open_profile[ri] + gap_extend_profile[ri]]
            state = ['M', 'Ix', 'Iy'][np.argmax(scores)]
            j -= 1
        else:
            aligned_ref.append(ref_seq[i - 1])
            aligned_query.append('-')
            ri = i - 1
            scores = [M[i - 1, j] + gap_open_profile[ri] + gap_extend_profile[ri],
                      Ix[i - 1, j] + gap_open_profile[ri] + gap_extend_profile[ri],
                      Iy[i - 1, j] + gap_extend_profile[ri]]
            state = ['M', 'Ix', 'Iy'][np.argmax(scores)]
            i -= 1
    while i > 0:
        aligned_ref.append(ref_seq[i - 1]); aligned_query.append('-'); i -= 1
    while j > 0:
        aligned_ref.append('-'); aligned_query.append(query_seq[j - 1]); j -= 1
    aligned_ref = ''.join(reversed(aligned_ref))
    aligned_query = ''.join(reversed(aligned_query))
    stats = calculate_alignment_stats(aligned_ref, aligned_query)
    stats['score'] = max_score
    return max_score, aligned_ref, aligned_query, stats
