"""Core sequence alignment algorithms — standard Gotoh and position-aware DP."""

import numpy as np
from typing import Tuple, Dict, List, Optional


def count_gap_blocks(seq: str) -> List[int]:
    """
    统计序列中的连续gap区块长度
    
    参数：
        seq: 比对序列（可能包含gap字符'-'）
    
    返回：
        连续gap区块长度列表
    """
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
    """
    计算比对统计信息
    
    参数：
        aligned_ref: 比对后的reference序列
        aligned_query: 比对后的query序列
    
    返回：
        统计信息字典
    """
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
    
    stats = {
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
        'identity': matches / (matches + mismatches) if (matches + mismatches) > 0 else 0
    }
    
    return stats


def affine_gap_alignment(ref_seq: str, query_seq: str,
                         match_score: float = 2.0,
                         mismatch_penalty: float = -3.0,
                         gap_open: float = -2.0,
                         gap_extend: float = -0.1) -> Tuple[float, str, str, Dict]:
    """
    仿射gap惩罚的全局序列比对（两端严格对齐）
    """
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
    for i in range(1, m + 1):
        Mi_1 = M[i - 1]; Mi = M[i]
        Ixi_1 = Ix[i - 1]; Ixi = Ix[i]
        Iyi_1 = Iy[i - 1]; Iyi = Iy[i]
        r_char = ref_seq[i - 1]
        for j in range(1, n + 1):
            s = match_score if r_char == query_seq[j - 1] else mismatch_penalty
            a, b, c = Mi_1[j - 1], Ixi_1[j - 1], Iyi_1[j - 1]
            Mi[j] = s + (a if a >= b and a >= c else (b if b >= c else c))
            a, b, c = Mi[j - 1] + go_ge, Ixi[j - 1] + gap_extend, Iyi[j - 1] + go_ge
            Ixi[j] = a if a >= b and a >= c else (b if b >= c else c)
            a, b, c = Mi_1[j] + go_ge, Ixi_1[j] + go_ge, Iyi_1[j] + gap_extend
            Iyi[j] = a if a >= b and a >= c else (b if b >= c else c)

    max_score = max(M[m, n], Ix[m, n], Iy[m, n])
    max_i, max_j = m, n
    max_state = 'M' if M[m, n] >= max(Ix[m, n], Iy[m, n]) else \
               ('Ix' if Ix[m, n] >= Iy[m, n] else 'Iy')

    i, j = max_i, max_j
    state = max_state
    aligned_ref, aligned_query = [], []

    while i > 0 and j > 0:
        if state == 'M':
            if i > 0 and j > 0:
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append(query_seq[j - 1])
                scores = [M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1]]
                prev = np.argmax(scores)
                state = ['M', 'Ix', 'Iy'][prev]
                i -= 1
                j -= 1
            elif i > 0:
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append('-')
                i -= 1
                state = 'Iy'
            else:
                break
        elif state == 'Ix':
            aligned_ref.append('-')
            aligned_query.append(query_seq[j - 1])
            scores = [
                M[i, j - 1] + gap_open + gap_extend,
                Ix[i, j - 1] + gap_extend,
                Iy[i, j - 1] + gap_open + gap_extend
            ]
            prev = np.argmax(scores)
            state = ['M', 'Ix', 'Iy'][prev]
            j -= 1
        else:
            aligned_ref.append(ref_seq[i - 1])
            aligned_query.append('-')
            scores = [
                M[i - 1, j] + gap_open + gap_extend,
                Ix[i - 1, j] + gap_open + gap_extend,
                Iy[i - 1, j] + gap_extend
            ]
            prev = np.argmax(scores)
            state = ['M', 'Ix', 'Iy'][prev]
            i -= 1

    while i > 0:
        aligned_ref.append(ref_seq[i - 1])
        aligned_query.append('-')
        i -= 1
    while j > 0:
        aligned_ref.append('-')
        aligned_query.append(query_seq[j - 1])
        j -= 1

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
) -> Tuple[float, str, str, Dict]:
    """
    位置感知的仿射gap全局序列比对（两端严格对齐）。

    使用位置依赖的gap惩罚数组，允许在cutsite区域优先开启gap。
    gap_open_profile[i] 和 gap_extend_profile[i] 对应参考序列第i位的惩罚。

    返回: (score, aligned_ref, aligned_query, stats)
    """
    m, n = len(ref_seq), len(query_seq)

    M = np.full((m + 1, n + 1), -np.inf, dtype=float)
    Ix = np.full((m + 1, n + 1), -np.inf, dtype=float)
    Iy = np.full((m + 1, n + 1), -np.inf, dtype=float)

    M[0, 0] = 0.0
    Ix[0, 0] = gap_open_profile[0]
    Iy[0, 0] = gap_open_profile[0]
    for i in range(1, m + 1):
        M[i, 0] = -np.inf
        Ix[i, 0] = -np.inf
        Iy[i, 0] = gap_open_profile[0] + gap_extend_profile[0] * (i - 1)
    for j in range(1, n + 1):
        M[0, j] = -np.inf
        Ix[0, j] = gap_open_profile[0] + gap_extend_profile[0] * (j - 1)
        Iy[0, j] = -np.inf

    for i in range(1, m + 1):
        Mi_1 = M[i - 1]; Mi = M[i]
        Ixi_1 = Ix[i - 1]; Ixi = Ix[i]
        Iyi_1 = Iy[i - 1]; Iyi = Iy[i]
        r_char = ref_seq[i - 1]

        gp_open_i_1 = gap_open_profile[i - 1]
        gp_extend_i_1 = gap_extend_profile[i - 1]
        gp_open_i = gap_open_profile[i - 1]
        gp_extend_i = gap_extend_profile[i - 1]
        go_ge_i = gp_open_i + gp_extend_i
        go_ge_i_1 = gp_open_i_1 + gp_extend_i_1

        for j in range(1, n + 1):
            s = match_score if r_char == query_seq[j - 1] else (
                mismatch_penalty_profile[i - 1] if mismatch_penalty_profile is not None else mismatch_penalty)
            a, b, c = Mi_1[j - 1], Ixi_1[j - 1], Iyi_1[j - 1]
            Mi[j] = s + (a if a >= b and a >= c else (b if b >= c else c))
            a, b, c = Mi[j - 1] + go_ge_i, Ixi[j - 1] + gp_extend_i, Iyi[j - 1] + go_ge_i
            Ixi[j] = a if a >= b and a >= c else (b if b >= c else c)
            a, b, c = Mi_1[j] + go_ge_i_1, Ixi_1[j] + go_ge_i_1, Iyi_1[j] + gp_extend_i_1
            Iyi[j] = a if a >= b and a >= c else (b if b >= c else c)

    max_score = max(M[m, n], Ix[m, n], Iy[m, n])
    max_i, max_j = m, n
    max_state = 'M' if M[m, n] >= max(Ix[m, n], Iy[m, n]) else \
                ('Ix' if Ix[m, n] >= Iy[m, n] else 'Iy')

    i, j = max_i, max_j
    state = max_state
    aligned_ref, aligned_query = [], []

    while i > 0 and j > 0:
        if state == 'M':
            if i > 0 and j > 0:
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append(query_seq[j - 1])
                scores = [M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1]]
                prev = np.argmax(scores)
                state = ['M', 'Ix', 'Iy'][prev]
                i -= 1
                j -= 1
            elif i > 0:
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append('-')
                i -= 1
                state = 'Iy'
            else:
                break
        elif state == 'Ix':
            if j > 0:
                aligned_ref.append('-')
                aligned_query.append(query_seq[j - 1])
                ri = max(0, i - 1)
                scores = [
                    M[i, j - 1] + gap_open_profile[ri] + gap_extend_profile[ri],
                    Ix[i, j - 1] + gap_extend_profile[ri],
                    Iy[i, j - 1] + gap_open_profile[ri] + gap_extend_profile[ri]
                ]
                prev = np.argmax(scores)
                state = ['M', 'Ix', 'Iy'][prev]
                j -= 1
            else:
                break
        else:
            if i > 0:
                aligned_ref.append(ref_seq[i - 1])
                aligned_query.append('-')
                ri = i - 1
                scores = [
                    M[i - 1, j] + gap_open_profile[ri] + gap_extend_profile[ri],
                    Ix[i - 1, j] + gap_open_profile[ri] + gap_extend_profile[ri],
                    Iy[i - 1, j] + gap_extend_profile[ri]
                ]
                prev = np.argmax(scores)
                state = ['M', 'Ix', 'Iy'][prev]
                i -= 1
            else:
                break

    while i > 0:
        aligned_ref.append(ref_seq[i - 1])
        aligned_query.append('-')
        i -= 1
    while j > 0:
        aligned_ref.append('-')
        aligned_query.append(query_seq[j - 1])
        j -= 1

    aligned_ref = ''.join(reversed(aligned_ref))
    aligned_query = ''.join(reversed(aligned_query))

    stats = calculate_alignment_stats(aligned_ref, aligned_query)
    stats['score'] = max_score

    return max_score, aligned_ref, aligned_query, stats
