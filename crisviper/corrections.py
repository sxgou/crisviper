"""Post-alignment correction functions for repetitive misalignment and artifact removal."""

from typing import Tuple, List, Optional
from crisviper.config import CutsiteRegion, AmpliconConfig


def _adjacent_has_bases(aligned_query: str, col: int, window: int = 3) -> bool:
    """Check if position `col` in the aligned query has non-gap bases within `window` bp on either side.

    This prevents moving bases into a completely deleted (fully gapped) target region
    that has no adjacent non-gap bases — i.e., the gap spans the entire target region.
    """
    alen = len(aligned_query)
    # Check left side
    left_start = max(0, col - window)
    for k in range(left_start, col):
        if aligned_query[k] != '-':
            return True
    # Check right side
    right_end = min(alen, col + window + 1)
    for k in range(col + 1, right_end):
        if aligned_query[k] != '-':
            return True
    return False


def _find_repeat_occurrences(ref_seq: str, repeat: str) -> List[int]:
    """查找 repeat motif 在 ref_seq 中的所有出现位置（允许重叠）。"""
    positions = []
    p = 0
    while True:
        p = ref_seq.find(repeat, p)
        if p < 0:
            break
        positions.append(p)
        p += 1
    return positions


def _generate_repeat_pairs(
    positions: List[int],
    strategy: str = "forward",
) -> List[Tuple[int, int]]:
    """根据出现位置生成 (correct_rp, wrong_rp) 矫正对。

    "forward" — 靠前的 occurrence 为"correct"，靠后的为"wrong"
                [a, b, c] → [(a,b), (a,c), (b,c)]
    "reverse" — 靠后的 occurrence 为"correct"，靠前的为"wrong"
                [a, b, c] → [(c,a), (c,b), (b,a)]
    """
    if len(positions) < 2:
        return []
    if strategy == "reverse":
        return [
            (positions[j], positions[i])
            for i in range(len(positions))
            for j in range(i + 1, len(positions))
        ]
    else:  # forward
        return [
            (positions[i], positions[j])
            for i in range(len(positions))
            for j in range(i + 1, len(positions))
        ]


# 动态模式下的重复序列矫正策略
# "forward": 靠前者正确, 靠后者错误（默认）
# "reverse": 靠后者正确, 靠前者错误（如 ACTA）
_DYNAMIC_REPEAT_STRATEGIES = {
    "ACTGCACGACAGTCG": "forward",
    "ACTCGCG": "forward",
    "GAGCGC": "forward",
    "GCGACT": "forward",
    "GATACG": "forward",
    "ACGCAC": "forward",
    "CGCGCA": "forward",
    "CGACTA": "forward",
    "ACAGTCG": "forward",
    "GACGA": "forward",
    "ACTA": "reverse",
}

# 原始硬编码配置（向后兼容 fallback）
_HARDCODED_REPEATS_CONFIG: List[Tuple[str, Optional[List[Tuple[int, int]]]]] = [
    ("ACTGCACGACAGTCG", None),
    ("ACTCGCG", None),
    ("GAGCGC", None),
    ("GCGACT", None),
    ("GATACG", None),
    ("ACGCAC", None),
    ("CGCGCA", None),
    ("CGACTA", [(109, 258)]),
    ("ACAGTCG", [(37, 86), (37, 253), (86, 253)]),
    ("GACGA", None),
    ("ACTA", [(260, 83), (260, 125), (125, 83)]),
]


def convert_dense_mismatch_to_indel(
    aligned_ref: str,
    aligned_query: str,
    threshold: float = 0.34,
    min_region: int = 3
) -> Tuple[str, str, bool]:
    alen = len(aligned_ref)
    if alen == 0:
        return aligned_ref, aligned_query, False
    is_mismatch = [
        (ar != '-' and aq != '-' and ar != aq)
        for ar, aq in zip(aligned_ref, aligned_query)
    ]
    window_size = max(6, min_region)
    half = window_size // 2
    in_dense = False
    dense_regions = []
    region_start = 0
    for i in range(alen):
        start = max(0, i - half)
        end = min(alen, i + half + 1)
        n_bases = sum(1 for k in range(start, end)
                      if aligned_ref[k] != '-' and aligned_query[k] != '-')
        n_mismatch = sum(is_mismatch[start:end])
        density = n_mismatch / n_bases if n_bases > 0 else 0
        is_dense = density > threshold and n_mismatch >= 2
        if is_dense:
            if not in_dense:
                in_dense = True
                region_start = i
        else:
            if in_dense:
                if i - region_start >= min_region:
                    dense_regions.append((region_start, i - 1))
                in_dense = False
    if in_dense:
        if alen - region_start >= min_region:
            dense_regions.append((region_start, alen - 1))
    if not dense_regions:
        return aligned_ref, aligned_query, False
    ref_list = list(aligned_ref)
    was_modified = False
    for start, end in dense_regions:
        for k in range(start, end + 1):
            if ref_list[k] != '-':
                ref_list[k] = '-'
                was_modified = True
    return ''.join(ref_list), aligned_query, was_modified


def filter_point_mutations(
    aligned_ref: str,
    aligned_query: str,
    cutsites: List[CutsiteRegion],
    window: int = 3
) -> Tuple[str, str, int]:
    if len(aligned_ref) != len(aligned_query):
        return aligned_ref, aligned_query, 0
    result_query = list(aligned_query)
    n_corrected = 0
    ref_pos = 0
    for i, (ar, aq) in enumerate(zip(aligned_ref, aligned_query)):
        if ar == '-' or aq == '-':
            if ar != '-':
                ref_pos += 1
            continue
        if ar != aq:
            in_window = any(cs.start - window <= ref_pos <= cs.end + window for cs in cutsites)
            if not in_window:
                adj_to_gap = False
                if i > 0 and (aligned_ref[i - 1] == '-' or aligned_query[i - 1] == '-'):
                    adj_to_gap = True
                if i < len(aligned_ref) - 1 and (aligned_ref[i + 1] == '-' or aligned_query[i + 1] == '-'):
                    adj_to_gap = True
                if not adj_to_gap:
                    result_query[i] = ar
                    n_corrected += 1
        if ar != '-':
            ref_pos += 1
    return aligned_ref, ''.join(result_query), n_corrected


def correct_repetitive_misalignment(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: Optional[str] = None,
    repeats: Optional[List[Tuple[str, Optional[List[Tuple[int, int]]]]]] = None,
) -> Tuple[str, str, bool]:
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return aligned_ref, aligned_query, False
    ar, aq = aligned_ref, aligned_query
    alen = len(ar)
    pos_map = []
    cur = 0
    for c in ar:
        pos_map.append(cur if c != '-' else -1)
        if c != '-':
            cur += 1

    # repeats 参数 > ref_seq 动态检测 > 硬编码 fallback
    if repeats is not None:
        repeats_config = repeats
    elif ref_seq is not None:
        repeats_config = []
        for repeat, strategy in _DYNAMIC_REPEAT_STRATEGIES.items():
            rlen = len(repeat)
            if rlen > alen:
                continue
            positions = _find_repeat_occurrences(ref_seq, repeat)
            pairs = _generate_repeat_pairs(positions, strategy)
            if pairs:
                repeats_config.append((repeat, pairs))
    else:
        repeats_config = _HARDCODED_REPEATS_CONFIG

    new_aq = list(aq)
    overall_modified = False
    for repeat, specific_pairs in repeats_config:
        rlen = len(repeat)
        if rlen > alen:
            continue
        pairs = []
        if specific_pairs is not None:
            pairs.extend(specific_pairs)
        elif ref_seq is not None:
            copies = []
            p = 0
            while True:
                p = ref_seq.find(repeat, p)
                if p < 0:
                    break
                copies.append(p)
                p += 1
            if len(copies) < 2:
                continue
            for wp in copies[1:]:
                pairs.append((copies[0], wp))
        for correct_rp, wrong_rp in pairs:
            wrong_cols = [next((ci for ci, p in enumerate(pos_map) if p == wrong_rp + k), -1) for k in range(rlen)]
            correct_cols = [next((ci for ci, p in enumerate(pos_map) if p == correct_rp + k), -1) for k in range(rlen)]
            if any(c < 0 for c in wrong_cols) or any(c < 0 for c in correct_cols):
                continue
            if not all(ar[wrong_cols[k]] == repeat[k] for k in range(rlen)):
                continue
            match_cnt = sum(1 for k in range(rlen) if new_aq[wrong_cols[k]] == repeat[k] and new_aq[wrong_cols[k]] != '-')
            if match_cnt < rlen * 0.5:
                continue
            gaps_at_correct = sum(1 for k in range(rlen) if new_aq[correct_cols[k]] == '-')
            if gaps_at_correct < rlen * 0.5:
                continue
            # Skip correction if gap block has no adjacent bases (completely deleted region)
            if not _adjacent_has_bases(''.join(new_aq), correct_cols[0]) and \
               not _adjacent_has_bases(''.join(new_aq), correct_cols[-1]):
                continue
            modified = False
            for k in range(rlen):
                if new_aq[wrong_cols[k]] != '-' and new_aq[correct_cols[k]] == '-':
                    new_aq[correct_cols[k]] = new_aq[wrong_cols[k]]
                    new_aq[wrong_cols[k]] = '-'
                    modified = True
            if modified and repeat == "ACTGCACGACAGTCG":
                wrong_start_col = wrong_cols[0]
                g_src = None
                for sc in range(wrong_start_col - 1, max(0, wrong_start_col - 30) - 1, -1):
                    if aq[sc] != '-' and aq[sc] == 'G':
                        if ar[sc] == 'G' or ar[sc] == '-':
                            if pos_map[sc] < 0 or abs(pos_map[sc] - correct_rp) >= rlen:
                                if g_src is None:
                                    g_src = sc
                if g_src is not None and correct_rp > 0 and ref_seq is not None and ref_seq[correct_rp - 1] == 'G':
                    if new_aq[correct_cols[0] - 1] == '-':
                        new_aq[correct_cols[0] - 1] = aq[g_src]
                        new_aq[g_src] = '-'
                        modified = True
            if modified:
                overall_modified = True
    if overall_modified:
        return ar, ''.join(new_aq), True
    return aligned_ref, aligned_query, False


def correct_target_misalignments(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: str,
    cutsite_region_start: Optional[int] = None,
    cutsite_region_end: Optional[int] = None,
    shift_region_size: int = 90,
) -> Tuple[str, str, bool]:
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return aligned_ref, aligned_query, False
    ar, aq = aligned_ref, aligned_query
    alen = len(ar)
    pos_map = []
    cur = 0
    for c in ar:
        pos_map.append(cur if c != '-' else -1)
        if c != '-':
            cur += 1
    new_aq = list(aq)
    modified = False
    tag_seq = 'TAGTAT'
    t_len = 6
    t8_tag_start = ref_seq.find(tag_seq)
    if t8_tag_start >= 0:
        t8_cs = next((ci for ci in range(alen) if pos_map[ci] == t8_tag_start), None)
        t8_ce = next((ci for ci in range(alen) if pos_map[ci] == t8_tag_start + t_len - 1), None)
        if t8_cs is not None and t8_ce is not None:
            gaps = sum(1 for k in range(t8_cs, t8_ce + 1) if new_aq[k] == '-')
            if gaps >= t_len - 2:
                # 动态模式(cutsite_region_start提供): 从tag结束位置向后搜索; 向后兼容: 使用硬编码的236
                search_anchor = (t8_tag_start + t_len) if cutsite_region_start is not None else 236
                t8_end_col = next((ci for ci in range(alen) if pos_map[ci] == search_anchor), None)
                if t8_end_col is not None:
                    srch_start = t8_end_col + 1
                    srch_end = min(alen, srch_start + 50)
                    qb = [(ci, new_aq[ci]) for ci in range(srch_start, srch_end) if new_aq[ci] != '-']
                    cols = None
                    for wi in range(len(qb) - t_len + 1):
                        mc = sum(1 for k in range(t_len) if qb[wi + k][1] == tag_seq[k])
                        if mc >= t_len - 1:
                            cols = [qb[wi + k][0] for k in range(t_len)]
                            break
                    if cols:
                        for k in range(t_len):
                            dst = t8_cs + k
                            src = cols[k]
                            if new_aq[src] == tag_seq[k]:
                                new_aq[dst] = new_aq[src]
                                new_aq[src] = '-'
                                modified = True
    # 动态模式: 从传入参数推导坐标; 向后兼容: 使用原始硬编码
    if cutsite_region_start is not None and cutsite_region_end is not None:
        pos_a = cutsite_region_start
        pos_b = cutsite_region_end
    else:
        pos_a = 44
        pos_b = 260

    a_src_col = next((ci for ci in range(alen) if pos_map[ci] == pos_b), None)
    a_dst_col = next((ci for ci in range(alen) if pos_map[ci] == pos_a), None)
    if a_dst_col is not None and a_src_col is not None and new_aq[a_dst_col] == '-':
        srch_start = max(0, a_src_col - 3)
        srch_end = min(alen, a_src_col + 4)
        for ci in range(srch_start, srch_end):
            if ar[ci] == '-' and new_aq[ci] == 'A':
                new_aq[a_dst_col] = 'A'
                new_aq[ci] = '-'
                modified = True
                break
            elif ar[ci] == 'A' and new_aq[ci] == 'A':
                lg = ci > 0 and (ar[ci-1] == '-' or new_aq[ci-1] == '-')
                rg = ci < alen - 1 and (ar[ci+1] == '-' or new_aq[ci+1] == '-')
                if lg and rg:
                    new_aq[a_dst_col] = 'A'
                    new_aq[ci] = '-'
                    modified = True
                    break
    if modified:
        return ar, ''.join(new_aq), True
    return aligned_ref, aligned_query, False


def remove_isolated_matches(
    aligned_ref: str,
    aligned_query: str,
) -> Tuple[str, str, bool]:
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return aligned_ref, aligned_query, False
    ar, aq = aligned_ref, aligned_query
    alen = len(ar)
    new_aq = list(aq)
    modified = False
    for i in range(1, alen - 1):
        if ar[i] != '-' and new_aq[i] != '-' and ar[i] == new_aq[i]:
            left_gap = ar[i - 1] == '-' or new_aq[i - 1] == '-'
            right_gap = ar[i + 1] == '-' or new_aq[i + 1] == '-'
            if left_gap and right_gap:
                new_aq[i] = '-'
                modified = True
    if modified:
        return ar, ''.join(new_aq), modified
    return aligned_ref, aligned_query, False
