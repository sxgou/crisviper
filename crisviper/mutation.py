"""crisviper/mutation.py — 突变识别独立模块

从比对结果中显式提取结构化突变事件。
这是管道的核心步骤——"突变识别"。
"""

from typing import List, Optional, Tuple
from crisviper.models import MutationEvent, MutationType, AlignmentStats
from crisviper.config import CutsiteRegion


# ═══════════════════════════════════════════════════════════════
# MATLAB兼容的突变识别（identify_sequence_events / identify_cas9_events）
#
# 这些函数实现与MATLAB @Mutation类相同的事件识别逻辑，
# 产生一致的事件结构，用于HGVS注释输出。
# ═══════════════════════════════════════════════════════════════

def _build_ref_pos_map_full(aligned_ref: str):
    """Build alignment-column → ref-position map (0-based ref coords).

    Returns (pos_map, total) where:
      - pos_map: list length=len(aligned_ref), -1 for gap columns
      - total: number of non-gap ref bases
    """
    pos_map = []
    cur = 0
    for c in aligned_ref:
        if c != '-':
            pos_map.append(cur)
            cur += 1
        else:
            pos_map.append(-1)
    return pos_map, cur


def _find_site_for_bp(bp: int, cutsites: List[CutsiteRegion]) -> Tuple[int, str]:
    """Find which motif site a base position belongs to.

    Returns (site_index, site_type) where site_type is 'cutsite' or 'consite'.
    MATLAB equivalent: CARLIN_def.locate(bp) → struct with .abs and .type

    Args:
        bp: 0-based reference position
        cutsites: list of CutsiteRegion defining cutsite boundaries
    """
    for i, cs in enumerate(cutsites):
        if cs.start <= bp <= cs.end:
            return (i, 'cutsite')
        # consite is the conserved part before the cutsite
        consite_end = cs.start - 1
        # Find which consite this bp belongs to
        if i == 0:
            if bp < cs.start:
                return (i, 'prefix')
        else:
            prev_cs_end = cutsites[i - 1].end
            if prev_cs_end < bp < cs.start:
                return (i, 'consite')
    # After last cutsite → postfix
    return (len(cutsites), 'postfix')


def identify_sequence_events(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: str,
    cutsites: Optional[List[CutsiteRegion]] = None,
    mutation_window: int = 3,
) -> List[MutationEvent]:
    """MATLAB-compatible identify_sequence_events — naive mutation extraction.

    MATLAB equivalent: @Mutation/identify_sequence_events

    Identifies mutations from aligned sequences using MATLAB's approach:
    1. Deletions (query gap, ref base)
    2. Insertions (ref gap, query base) — with left/right alignment by site
    3. Substitutions (mismatched bases)

    Insertion position assignment follows MATLAB convention:
    - After cutsite: position assigned to cutsite+1 (left-aligned to cutsite)
    - Elsewhere: position assigned to preceding base (right-aligned)

    Args:
        aligned_ref: Gapped reference sequence.
        aligned_query: Gapped query sequence.
        ref_seq: Full ungapped reference sequence.
        cutsites: Optional list of CutsiteRegion for insertion positioning.
        mutation_window: Window radius for cutsite proximity detection.

    Returns:
        List of MutationEvent objects in MATLAB order (D, I, M sorted by position).
    """
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return []

    N = len(ref_seq)
    alen = len(aligned_ref)
    ref_pos_map, _ = _build_ref_pos_map_full(aligned_ref)

    # Step 1: Build ref_mask = alignment column indices of non-gap ref bases
    # MATLAB equivalent: ref_mask = find(~ref_gap)
    ref_mask = [i for i, p in enumerate(ref_pos_map) if p >= 0]
    if not ref_mask:
        return []

    # bp_event tracks per-base-pair status: 'N' (normal), 'D' (deleted), 'I' (inserted)
    bp_event = ['N'] * N

    events: List[MutationEvent] = []

    # ── Step 2: Find deletions ──
    del_regions = []
    i = 0
    while i < alen:
        if aligned_ref[i] != '-' and aligned_query[i] == '-':
            start = i
            while i < alen and aligned_ref[i] != '-' and aligned_query[i] == '-':
                i += 1
            end = i - 1
            del_regions.append((start, end))
        else:
            i += 1

    for (del_start, del_end) in del_regions:
        bp_start = ref_pos_map[del_start]
        bp_end = ref_pos_map[del_end]
        length = bp_end - bp_start + 1
        seg_ref = aligned_ref[del_start:del_end + 1]

        for bp in range(bp_start, bp_end + 1):
            bp_event[bp] = 'D'

        event = MutationEvent(
            type=MutationType.DELETION,
            ref_pos=bp_start,
            ref_base=seg_ref,
            query_base='-' * length,
            length=length,
            in_cutsite_window=_in_any_cutsite_window(bp_start, cutsites, mutation_window),
            raw_ref_segment=seg_ref,
            raw_query_segment='-' * length,
        )
        events.append(event)

    # ── Step 3: Find insertions ──
    # Gap between ref_mask[i] and ref_mask[i+1] in the alignment = insertion
    for i in range(len(ref_mask) - 1):
        aln_i = ref_mask[i]       # alignment column of i-th ref base
        aln_j = ref_mask[i + 1]   # alignment column of (i+1)-th ref base
        gap_cols = aln_j - aln_i - 1  # number of gap columns between

        if gap_cols > 0:
            cur_bp = ref_pos_map[aln_i]  # 0-based ref position
            # Determine insertion position (MATLAB convention)
            is_after_cutsite = False
            if cutsites:
                for cs in cutsites:
                    if cs.end == cur_bp:
                        is_after_cutsite = True
                        break

            loc = cur_bp + 1 if is_after_cutsite else cur_bp

            ins_seq = aligned_query[aln_i + 1:aln_j]
            bp_event[loc] = 'I'

            events.append(MutationEvent(
                type=MutationType.INSERTION,
                ref_pos=loc,
                ref_base='-' * gap_cols,
                query_base=ins_seq,
                length=gap_cols,
                in_cutsite_window=_in_any_cutsite_window(loc, cutsites, mutation_window),
                raw_ref_segment='-' * gap_cols,
                raw_query_segment=ins_seq,
            ))

    # Leading insertion (before first ref base)
    if ref_mask[0] > 0:
        ins_seq = aligned_query[:ref_mask[0]]
        loc = 0
        bp_event[0] = 'I'
        events.append(MutationEvent(
            type=MutationType.INSERTION,
            ref_pos=loc,
            ref_base='-' * len(ins_seq),
            query_base=ins_seq,
            length=len(ins_seq),
            in_cutsite_window=_in_any_cutsite_window(loc, cutsites, mutation_window),
            raw_ref_segment='-' * len(ins_seq),
            raw_query_segment=ins_seq,
        ))

    # Trailing insertion (after last ref base)
    last_col = ref_mask[-1]
    if last_col < alen - 1:
        ins_seq = aligned_query[last_col + 1:]
        loc = N - 1
        bp_event[N - 1] = 'I'
        events.append(MutationEvent(
            type=MutationType.INSERTION,
            ref_pos=loc,
            ref_base='-' * len(ins_seq),
            query_base=ins_seq,
            length=len(ins_seq),
            in_cutsite_window=_in_any_cutsite_window(N - 1, cutsites, mutation_window),
            raw_ref_segment='-' * len(ins_seq),
            raw_query_segment=ins_seq,
        ))

    # ── Step 4: Find substitutions (M) ──
    for col_idx in ref_mask:
        ar = aligned_ref[col_idx]
        aq = aligned_query[col_idx]
        if ar != '-' and aq != '-' and ar != aq:
            bp = ref_pos_map[col_idx]
            if bp_event[bp] == 'N':
                bp_event[bp] = 'M'
                events.append(MutationEvent(
                    type=MutationType.SUBSTITUTION,
                    ref_pos=bp,
                    ref_base=ar,
                    query_base=aq,
                    length=1,
                    in_cutsite_window=_in_any_cutsite_window(bp, cutsites, mutation_window),
                    raw_ref_segment=ar,
                    raw_query_segment=aq,
                ))

    events.sort(key=lambda e: e.ref_pos)
    return events


def identify_cas9_events(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: str,
    cutsites: List[CutsiteRegion],
    mutation_window: int = 3,
) -> List[MutationEvent]:
    """MATLAB-compatible identify_cas9_events — compound event merging.

    MATLAB equivalent: @Mutation/identify_cas9_events

    Takes naive events from identify_sequence_events and merges adjacent
    events into compound (COMPLEX) events based on Cas9 edit patterns:
    - Adjacent base positions → merge
    - Events in the same motif site → merge
    - Events in consecutive sites where first is cutsite → merge
    - Overlapping events → merge

    Args:
        aligned_ref: Gapped reference sequence.
        aligned_query: Gapped query sequence.
        ref_seq: Full ungapped reference sequence.
        cutsites: List of CutsiteRegion defining cutsite boundaries.
        mutation_window: Window radius for cutsite proximity detection.

    Returns:
        List of MutationEvent objects with compound events merged.
    """
    naive_events = identify_sequence_events(
        aligned_ref, aligned_query, ref_seq, cutsites, mutation_window
    )

    if len(naive_events) < 2:
        return naive_events

    # ── Determine compound event boundaries ──
    # MATLAB logic: merge events that are adjacent, in same site, or
    # in consecutive sites where first is a cutsite
    compound_starts = [0]
    compound_end = 0

    for i in range(1, len(naive_events)):
        cur = naive_events[i]
        prev_end = naive_events[compound_end]

        cur_start_bp = cur.ref_pos
        prev_end_bp = prev_end.ref_pos + prev_end.length - 1

        # Check adjacency
        is_adjacent = (cur_start_bp == prev_end_bp + 1)

        # Check if in same site
        cur_site = _find_site_for_bp(cur_start_bp, cutsites)
        prev_site = _find_site_for_bp(prev_end_bp, cutsites)
        same_site = (cur_site == prev_site)

        # Check if in consecutive sites (prev is cutsite, cur is next site)
        next_site = False
        if not same_site:
            prev_site_info = _find_site_for_bp(prev_end_bp, cutsites)
            cur_site_info = _find_site_for_bp(cur_start_bp, cutsites)
            next_site = (prev_site_info[0] + 1 == cur_site_info[0] and prev_site_info[1] == 'cutsite')

        # Check overlapping
        is_overlap = (cur_start_bp <= prev_end_bp)

        if is_adjacent or same_site or next_site or is_overlap:
            compound_end = i
        else:
            compound_starts.append(i)
            compound_end = i

    # Build merged events
    merged_events = []
    for idx in range(len(compound_starts)):
        start = compound_starts[idx]
        end = compound_starts[idx + 1] if idx + 1 < len(compound_starts) else len(naive_events)

        if end - start <= 1:
            # Single event, keep as-is
            merged_events.append(naive_events[start])
        else:
            # Merge compound events
            sub = naive_events[start:end]

            ref_start = min(e.ref_pos for e in sub)
            ref_end = max(e.ref_pos + e.length - 1 for e in sub)

            # Build compound sequences
            all_ref_bases = ''.join(e.ref_base for e in sub)
            all_query_bases = ''.join(e.query_base for e in sub)

            compound_event = MutationEvent(
                type=MutationType.COMPLEX,
                ref_pos=ref_start,
                ref_base=all_ref_bases,
                query_base=all_query_bases,
                length=ref_end - ref_start + 1,
                in_cutsite_window=any(e.in_cutsite_window for e in sub),
                raw_ref_segment=all_ref_bases,
                raw_query_segment=all_query_bases,
            )
            merged_events.append(compound_event)

    # Ensure single events that are insertions with both ends non-matching
    # are promoted to COMPLEX (MATLAB: insertion where first and last base
    # of seq_new don't match first and last base of seq_old)
    for i, e in enumerate(merged_events):
        if e.type == MutationType.INSERTION:
            if e.ref_base and e.query_base:
                if e.ref_base[0] != e.query_base[0] and e.ref_base[-1] != e.query_base[-1]:
                    merged_events[i] = MutationEvent(
                        type=MutationType.COMPLEX,
                        ref_pos=e.ref_pos,
                        ref_base=e.ref_base,
                        query_base=e.query_base,
                        length=e.length,
                        in_cutsite_window=e.in_cutsite_window,
                        raw_ref_segment=e.raw_ref_segment,
                        raw_query_segment=e.raw_query_segment,
                    )

    merged_events.sort(key=lambda e: e.ref_pos)
    return merged_events


def classify_bp_event(
    events: List[MutationEvent],
    ref_length: int,
    cutsites: Optional[List[CutsiteRegion]] = None,
) -> Tuple[List[str], List[bool], List[bool]]:
    """MATLAB-compatible classify_bp_event — per-base event classification.

    MATLAB equivalent: @Mutation/classify_bp_event

    Classifies each base pair in the reference as:
    - 'N': Normal (no mutation)
    - 'D': Deleted
    - 'I': Insertion at this bp (or within insertion span)

    Args:
        events: List of MutationEvent (from identify_cas9_events or identify_sequence_events).
        ref_length: Length of the reference sequence.
        cutsites: Optional list of cutsites for insertion span computation.

    Returns:
        Tuple of (bp_event, del_event, ins_event):
        - bp_event: list of chars ('N', 'D', 'I') for each ref position
        - del_event: list of bools, True where deleted
        - ins_event: list of bools, True where insertion occurred
    """
    bp_event = ['N'] * ref_length
    del_event = [False] * ref_length
    ins_event = [False] * ref_length

    for ev in events:
        start = ev.ref_pos
        end = ev.ref_pos + ev.length - 1

        if ev.type == MutationType.SUBSTITUTION:
            # M in MATLAB: both del and ins at this position → 'I'
            del_event[start] = True
            ins_event[start] = True
            bp_event[start] = 'I'

        elif ev.type == MutationType.DELETION:
            for bp in range(start, end + 1):
                del_event[bp] = True
                bp_event[bp] = 'D'

        elif ev.type == MutationType.INSERTION:
            ins_event[start] = True
            # Determine insertion span end
            ins_len = len(ev.query_base.replace('-', ''))
            loc_end = _get_insertion_end(start, ins_len, cutsites, ref_length)
            for bp in range(start, loc_end + 1):
                bp_event[bp] = 'I'

        elif ev.type == MutationType.COMPLEX:
            # Deletion span
            for bp in range(start, end + 1):
                del_event[bp] = True
                bp_event[bp] = 'D'
            # Insertion span (overlapping with deletion)
            ins_len = len(ev.query_base.replace('-', ''))
            ins_event[start] = True
            loc_end = _get_insertion_end(start, ins_len, cutsites, ref_length)
            for bp in range(start, loc_end + 1):
                bp_event[bp] = 'I'

    return bp_event, del_event, ins_event


def _get_insertion_end(
    loc_start: int,
    ins_length: int,
    cutsites: Optional[List[CutsiteRegion]],
    ref_length: int,
) -> int:
    """Determine the end position of an insertion span.

    MATLAB equivalent: get_insertion_end helper in classify_bp_event.
    An insertion extends up to min(loc_start + ins_length - 1, next_cutsite - 1, ref_length).
    """
    loc_end = loc_start + ins_length - 1
    if cutsites:
        for cs in sorted(cutsites, key=lambda c: c.start):
            if cs.start > loc_start:
                loc_end = min(loc_end, cs.start - 1)
                break
    return min(loc_end, ref_length - 1)


def annotate_mutation(mutation: MutationEvent, full: bool = False) -> str:
    """Format a single mutation event as an HGVS-like annotation string.

    Follows the same convention as MATLAB's Mutation.annotate().

    Substitution:     "41A>T"
    Deletion:         "41_43del"
    Insertion (short): "42_43ins3"
    Insertion (full):  "42_43insXYZ"
    Complex (short):   "41_43delins5"
    Complex (full):    "41_43delinsXYZ"

    Args:
        mutation: The mutation event to annotate.
        full: If True, include inserted sequence rather than just length.

    Returns:
        HGVS-like annotation string.
    """
    pos = mutation.ref_pos
    if mutation.type == MutationType.SUBSTITUTION:
        return f"{pos + 1}{mutation.ref_base}>{mutation.query_base}"
    elif mutation.type == MutationType.DELETION:
        end_pos = pos + mutation.length
        if mutation.length == 1:
            return f"{pos + 1}del"
        return f"{pos + 1}_{end_pos}del"
    elif mutation.type == MutationType.INSERTION:
        ref_base = mutation.ref_base
        query_base = mutation.query_base
        if full:
            inserted = query_base.replace("-", "")
            return f"{pos + 1}_{pos + 2}ins{inserted}"
        else:
            inserted_len = len(query_base.replace("-", ""))
            return f"{pos + 1}_{pos + 2}ins{inserted_len}"
    elif mutation.type == MutationType.COMPLEX:
        if full:
            inserted = mutation.query_base.replace("-", "")
            end_pos = pos + mutation.length
            return f"{pos + 1}_{end_pos}delins{inserted}"
        else:
            inserted_len = len(mutation.query_base.replace("-", ""))
            end_pos = pos + mutation.length
            return f"{pos + 1}_{end_pos}delins{inserted_len}"
    return ""


def annotate_mutations(mutations: List[MutationEvent], full: bool = False) -> str:
    """Format a list of mutations as a compact annotation string.

    Multiple events are joined with ";".
    """
    return ";".join(annotate_mutation(m, full) for m in mutations)


def extract_mutations(
    aligned_ref: str,
    aligned_query: str,
    cutsites: Optional[List[CutsiteRegion]] = None,
    mutation_window: int = 3,
) -> List[MutationEvent]:
    """从比对结果中提取所有突变事件

    遍历比对（aligned_ref ↔ aligned_query），识别：
      1. 替换（substitution）: 两列都有碱基但不同
      2. 删除（deletion）: query 为 '-'，ref 有碱基
      3. 插入（insertion）: ref 为 '-'，query 有碱基
      4. 复合（complex）: 相邻的插入+删除合并为复合事件

    参数:
        aligned_ref: 比对后的参考序列（含 '-'）
        aligned_query: 比对后的查询序列（含 '-'）
        cutsites: cutsite区域列表（用于判断突变是否在窗口内）
        mutation_window: cutsite窗口半径（bp）

    返回:
        MutationEvent 列表，按 ref 位置排序
    """
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return []

    events: List[MutationEvent] = []
    i = 0
    alen = len(aligned_ref)

    # 构建 ref 位置映射：aligned列 → ref序列上的坐标
    ref_pos_map, _ = _build_ref_pos_map_full(aligned_ref)

    while i < alen:
        ar = aligned_ref[i]
        aq = aligned_query[i]

        # ── 情况1: 替换（两列都有碱基但不相同）──
        if ar != '-' and aq != '-' and ar != aq:
            ref_pos = ref_pos_map[i]
            event = MutationEvent(
                type=MutationType.SUBSTITUTION,
                ref_pos=ref_pos,
                ref_base=ar,
                query_base=aq,
                length=1,
                in_cutsite_window=_in_any_cutsite_window(ref_pos, cutsites, mutation_window) if cutsites else True,
                raw_ref_segment=ar,
                raw_query_segment=aq,
            )
            events.append(event)
            i += 1

        # ── 情况2: 删除（query 为 '-'）──
        elif ar != '-' and aq == '-':
            ref_start = ref_pos_map[i]
            j = i
            while j < alen and aligned_ref[j] != '-' and aligned_query[j] == '-':
                j += 1
            length = j - i
            ref_end = ref_start + length
            seg_ref = aligned_ref[i:j]
            event = MutationEvent(
                type=MutationType.DELETION,
                ref_pos=ref_start,
                ref_base=seg_ref,
                query_base='-' * length,
                length=length,
                in_cutsite_window=_in_any_cutsite_window(ref_start, cutsites, mutation_window) or \
                              _in_any_cutsite_window(ref_end - 1, cutsites, mutation_window) if cutsites else True,
                raw_ref_segment=seg_ref,
                raw_query_segment='-' * length,
            )
            events.append(event)
            i = j

        # ── 情况3: 插入（ref 为 '-'）──
        elif ar == '-' and aq != '-':
            ref_pos = ref_pos_map[i] if i > 0 else 0
            # 如果是首列没有 ref_pos_map，用第一个非gap的位置
            if ref_pos < 0:
                for k in range(i, alen):
                    p = ref_pos_map[k]
                    if p >= 0:
                        ref_pos = p
                        break
                    elif k < alen - 1 and ref_pos_map[k + 1] >= 0:
                        ref_pos = ref_pos_map[k + 1]
                        break
            j = i
            while j < alen and aligned_ref[j] == '-' and aligned_query[j] != '-':
                j += 1
            length = j - i
            seg_query = aligned_query[i:j]
            event = MutationEvent(
                type=MutationType.INSERTION,
                ref_pos=max(0, ref_pos),
                ref_base='-' * length,
                query_base=seg_query,
                length=length,
                in_cutsite_window=_in_any_cutsite_window(max(0, ref_pos), cutsites, mutation_window) if cutsites else True,
                raw_ref_segment='-' * length,
                raw_query_segment=seg_query,
            )
            events.append(event)
            i = j

        else:
            # 匹配列，跳过
            i += 1

    # 合并相邻的插入+删除为复合事件
    events = _merge_adjacent_indels(events)

    return events


def _in_any_cutsite_window(
    ref_pos: int,
    cutsites: List[CutsiteRegion],
    window: int = 3,
) -> bool:
    """判断一个参考序列坐标是否在任何cutsite的窗口内"""
    if not cutsites:
        return True  # 没有cutsite定义时视为全在窗口内
    for cs in cutsites:
        if cs.start - window <= ref_pos <= cs.end + window:
            return True
    return False


def _merge_adjacent_indels(events: List[MutationEvent]) -> List[MutationEvent]:
    """合并相邻的插入和删除为复合事件

    例如：[DEL(3-5), INS(5-6)] → COMPLEX(3-6)
    如果插入和删除的位置相邻且方向兼容，合并为一个复合突变事件。
    """
    if len(events) < 2:
        return events

    merged = []
    i = 0
    while i < len(events):
        if i + 1 < len(events):
            e1 = events[i]
            e2 = events[i + 1]
            # 插入相邻删除 → 复合
            if (e1.type == MutationType.INSERTION and e2.type == MutationType.DELETION) or \
               (e1.type == MutationType.DELETION and e2.type == MutationType.INSERTION):
                # 检查是否相邻
                adj = False
                if e1.type == MutationType.INSERTION:
                    # INS在ref_pos X，DEL在ref_pos Y，如果X≈Y则是相邻
                    adj = abs(e1.ref_pos - e2.ref_pos) <= max(e1.length, e2.length) + 1
                else:
                    adj = abs(e1.ref_pos - e2.ref_pos) <= max(e1.length, e2.length)

                if adj:
                    merged.append(MutationEvent(
                        type=MutationType.COMPLEX,
                        ref_pos=min(e1.ref_pos, e2.ref_pos),
                        length=e1.length + e2.length,
                        in_cutsite_window=e1.in_cutsite_window or e2.in_cutsite_window,
                        ref_base=e1.ref_base + e2.ref_base,
                        query_base=e1.query_base + e2.query_base,
                        raw_ref_segment=e1.raw_ref_segment + e2.raw_ref_segment,
                        raw_query_segment=e1.raw_query_segment + e2.raw_query_segment,
                    ))
                    i += 2
                    continue
        merged.append(events[i])
        i += 1
    return merged


def classify_mutation_type(
    stats: AlignmentStats,
) -> str:
    """根据比对统计信息分类突变类型（用于报告）

    返回字符串标签如: "only_substitution", "insertion_and_deletion" 等
    """
    has_ins = stats.gaps_in_ref > 0
    has_del = stats.gaps_in_query > 0
    has_sub = stats.mismatches > 0

    if has_ins and not has_del and not has_sub:
        return "only_insertion"
    elif has_del and not has_ins and not has_sub:
        return "only_deletion"
    elif has_sub and not has_ins and not has_del:
        return "only_substitution"
    elif has_ins and has_del and not has_sub:
        return "insertion_and_deletion"
    elif has_ins and has_sub and not has_del:
        return "insertion_and_substitution"
    elif has_del and has_sub and not has_ins:
        return "deletion_and_substitution"
    elif has_ins and has_del and has_sub:
        return "insertion_deletion_substitution"
    return "unmutated"


def build_mutation_summary(
    results: List['AlignmentResult'],
) -> dict:
    """构建突变类型汇总统计

    统计各类突变的数量（序列数 + Reads数），
    以及 indel 长度分布和点突变总数。
    """
    from collections import defaultdict

    type_counts = defaultdict(lambda: {"sequences": 0, "reads": 0})
    ins_lengths = []
    del_lengths = []
    total_mismatches = 0

    for r in results:
        if not r.success or r.stats is None:
            continue
        rc = r.query.readCount
        mtype = classify_mutation_type(r.stats)
        type_counts[mtype]["sequences"] += 1
        type_counts[mtype]["reads"] += rc

        # indel 长度收集
        for block in r.stats.gap_blocks_ref:
            ins_lengths.append(block)
        for block in r.stats.gap_blocks_query:
            del_lengths.append(block)
        total_mismatches += r.stats.mismatches

    return {
        "type_counts": dict(type_counts),
        "ins_lengths": ins_lengths,
        "del_lengths": del_lengths,
        "total_mismatches": total_mismatches,
    }


def format_mutations_for_display(
    mutations: List[MutationEvent],
    max_events: int = 20,
) -> str:
    """将突变事件列表格式化为可读字符串（用于调试/日志）"""
    if not mutations:
        return "未检测到突变"

    lines = []
    for i, m in enumerate(mutations[:max_events]):
        pos_info = f"@pos{m.ref_pos}" if m.ref_pos >= 0 else ""
        cut_info = " [窗口内]" if m.in_cutsite_window else " [窗口外]"
        if m.type == MutationType.SUBSTITUTION:
            lines.append(f"  {i+1}. 替换 {pos_info}: {m.ref_base}→{m.query_base}{cut_info}")
        elif m.type == MutationType.DELETION:
            lines.append(f"  {i+1}. 删除 {pos_info}: -{m.length}bp ({m.ref_base}){cut_info}")
        elif m.type == MutationType.INSERTION:
            lines.append(f"  {i+1}. 插入 {pos_info}: +{m.length}bp ({m.query_base}){cut_info}")
        elif m.type == MutationType.COMPLEX:
            lines.append(f"  {i+1}. 复合 {pos_info}: {m.length}bp{cut_info}")

    if len(mutations) > max_events:
        lines.append(f"  ... 还有 {len(mutations) - max_events} 个事件未显示")

    return "\n".join(lines)
