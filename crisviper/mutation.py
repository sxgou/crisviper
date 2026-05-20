"""crisviper/mutation.py — 突变识别独立模块

从比对结果中显式提取结构化突变事件。
这是管道的核心步骤——"突变识别"。
"""

from typing import List, Optional, Tuple
from crisviper.models import MutationEvent, MutationType, AlignmentStats
from crisviper.config import CutsiteRegion


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

        elif ev.type == MutationType.INDEL:
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
    elif mutation.type == MutationType.INDEL:
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
    sub_window: int = 3,
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
        sub_window: cutsite窗口半径（bp）

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
                in_cutsite_window=_in_any_cutsite_window(ref_pos, cutsites, sub_window) if cutsites else True,
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
                in_cutsite_window=_in_any_cutsite_window(ref_start, cutsites, sub_window) or \
                              _in_any_cutsite_window(ref_end - 1, cutsites, sub_window) if cutsites else True,
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
                in_cutsite_window=_in_any_cutsite_window(max(0, ref_pos), cutsites, sub_window) if cutsites else True,
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
    """Greedy grouping merge of adjacent events into INDEL.

    Scans sorted events, collecting consecutive adjacent events into groups.
    Groups containing at least one INS or DEL are merged into a single INDEL.
    Pure-SUB groups stay as individual events.

    Adjacency: two events are adjacent if their ref intervals overlap or
    are within 1bp of each other.

    For DEL:  interval = [ref_pos, ref_pos + length)
    For INS:  interval = [ref_pos, ref_pos + 1)     (no ref span)
    For SUB:  interval = [ref_pos, ref_pos + 1)
    """
    if len(events) < 2:
        return events

    def exclusive_end(e: MutationEvent) -> int:
        if e.type == MutationType.DELETION:
            return e.ref_pos + e.length
        return e.ref_pos + 1  # INS and SUB don't span ref positions

    def intervals_adjacent(a: MutationEvent, b: MutationEvent) -> bool:
        a_end = exclusive_end(a)
        b_end = exclusive_end(b)
        return max(a.ref_pos, b.ref_pos) <= min(a_end, b_end) + 1

    merged = []
    i = 0
    while i < len(events):
        group = [events[i]]
        j = i + 1
        while j < len(events) and intervals_adjacent(group[-1], events[j]):
            group.append(events[j])
            j += 1

        has_indel = any(
            e.type in (MutationType.INSERTION, MutationType.DELETION)
            for e in group
        )

        if has_indel:
            ref_start = min(e.ref_pos for e in group)
            ref_end = max(exclusive_end(e) for e in group)
            length = ref_end - ref_start
            all_ref = ''.join(e.ref_base for e in group)
            all_query = ''.join(e.query_base for e in group)
            merged.append(MutationEvent(
                type=MutationType.INDEL,
                ref_pos=ref_start,
                ref_base=all_ref,
                query_base=all_query,
                length=length,
                in_cutsite_window=any(e.in_cutsite_window for e in group),
                raw_ref_segment=all_ref,
                raw_query_segment=all_query,
            ))
        else:
            merged.extend(group)

        i = j

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
        elif m.type == MutationType.INDEL:
            lines.append(f"  {i+1}. 复合 {pos_info}: {m.length}bp{cut_info}")

    if len(mutations) > max_events:
        lines.append(f"  ... 还有 {len(mutations) - max_events} 个事件未显示")

    return "\n".join(lines)
