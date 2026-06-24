"""crisviper/mutation.py — Independent mutation extraction and classification module.

Extracts structured mutation events from pairwise alignments produced by the
alignment pipeline. This module is the core "mutation identification" step:
it converts raw aligned sequences into typed, position-annotated mutation events
(substitutions, deletions, insertions, and complex indels).
"""

from typing import List, Optional, Tuple
from crisviper.models import MutationEvent, MutationType, AlignmentStats
from crisviper.config import CutsiteRegion


def _build_ref_pos_map_full(aligned_ref: str):
    """Build alignment-column-to-reference-position mapping.

    Maps each column in the gapped alignment to its corresponding 0-based
    position on the ungapped reference sequence. Gap columns (dashes in the
    reference sequence) map to -1.

    Returns:
        Tuple of (pos_map, total):
        - pos_map: list of length len(aligned_ref), -1 for gap columns
        - total: number of non-gap reference bases
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
    """Extract all mutation events from a pairwise alignment.

    Walks through the aligned reference and query sequences column by column
    to identify four mutation types:
      1. Substitution: both columns have bases but they differ
      2. Deletion: query has a gap ('-'), reference has a base
      3. Insertion: reference has a gap ('-'), query has a base
      4. Complex (INDEL): adjacent insertion + deletion blocks are merged
         into a single composite event via greedy grouping

    Args:
        aligned_ref: Gapped reference sequence from alignment.
        aligned_query: Gapped query sequence from alignment.
        cutsites: List of cutsite regions for in-window annotation.
        sub_window: Radius around cutsites for window annotation (bp).

    Returns:
        List of MutationEvent objects sorted by reference position.
    """
    if len(aligned_ref) != len(aligned_query) or len(aligned_ref) == 0:
        return []

    events: List[MutationEvent] = []
    i = 0
    alen = len(aligned_ref)

    # Build reference position map: alignment column → ref sequence coordinate
    ref_pos_map, _ = _build_ref_pos_map_full(aligned_ref)

    while i < alen:
        ar = aligned_ref[i]
        aq = aligned_query[i]

        # ── Case 1: Substitution (both columns have bases, but they differ) ──
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

        # ── Case 2: Deletion (query has gap, reference has base) ──
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
                in_cutsite_window=(_in_any_cutsite_window(ref_start, cutsites, sub_window) or
                                  _in_any_cutsite_window(ref_end - 1, cutsites, sub_window)) if cutsites else True,
                raw_ref_segment=seg_ref,
                raw_query_segment='-' * length,
            )
            events.append(event)
            i = j

        # ── Case 3: Insertion (reference has gap, query has base) ──
        elif ar == '-' and aq != '-':
            # Insertion position: use the ref base immediately before the
            # insertion gap (i > 0 → ref_pos_map[i-1] is always >= 0 since
            # this is the first column of a new insertion block).
            if i > 0:
                ref_pos = max(0, ref_pos_map[i - 1])
            else:
                ref_pos = 0
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
            # Matching column — skip
            i += 1

    # Merge adjacent insertion+deletion blocks into composite INDEL events
    events = _merge_adjacent_indels(events)

    return events


def _in_any_cutsite_window(
    ref_pos: int,
    cutsites: List[CutsiteRegion],
    window: int = 3,
) -> bool:
    """Check whether a reference position falls within any cutsite's extended window.

    Args:
        ref_pos: 0-based reference position to check.
        cutsites: List of cutsite regions.
        window: Radius to extend beyond each cutsite boundary (bp).

    Returns:
        True if ref_pos is within any cutsite region ± window.
        When cutsites list is empty, returns True (all positions treated as in-window).
    """
    if not cutsites:
        return True  # No cutsites defined → treat all positions as in-window
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
        if e.type in (MutationType.DELETION, MutationType.INDEL):
            return e.ref_pos + e.length
        return e.ref_pos + 1  # INS, SUB: single ref position

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

        if has_indel and len(group) > 1:
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
        elif has_indel:
            # Single insertion or deletion — keep original type
            merged.append(group[0])
        else:
            merged.extend(group)

        i = j

    return merged


def classify_mutation_type(
    stats: AlignmentStats,
) -> str:
    """Classify the mutation type of an alignment based on its statistics.

    Uses the presence/absence of insertions (gaps_in_ref), deletions
    (gaps_in_query), and substitutions (mismatches) to produce a category
    label for reporting purposes.

    Returns:
        String label such as "only_substitution", "insertion_and_deletion",
        "unmutated", etc. Full list:
        - only_insertion / only_deletion / only_substitution
        - insertion_and_deletion / insertion_and_substitution / deletion_and_substitution
        - insertion_deletion_substitution / unmutated
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
    """Build aggregate mutation type summary across all alignment results.

    Counts sequences and reads per mutation category (using classify_mutation_type),
    collects indel length distributions from gap block statistics, and
    sums total point mutations (mismatches).

    Returns:
        dict with keys:
        - type_counts: dict of {mutation_label: {"sequences": N, "reads": N}}
        - ins_lengths: list of insertion lengths (bp), one per gap block in ref
        - del_lengths: list of deletion lengths (bp), one per gap block in query
        - total_mismatches: sum of all mismatches across successful alignments
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

        # Collect indel lengths for distribution analysis
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


