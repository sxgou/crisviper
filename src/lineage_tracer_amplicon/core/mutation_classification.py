"""
Mutation classification and correction for CRISPR amplicon analysis.
"""

from typing import List, Dict, Any, Optional, Tuple
from ..core.structures import (
    DeletionEvent, InsertionEvent, SNVEvent, Target, AmpliconStructure
)


def detect_microhomology(
    fragment: str,
    flank_left: str,
    flank_right: str,
    min_homology: int = 3,
    max_mismatch: int = 1
) -> Optional[Tuple[str, str]]:
    """
    Detect microhomology in fragment between left and right flanks.
    
    Args:
        fragment: Inserted fragment sequence
        flank_left: Left flanking sequence
        flank_right: Right flanking sequence
        min_homology: Minimum homology length
        max_mismatch: Maximum allowed mismatches
        
    Returns:
        Tuple of (left_homology, right_homology) or None
    """
    if len(fragment) < 2 * min_homology:
        return None
    
    for split in range(min_homology, len(fragment) - min_homology + 1):
        left_part = fragment[:split]
        right_part = fragment[split:]
        
        left_match = False
        if len(flank_left) >= len(left_part):
            left_suffix = flank_left[-len(left_part):]
            mismatches = sum(1 for a, b in zip(left_suffix, left_part) if a != b)
            left_match = mismatches <= max_mismatch
        
        right_match = False
        if len(flank_right) >= len(right_part):
            right_prefix = flank_right[:len(right_part)]
            mismatches = sum(1 for a, b in zip(right_prefix, right_part) if a != b)
            right_match = mismatches <= max_mismatch
        
        if left_match and right_match:
            return left_part, right_part
    
    return None


def endswith_mismatch(s: str, suffix: str, max_mismatch: int = 1) -> bool:
    """Check if string ends with suffix allowing mismatches."""
    if len(s) < len(suffix):
        return False
    actual_suffix = s[-len(suffix):]
    mismatches = sum(1 for a, b in zip(actual_suffix, suffix) if a != b)
    return mismatches <= max_mismatch


def startswith_mismatch(s: str, prefix: str, max_mismatch: int = 1) -> bool:
    """Check if string starts with prefix allowing mismatches."""
    if len(s) < len(prefix):
        return False
    actual_prefix = s[:len(prefix)]
    mismatches = sum(1 for a, b in zip(actual_prefix, prefix) if a != b)
    return mismatches <= max_mismatch


def classify_deletion_with_insertion(
    read_fragment: str,
    ref_seq: str,
    del_start: int,
    del_end: int,
    min_homology: int = 3
) -> Dict[str, Any]:
    """
    Classify deletion with insertion event.
    """
    if len(read_fragment) == 0:
        return {
            'type': 'pure_deletion',
            'del_length': del_end - del_start,
            'inserted_seq': '',
            'is_mmej': False,
            'mh_left': '',
            'mh_right': ''
        }
    
    flank_left_start = max(0, del_start - 20)
    flank_left = ref_seq[flank_left_start:del_start]
    
    flank_right_end = min(len(ref_seq), del_end + 20)
    flank_right = ref_seq[del_end:flank_right_end]
    
    mh = detect_microhomology(read_fragment, flank_left, flank_right, min_homology)
    
    if mh:
        left_homology, right_homology = mh
        return {
            'type': 'mmej',
            'del_length': del_end - del_start,
            'inserted_seq': read_fragment,
            'is_mmej': True,
            'mh_left': left_homology,
            'mh_right': right_homology
        }
    else:
        return {
            'type': 'deletion_with_insertion',
            'del_length': del_end - del_start,
            'inserted_seq': read_fragment,
            'is_mmej': False,
            'mh_left': '',
            'mh_right': ''
        }


def is_in_valid_snv_window(
    ref_pos: int,
    structure: AmpliconStructure,
    window_size: int = 3
) -> Tuple[bool, Optional[str]]:
    """
    Check if reference position is in valid SNV window of any target.
    """
    for target in structure.targets:
        window_start = max(0, target.cutsite_start - window_size)
        window_end = min(len(structure.reference), target.cutsite_end + window_size)
        
        if window_start <= ref_pos < window_end:
            return True, target.name
    
    return False, None


def process_point_mutations(
    read_seq: str,
    ref_seq: str,
    anchor_matches: List[Dict[str, Any]],
    structure: AmpliconStructure,
    window_size: int = 3
) -> List[SNVEvent]:
    """
    Process point mutations, filtering only those in valid windows.
    """
    valid_snvs = []
    
    for match in anchor_matches:
        if not isinstance(match, dict):
            continue
            
        read_start = match.get('read_start', 0)
        ref_start = match.get('ref_start', 0)
        length = match.get('length', 0)
        
        for i in range(length):
            read_pos = read_start + i
            ref_pos = ref_start + i
            
            if read_pos >= len(read_seq) or ref_pos >= len(ref_seq):
                continue
                
            read_base = read_seq[read_pos]
            ref_base = ref_seq[ref_pos]
            
            if read_base != ref_base:
                is_valid, target_name = is_in_valid_snv_window(ref_pos, structure, window_size)
                
                if is_valid:
                    snv = SNVEvent(
                        ref_pos=ref_pos,
                        ref_base=ref_base,
                        alt_base=read_base,
                        target_name=target_name
                    )
                    valid_snvs.append(snv)
    
    return valid_snvs


def calculate_coverage(
    target_start: int,
    target_end: int,
    event_start: int,
    event_end: int
) -> float:
    """Calculate coverage of target by event."""
    overlap_start = max(target_start, event_start)
    overlap_end = min(target_end, event_end)
    
    if overlap_end <= overlap_start:
        return 0.0
    
    overlap_length = overlap_end - overlap_start
    target_length = target_end - target_start
    
    return overlap_length / target_length if target_length > 0 else 0.0


def determine_target_state(
    target: Target,
    deletion_events: List[DeletionEvent],
    insertion_events: List[InsertionEvent],
    snv_events: List[SNVEvent],
    coverage_threshold: float = 0.9
) -> str:
    """
    Determine mutation state for a target.
    """
    # 1. Check if target is completely covered by a large deletion
    for del_event in deletion_events:
        target_coverage = calculate_coverage(
            target.ref_start, target.ref_end,
            del_event.ref_start, del_event.ref_end
        )
        if target_coverage >= coverage_threshold:
            return "DELETED"
    
    # 2. Check for local indels/SNVs in cutsite region
    cutsite_center = target.cutsite_start + 3  # PAM upstream 3bp
    cutsite_region_start = max(target.ref_start, cutsite_center - 5)
    cutsite_region_end = min(target.ref_end, cutsite_center + 5)
    
    # Find deletion events that have breakpoint in cutsite region
    local_del_event = None
    for del_event in deletion_events:
        if (cutsite_region_start <= del_event.ref_start <= cutsite_region_end or
            cutsite_region_start <= del_event.ref_end <= cutsite_region_end):
            local_del_event = del_event
            break
    
    # Find insertion events in cutsite region
    local_ins_event = None
    for ins_event in insertion_events:
        if cutsite_region_start <= ins_event.ref_pos <= cutsite_region_end:
            local_ins_event = ins_event
            break
    
    # Find SNV events in this target
    target_snvs = [s for s in snv_events if s.target_name == target.name]
    
    has_local_del = local_del_event is not None
    has_local_ins = local_ins_event is not None
    has_snv = len(target_snvs) > 0
    
    # 3. If no events at all, it's wild-type
    if not has_local_del and not has_local_ins and not has_snv:
        return "WT"
    
    # 4. Build complex state string
    if has_local_del and has_local_ins and has_snv:
        return "COMPLEX:del+ins+snv"
    elif has_local_del and has_local_ins:
        return "COMPLEX:del+ins"
    elif has_local_del and has_snv:
        return "COMPLEX:del+snv"
    elif has_local_ins and has_snv:
        return "COMPLEX:ins+snv"
    elif has_local_del:
        del_len = local_del_event.ref_end - local_del_event.ref_start
        return f"DEL:{del_len}"
    elif has_local_ins:
        return f"INS:{local_ins_event.inserted_seq}"
    elif has_snv:
        snv_strs = [f"{s.ref_base}>{s.alt_base}@{s.ref_pos-target.ref_start}" for s in target_snvs]
        return f"SNV:{','.join(snv_strs)}"
    
    return "WT"