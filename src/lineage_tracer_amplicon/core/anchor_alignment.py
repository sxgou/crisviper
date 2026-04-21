"""
Anchor-based alignment for CRISPR amplicon analysis.
"""

from typing import List, Dict, Any, Optional, Tuple, Set
import edlib
import re
from ..core.structures import Anchor, AnchorMatch


class AnchorIndex:
    """Index for fast anchor matching using k-mers."""
    
    def __init__(self, k: int = 9):
        """
        Initialize anchor index.
        
        Args:
            k: k-mer size for indexing
        """
        self.k = k
        self.index: Dict[str, List[Tuple[int, int]]] = {}  # kmer -> [(anchor_id, position)]
        self.anchors: List[Anchor] = []
    
    def build(self, anchors: List[Anchor]):
        """
        Build index from anchor list.
        
        Args:
            anchors: List of Anchor objects
        """
        self.anchors = anchors
        self.index = {}
        
        for anchor_id, anchor in enumerate(anchors):
            anchor_seq = anchor.seq
            # Extract all k-mers from anchor sequence
            for pos in range(len(anchor_seq) - self.k + 1):
                kmer = anchor_seq[pos:pos+self.k]
                if kmer not in self.index:
                    self.index[kmer] = []
                self.index[kmer].append((anchor_id, pos))
    
    def query(self, kmer: str) -> List[Tuple[int, int]]:
        """
        Query k-mer in index.
        
        Args:
            kmer: k-mer to query
            
        Returns:
            List of (anchor_id, position) tuples
        """
        return self.index.get(kmer, [])


def semi_global_align(
    query_seq: str,
    target_seq: str
) -> Optional[Dict[str, Any]]:
    """
    Perform semi-global alignment using edlib.
    Target is the anchor sequence (shorter), query is the read (longer).
    Mode 'HW' (half-way): query can start/end anywhere, target must be fully aligned.
    
    Args:
        query_seq: Query sequence (read)
        target_seq: Target sequence (anchor)
        
    Returns:
        Alignment result dictionary or None if alignment fails
    """
    # Use edlib with mode="HW" (semi-global: target full length, query partial)
    result = edlib.align(target_seq, query_seq, mode="HW", task="path")
    
    if result["editDistance"] == -1:
        return None
    
    edit_distance = result["editDistance"]
    cigar = result["cigar"]
    
    # Parse CIGAR to get query start/end positions
    query_pos = 0
    target_pos = 0
    query_start = None
    query_end = None
    
    ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    for length_str, op in ops:
        length = int(length_str)
        if op in ('M', '=', 'X'):
            if query_start is None:
                query_start = query_pos
            query_pos += length
            query_end = query_pos
            target_pos += length
        elif op == 'I':
            query_pos += length
        elif op == 'D':
            target_pos += length
        elif op in ('S', 'H'):
            query_pos += length
    
    if query_start is None:
        return None
    
    aligned_length = target_pos  # Should equal len(target_seq)
    matches = aligned_length - edit_distance
    identity = matches / len(target_seq) if len(target_seq) > 0 else 0.0
    
    return {
        'query_start': query_start,
        'query_end': query_end,
        'target_start': 0,
        'target_end': len(target_seq),
        'identity': identity,
        'edit_distance': edit_distance,
        'cigar': cigar,
        'aligned_length': aligned_length
    }


def match_anchors(
    read_seq: str,
    anchors: List[Anchor],
    anchor_index: AnchorIndex,
    min_identity: float = 0.8,
    k: int = 9
) -> List[AnchorMatch]:
    """
    Match anchors to read sequence.
    
    Args:
        read_seq: Read sequence
        anchors: List of Anchor objects
        anchor_index: Pre-built AnchorIndex
        min_identity: Minimum identity threshold
        k: k-mer size
        
    Returns:
        List of AnchorMatch objects
    """
    # Step 1: Extract read k-mers and find candidate anchors
    candidate_anchors: Set[int] = set()
    k = anchor_index.k
    for i in range(len(read_seq) - k + 1):
        kmer = read_seq[i:i+k]
        matches = anchor_index.query(kmer)
        for anchor_id, _ in matches:
            candidate_anchors.add(anchor_id)
    
    # Step 2: Perform semi-global alignment for each candidate anchor
    matches = []
    for anchor_id in candidate_anchors:
        anchor = anchors[anchor_id]
        anchor_seq = anchor.seq
        
        alignment = semi_global_align(read_seq, anchor_seq)
        
        if alignment is None:
            continue
        
        identity = alignment['identity']
        # For short anchors (<10bp), be more lenient
        if len(anchor_seq) < 10:
            min_identity_adj = min(0.7, min_identity)
        else:
            min_identity_adj = min_identity
            
        if identity >= min_identity_adj:
            # Check alignment covers most of anchor
            aligned_anchor_len = alignment['aligned_length']
            min_required_len = max(3, int(len(anchor_seq) * 0.7))
            if aligned_anchor_len >= min_required_len:
                match = AnchorMatch(
                    anchor=anchor,
                    read_start=alignment['query_start'],
                    read_end=alignment['query_end'],
                    identity=identity,
                    cigar=alignment['cigar']
                )
                matches.append(match)
    
    # Step 3: Filter overlapping matches (keep best non-overlapping)
    matches = resolve_overlapping_matches(matches)
    
    # Step 4: Sort by reference position
    matches.sort(key=lambda m: m.anchor.ref_start)
    
    return matches


def resolve_overlapping_matches(matches: List[AnchorMatch]) -> List[AnchorMatch]:
    """
    Resolve overlapping anchor matches, keeping best non-overlapping set.
    
    Args:
        matches: List of AnchorMatch objects
        
    Returns:
        Filtered list of non-overlapping AnchorMatch objects
    """
    if not matches:
        return []
    
    matches.sort(key=lambda m: m.read_start)
    
    filtered = []
    current = matches[0]
    
    for match in matches[1:]:
        overlap_start = max(current.read_start, match.read_start)
        overlap_end = min(current.read_end, match.read_end)
        
        if overlap_end <= overlap_start:
            filtered.append(current)
            current = match
        else:
            if match.identity > current.identity:
                current = match
            elif match.identity == current.identity and match.read_start < current.read_start:
                current = match
    
    filtered.append(current)
    
    return filtered


def infer_deletions(
    anchor_matches: List[AnchorMatch],
    ref_length: int,
    min_gap_size: int = 5
) -> List[Tuple[int, int]]:
    """
    Infer deletion intervals based on anchor matches.
    Only report gaps larger than min_gap_size as deletions.
    
    Args:
        anchor_matches: List of AnchorMatch objects sorted by ref_start
        ref_length: Total reference length
        min_gap_size: Minimum gap size to report as deletion
        
    Returns:
        List of (start, end) deletion intervals (0-based exclusive)
    """
    deletions = []
    prev_ref_end = 0
    
    sorted_matches = sorted(anchor_matches, key=lambda m: m.anchor.ref_start)
    
    for match in sorted_matches:
        anchor = match.anchor
        if anchor.ref_start > prev_ref_end:
            gap_size = anchor.ref_start - prev_ref_end
            if gap_size >= min_gap_size:
                deletions.append((prev_ref_end, anchor.ref_start))
        prev_ref_end = max(prev_ref_end, anchor.ref_end)
    
    if prev_ref_end < ref_length:
        gap_size = ref_length - prev_ref_end
        if gap_size >= min_gap_size:
            deletions.append((prev_ref_end, ref_length))
    
    # Merge adjacent deletion intervals
    merged = []
    for interval in deletions:
        if not merged:
            merged.append(interval)
        else:
            last_start, last_end = merged[-1]
            if interval[0] <= last_end:
                merged[-1] = (last_start, max(last_end, interval[1]))
            else:
                merged.append(interval)
    
    return merged


def find_closest_anchors(
    anchor_matches: List[AnchorMatch],
    ref_position: int
) -> Tuple[Optional[AnchorMatch], Optional[AnchorMatch]]:
    """
    Find closest anchor before and after a reference position.
    
    Args:
        anchor_matches: List of AnchorMatch objects
        ref_position: Reference position
        
    Returns:
        Tuple of (anchor_before, anchor_after)
    """
    anchor_before = None
    anchor_after = None
    min_before_dist = float('inf')
    min_after_dist = float('inf')
    
    for match in anchor_matches:
        anchor = match.anchor
        dist = ref_position - anchor.ref_end
        
        if dist >= 0 and dist < min_before_dist:
            anchor_before = match
            min_before_dist = dist
        
        dist = anchor.ref_start - ref_position
        if dist >= 0 and dist < min_after_dist:
            anchor_after = match
            min_after_dist = dist
    
    return anchor_before, anchor_after