"""
Primer validation and read filtering.
"""

from typing import Optional, Tuple, Dict, Any
import edlib
from ..core.structures import PrimerMatch


def detect_primer_at_end(
    read_seq: str, 
    primer_seq: str, 
    end: str = 'start', 
    max_mismatch_rate: float = 0.15,
    min_overlap_ratio: float = 0.8
) -> Optional[Dict[str, Any]]:
    """
    Detect primer sequence at specified end of read.
    
    Args:
        read_seq: Read sequence
        primer_seq: Primer sequence
        end: 'start' for 5' end, 'end' for 3' end
        max_mismatch_rate: Maximum allowed mismatch rate (0.0-1.0)
        min_overlap_ratio: Minimum overlap ratio with primer (0.0-1.0)
        
    Returns:
        Dictionary with match info or None if no match found
    """
    primer_len = len(primer_seq)
    max_mismatch_allowed = int(primer_len * max_mismatch_rate)
    
    # Determine search window
    if end == 'start':
        # Search at 5' end, allow a few extra bases before primer
        search_window = read_seq[:primer_len + 10]  # Allow up to 10 extra bases
        search_offset = 0
    else:  # 'end'
        # Search at 3' end
        search_window = read_seq[-(primer_len + 10):]
        search_offset = max(0, len(read_seq) - len(search_window))
    
    best_match = None
    best_score = 0
    
    # Slide window and calculate edit distance
    for i in range(len(search_window) - primer_len + 1):
        subseq = search_window[i:i+primer_len]
        
        # Use edlib for fast alignment
        result = edlib.align(primer_seq, subseq, task="path")
        
        # Calculate mismatches (substitutions + indels)
        edit_distance = result["editDistance"]
        
        # For primer detection, we want to know number of mismatches (substitutions)
        # but edlib gives edit distance (includes indels)
        # For now, use edit distance as mismatch approximation
        mismatches = edit_distance
        
        if mismatches <= max_mismatch_allowed:
            # Calculate identity
            identity = (primer_len - mismatches) / primer_len
            
            # Check minimum overlap
            if identity >= min_overlap_ratio and identity > best_score:
                best_score = identity
                read_start = search_offset + i
                read_end = read_start + primer_len
                
                # Get CIGAR for match details
                cigar = result.get("cigar", "")
                
                best_match = {
                    'read_start': read_start,
                    'read_end': read_end,
                    'identity': identity,
                    'mismatches': mismatches,
                    'cigar': cigar,
                    'primer_seq': primer_seq,
                    'matched_seq': subseq
                }
    
    return best_match


def validate_primers(
    read_seq: str,
    primer5_seq: str,
    primer3_seq: str,
    max_mismatch_rate: float = 0.15,
    min_overlap_ratio: float = 0.8
) -> Tuple[bool, Optional[Dict[str, Any]], bool, Optional[Dict[str, Any]]]:
    """
    Validate read contains both Primer5 and Primer3.
    
    Args:
        read_seq: Read sequence to validate
        primer5_seq: Primer5 reference sequence
        primer3_seq: Primer3 reference sequence (already in correct orientation)
        max_mismatch_rate: Maximum allowed mismatch rate
        min_overlap_ratio: Minimum overlap ratio with primer
        
    Returns:
        Tuple: (has_primer5, primer5_match_info, has_primer3, primer3_match_info)
    """
    # Detect Primer5 at 5' end
    p5_match = detect_primer_at_end(
        read_seq, primer5_seq, end='start',
        max_mismatch_rate=max_mismatch_rate,
        min_overlap_ratio=min_overlap_ratio
    )
    
    # Detect Primer3 at 3' end
    p3_match = detect_primer_at_end(
        read_seq, primer3_seq, end='end',
        max_mismatch_rate=max_mismatch_rate,
        min_overlap_ratio=min_overlap_ratio
    )
    
    return p5_match is not None, p5_match, p3_match is not None, p3_match


def create_primer_match(
    primer_name: str,
    match_info: Optional[Dict[str, Any]],
    is_valid: bool
) -> Optional[PrimerMatch]:
    """
    Create PrimerMatch object from match info.
    
    Args:
        primer_name: Name of primer ('Primer5' or 'Primer3')
        match_info: Match info dictionary from detect_primer_at_end
        is_valid: Whether primer is considered valid
        
    Returns:
        PrimerMatch object or None if match_info is None
    """
    if match_info is None:
        return None
    
    return PrimerMatch(
        primer_name=primer_name,
        read_start=match_info['read_start'],
        read_end=match_info['read_end'],
        identity=match_info['identity'],
        is_valid=is_valid
    )


def filter_reads_by_primers(
    reads: Dict[str, str],
    primer5_seq: str,
    primer3_seq: str,
    max_mismatch_rate: float = 0.15,
    min_overlap_ratio: float = 0.8
) -> Dict[str, Tuple[bool, Dict[str, Any], Dict[str, Any]]]:
    """
    Filter reads by primer validation.
    
    Args:
        reads: Dictionary of read_id -> read_seq
        primer5_seq: Primer5 sequence
        primer3_seq: Primer3 sequence
        max_mismatch_rate: Maximum allowed mismatch rate
        min_overlap_ratio: Minimum overlap ratio
        
    Returns:
        Dictionary of read_id -> (is_valid, p5_match_info, p3_match_info)
    """
    results = {}
    
    for read_id, read_seq in reads.items():
        has_p5, p5_match, has_p3, p3_match = validate_primers(
            read_seq, primer5_seq, primer3_seq,
            max_mismatch_rate=max_mismatch_rate,
            min_overlap_ratio=min_overlap_ratio
        )
        
        is_valid = has_p5 and has_p3
        results[read_id] = (is_valid, p5_match, p3_match)
    
    return results