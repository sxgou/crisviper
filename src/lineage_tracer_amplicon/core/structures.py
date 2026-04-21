"""
Core data structures for lineage tracing amplicon analysis.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Dict, Any


@dataclass
class Anchor:
    """Anchor definition."""
    name: str
    ref_start: int  # 0-based inclusive
    ref_end: int    # 0-based exclusive
    seq: str
    type: str       # "primer", "conserved", "pam", "prefix", "postfix", etc.


@dataclass
class AnchorMatch:
    """Anchor match result."""
    anchor: Anchor
    read_start: int
    read_end: int
    identity: float
    cigar: str      # CIGAR relative to anchor


@dataclass
class PrimerMatch:
    """Primer match information."""
    primer_name: str    # "Primer5" or "Primer3"
    read_start: int
    read_end: int
    identity: float
    is_valid: bool


@dataclass
class DeletionEvent:
    """Deletion event."""
    ref_start: int
    ref_end: int
    inserted_seq: str = ""      # empty for pure deletion
    is_mmej: bool = False
    mh_left: str = ""
    mh_right: str = ""


@dataclass
class InsertionEvent:
    """Insertion event (small insertion)."""
    ref_pos: int          # insertion point reference coordinate
    inserted_seq: str


@dataclass
class SNVEvent:
    """Point mutation event (only valid SNVs)."""
    ref_pos: int
    ref_base: str
    alt_base: str
    target_name: str     # target it belongs to


@dataclass
class Target:
    """Target definition."""
    name: str
    ref_start: int
    ref_end: int
    conserved_start: int
    conserved_end: int
    cutsite_start: int
    cutsite_end: int


@dataclass
class AmpliconStructure:
    """Amplicon structure definition."""
    reference: str
    features: List[Dict[str, Any]]
    anchors: List[Anchor] = field(default_factory=list)
    targets: List[Target] = field(default_factory=list)
    
    def __post_init__(self):
        """Parse features to extract anchors and targets."""
        self._parse_features()
    
    def _parse_features(self):
        """Parse features list to create anchors and targets."""
        for feature in self.features:
            feat_type = feature.get("type", "")
            name = feature.get("name", "")
            start = feature.get("start", 0)
            end = feature.get("end", 0)
            
            # Extract anchor based on type
            if feat_type in ["primer", "prefix", "postfix", "conserved", "pam"]:
                seq = self.reference[start:end]
                anchor = Anchor(
                    name=name,
                    ref_start=start,
                    ref_end=end,
                    seq=seq,
                    type=feat_type
                )
                self.anchors.append(anchor)
            
            # Extract target
            if feat_type == "target":
                conserved = feature.get("conserved", [start, start])
                cutsite = feature.get("cutsite", [start, end])
                target = Target(
                    name=name,
                    ref_start=start,
                    ref_end=end,
                    conserved_start=conserved[0],
                    conserved_end=conserved[1],
                    cutsite_start=cutsite[0],
                    cutsite_end=cutsite[1]
                )
                self.targets.append(target)


@dataclass
class ReadResult:
    """Complete read analysis result."""
    read_id: str
    is_valid: bool                     # passed primer validation
    primer5_match: Optional[PrimerMatch] = None
    primer3_match: Optional[PrimerMatch] = None
    anchor_matches: List[AnchorMatch] = field(default_factory=list)
    deletion_events: List[DeletionEvent] = field(default_factory=list)
    insertion_events: List[InsertionEvent] = field(default_factory=list)
    snv_events: List[SNVEvent] = field(default_factory=list)
    target_states: Dict[str, str] = field(default_factory=dict)  # target name -> state