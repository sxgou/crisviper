"""crisviper/models.py — Type-safe data models for lineage tracing analysis.

All data models that flow through the pipeline stages are defined here.
Each model is a dataclass with serialization support (to_dict/from_dict).
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import List, Dict, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from crisviper.config import AmpliconConfig, CutsiteRegion


# ═══════════════════════════════════════════════════════════════
# Mutation type enumeration
# ═══════════════════════════════════════════════════════════════

class MutationType(Enum):
    """Enumeration of mutation types detected from alignments."""
    SUBSTITUTION = "substitution"       # Point mutation (base replacement)
    DELETION = "deletion"               # Deletion (gap in query sequence)
    INSERTION = "insertion"             # Insertion (gap in reference sequence)
    INDEL = "indel"                     # Complex event with adjacent insertion+deletion


# ═══════════════════════════════════════════════════════════════
# Single mutation event
# ═══════════════════════════════════════════════════════════════

@dataclass
class MutationEvent:
    """A single mutation event extracted from an alignment.

    Represents one independent mutation (substitution, deletion, insertion,
    or combined indel) parsed from the pairwise alignment.
    One alignment result may contain multiple MutationEvents.

    Attributes:
        type: Mutation type (SUBSTITUTION, DELETION, INSERTION, or INDEL).
        ref_pos: Start position on the reference sequence (0-indexed).
        ref_base: Reference bases involved (for substitution/del).
        query_base: Query bases involved (for substitution/ins).
        length: Length of the mutation in bp.
        in_cutsite_window: Whether the mutation falls within a cutsite window.
        raw_ref_segment: Raw reference segment from the alignment.
        raw_query_segment: Raw query segment from the alignment.
        score: Confidence score for this event.
    """
    type: MutationType
    ref_pos: int
    ref_base: str = ""
    query_base: str = ""
    length: int = 1
    in_cutsite_window: bool = False
    raw_ref_segment: str = ""
    raw_query_segment: str = ""
    score: float = 0.0

    def to_dict(self) -> Dict:
        return {
            "type": self.type.value if hasattr(self.type, 'value') else str(self.type),
            "ref_pos": self.ref_pos,
            "ref_base": self.ref_base,
            "query_base": self.query_base,
            "length": self.length,
            "in_cutsite_window": self.in_cutsite_window,
            "score": self.score,
            "raw_ref_segment": self.raw_ref_segment,
            "raw_query_segment": self.raw_query_segment,
        }


# ═══════════════════════════════════════════════════════════════
# Input query record
# ═══════════════════════════════════════════════════════════════

@dataclass
class QueryRecord:
    """A single query sequence record with metadata.

    Represents one unique sequence read from FASTQ/TSV input,
    after deduplication by sequence. When keep_read_names is enabled,
    original_read_names stores the original FASTQ read identifiers
    that map to this deduplicated sequence.
    """
    readName: str
    cellBC: str = "unknown"
    UMI: str = "unknown"
    readCount: int = 1
    seq: str = ""
    original_read_names: List[str] = field(default_factory=list)


# ═══════════════════════════════════════════════════════════════
# Alignment statistics
# ═══════════════════════════════════════════════════════════════

@dataclass
class AlignmentStats:
    """Detailed alignment statistics for a single sequence.

    Tracks matches, mismatches, gap counts, gap block lengths,
    and derived metrics like similarity and identity.
    """
    matches: int = 0
    mismatches: int = 0
    gaps_in_ref: int = 0
    gaps_in_query: int = 0
    gap_blocks_ref: List[int] = field(default_factory=list)
    gap_blocks_query: List[int] = field(default_factory=list)
    avg_gap_len_ref: float = 0.0
    avg_gap_len_query: float = 0.0
    alignment_length: int = 0
    similarity: float = 0.0
    identity: float = 0.0
    score: float = 0.0

    @classmethod
    def from_dict(cls, d: Dict) -> "AlignmentStats":
        """Create from a dict (for backward compatibility)."""
        return cls(
            matches=d.get("matches", 0),
            mismatches=d.get("mismatches", 0),
            gaps_in_ref=d.get("gaps_in_ref", 0),
            gaps_in_query=d.get("gaps_in_query", 0),
            gap_blocks_ref=d.get("gap_blocks_ref", []),
            gap_blocks_query=d.get("gap_blocks_query", []),
            avg_gap_len_ref=d.get("avg_gap_len_ref", 0.0),
            avg_gap_len_query=d.get("avg_gap_len_query", 0.0),
            alignment_length=d.get("alignment_length", 0),
            similarity=d.get("similarity", 0.0),
            identity=d.get("identity", 0.0),
            score=d.get("score", 0.0),
        )

    def to_dict(self) -> Dict:
        """Convert to dict (for backward compatibility)."""
        return {
            "matches": self.matches,
            "mismatches": self.mismatches,
            "gaps_in_ref": self.gaps_in_ref,
            "gaps_in_query": self.gaps_in_query,
            "gap_blocks_ref": self.gap_blocks_ref,
            "gap_blocks_query": self.gap_blocks_query,
            "avg_gap_len_ref": self.avg_gap_len_ref,
            "avg_gap_len_query": self.avg_gap_len_query,
            "alignment_length": self.alignment_length,
            "similarity": self.similarity,
            "identity": self.identity,
            "score": self.score,
        }

    @property
    def has_indel(self) -> bool:
        """Whether this alignment contains insertions or deletions."""
        return self.gaps_in_ref > 0 or self.gaps_in_query > 0

    @property
    def has_mutation(self) -> bool:
        """Whether this alignment has any mutation (substitution or indel)."""
        return self.mismatches > 0 or self.has_indel


# ═══════════════════════════════════════════════════════════════
# Alignment result
# ═══════════════════════════════════════════════════════════════

@dataclass
class AlignmentResult:
    """Complete alignment result for a single query sequence.

    Includes the raw pairwise alignment, statistics, and extracted
    mutation events. Used as the primary data structure flowing
    through the pipeline.
    """
    query: QueryRecord           # Original query record
    success: bool = True         # Whether alignment succeeded
    score: float = 0.0           # Alignment score
    aligned_ref: str = ""        # Aligned reference sequence (with gaps)
    aligned_query: str = ""      # Aligned query sequence (with gaps)
    stats: Optional[AlignmentStats] = None  # Alignment statistics
    error: str = ""              # Error message (when success=False)
    mutations: List[MutationEvent] = field(default_factory=list)  # Extracted mutation events
    mode: str = "standard"       # Alignment mode (standard / lineage)
    failure_category: str = ""   # Failure category: "anchor", "noise", "alignment", "extraction"

    @classmethod
    def error_result(cls, query: QueryRecord, error_msg: str, category: str = "alignment") -> "AlignmentResult":
        """Create an error result for a failed alignment."""
        return cls(query=query, success=False, error=error_msg, failure_category=category)

    def to_dict(self) -> Dict:
        """Convert to dict (for backward compatibility with old output format)."""
        base = {
            "readName": self.query.readName,
            "cellBC": self.query.cellBC,
            "UMI": self.query.UMI,
            "readCount": self.query.readCount,
        }
        if self.query.original_read_names:
            base["original_read_names"] = self.query.original_read_names
        if not self.success or self.stats is None:
            base.update({
                "error": self.error,
                "score": None,
                "aligned_ref": None,
                "aligned_query": None,
                "stats": None,
                "mutations": [],
            })
        else:
            base.update({
                "score": self.score,
                "aligned_ref": self.aligned_ref,
                "aligned_query": self.aligned_query,
                "stats": self.stats.to_dict(),
                "mutations": [m.to_dict() for m in self.mutations],
            })
        return base


# ═══════════════════════════════════════════════════════════════
# Pipeline configuration
# ═══════════════════════════════════════════════════════════════

@dataclass
class PipelineConfig:
    """Centralized pipeline configuration with sensible defaults.

    All configurable pipeline parameters are collected here.
    Override any field to adjust pipeline behavior without modifying
    function code.
    """
    # ── Alignment scoring ──
    match_score: float = 2.0
    mismatch_penalty: float = -3.0
    gap_open: float = -2.0
    gap_extend: float = -0.1

    # ── Lineage tracer mode ──
    lineage_mode: bool = False
    # Gradient penalty: smoothly varying penalties centered on cutsites
    gradient_mode: bool = True      # Enable position-aware penalties in standard mode (requires cutsites)
    min_scale: float = 1.0          # Minimum penalty scale at cut center
    max_scale: float = 6.0          # Maximum penalty scale in conserved regions
    cutsite_edge_scale: float = 2.0  # Penalty scale at cutsite boundary
    gradient_radius: Optional[float] = None  # Gradient half-radius: None=auto(multi-cutsite)/30bp(single)
    mismatch_density_threshold: float = 0.34
    sub_window: int = 3

    # ── Gap exit penalty ──
    # gap_exit_strength: How strongly to penalize exiting a gap (<=0, 0=disabled).
    #   Automatically gradient-scaled by cutsite position for strongest suppression
    #   at cutsite centers. Default 0.0=disabled, -3.5 recommended (~-21.0 peak at center).
    gap_exit_strength: float = 0.0

    # ── Short-match discount ──
    # short_match_window: Threshold for short match regions (bp), 0=disabled.
    # short_match_discount: Match score discount factor (0~1), 1.0=no discount.
    #   E.g., window=3, discount=0.5 halves match_score for <=3bp match runs.
    short_match_window: int = 0
    short_match_discount: float = 1.0

    # ── Dense mismatch region penalty ──
    # dense_mismatch_window: Window size for dense mismatch detection (bp).
    # dense_mismatch_penalty: Extra penalty for dense mismatch regions (<=0, 0=disabled).
    #   Negative values bias the DP toward insertion paths in dense mismatch areas.
    dense_mismatch_window: int = 6
    dense_mismatch_penalty: float = 0.0

    # ── Homology repeat penalty (cross-target protection) ──
    # homology_window: Window size for homology detection (bp).
    # homology_penalty: <=0, subtracted from match_score at homologous positions, 0=disabled.
    #   Negative values make the DP less likely to match in repetitive ref regions.
    homology_window: int = 8
    homology_penalty: float = 0.0

    # ── Isolated base endpoint consolidation ──
    # Absorbs isolated single-base matches into adjacent gap endpoints.
    #   Works with gap_exit_strength: gap_exit_strength penalizes gap→M transitions,
    #   isolated_base_penalty penalizes when only 1bp matches after the transition.
    isolated_base_penalty: float = 0.0

    # ── Primer parameters ──
    primer5_len: int = 23
    primer3_len: int = 33
    primer5_threshold: int = 19
    primer3_threshold: int = 29

    # ── Allele filtering (inclusive thresholds: >=threshold passes) ──
    min_reads_sub: int = 5       # Minimum readCount for pure substitution alleles
    min_reads_indel: int = 0     # Minimum readCount for indel-containing alleles (0=no filter)

    # ── Background substitution correction ──
    correct_bg_sub: bool = True           # Enable background substitution correction
    keep_sub_indel_window: int = 3        # Regions to keep near indels (bp)

    # ── Multi-threading ──
    threads: int = 1
    chunk_size: int = 500

    # ── Reporting ──
    report_format: Optional[str] = None   # json / html
    allele_top_n: int = 50
    allele_window_start: int = 0
    allele_window_end: Optional[int] = None

    # ── Cutsite configuration ──
    cutsites_path: Optional[str] = None    # Path to JSON cutsite config file
    auto_detect_cutsites: bool = True

    # ── YAML config extension fields (target/amplicon structure) ──
    # Populated by cli.py when loading YAML configuration
    amplicon_config: Optional["AmpliconConfig"] = None   # Amplicon structure config
    explicit_cutsites: Optional[List["CutsiteRegion"]] = None  # Explicit cutsites from YAML

    # ── Denoising and allele calling ──
    denoise_enabled: bool = False          # Enable UMI/CB denoising
    call_alleles_enabled: bool = False     # Enable allele calling
    call_alleles_mode: str = "coarse"      # "coarse" or "exact"
    dominant_frac: float = 0.5             # Dominant allele fraction threshold


# ═══════════════════════════════════════════════════════════════
# Pipeline statistics
# ═══════════════════════════════════════════════════════════════

@dataclass
class PipelineStats:
    """Aggregate pipeline statistics for report generation."""
    total_queries: int = 0
    successful: int = 0
    failed: int = 0
    total_reads: int = 0
    mutated_sequences: int = 0
    unmutated_sequences: int = 0
    mutated_reads: int = 0
    n_anchor_failed: int = 0
    n_noise_filtered: int = 0

    @property
    def editing_efficiency_pct(self) -> float:
        """Editing efficiency as a percentage (sequence-level)."""
        if self.successful == 0:
            return 0.0
        return self.mutated_sequences / self.successful * 100

    @property
    def editing_efficiency_reads_pct(self) -> float:
        """Editing efficiency as a percentage (read-level, weighted by read count)."""
        if self.total_reads == 0:
            return 0.0
        return self.mutated_reads / self.total_reads * 100


@dataclass
class PipelineResult:
    """Final result from a complete pipeline run."""
    results: List[AlignmentResult]         # All alignment results
    config: PipelineConfig                 # Configuration used
    stats: PipelineStats                   # Pipeline statistics
    ref_length: int = 0                    # Reference sequence length
    mutation_type_counts: Dict = field(default_factory=dict)  # Mutation type counts
    total_mismatches: int = 0              # Total point mutations
    insertion_lengths: List[int] = field(default_factory=list)
    deletion_lengths: List[int] = field(default_factory=list)
    called_alleles: List = field(default_factory=list)  # CalledAllele list (optional)

    def get_successful(self) -> List[AlignmentResult]:
        """Get all successfully aligned results."""
        return [r for r in self.results if r.success]

    def get_failed(self) -> List[AlignmentResult]:
        """Get all failed alignment results."""
        return [r for r in self.results if not r.success]

    def get_mutated(self) -> List[AlignmentResult]:
        """Get all results containing mutations."""
        return [r for r in self.get_successful()
                if r.stats and r.stats.has_mutation]

    def get_unmutated(self) -> List[AlignmentResult]:
        """Get all results without mutations."""
        return [r for r in self.get_successful()
                if r.stats and not r.stats.has_mutation]
