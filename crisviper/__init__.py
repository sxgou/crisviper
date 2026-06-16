"""
crisviper — Lineage tracing analysis toolkit
=============================================
Core functionality for CRISPR-Cas9 lineage tracer amplicon analysis.

Module overview (in pipeline order):
  - io:           Data conversion (FASTQ → TSV/FASTA), single-cell 10x/InDrops parsing
  - alignment:    Core sequence alignment algorithms (Gotoh, position-aware DP)
  - lineage:      Lineage-tracer alignment mode with structure-aware penalties
  - mutation:     Mutation event extraction and classification from alignments
  - pipeline:     Full analysis pipeline orchestration and per-sequence workflow
  - reporting:    Output serialization JSON/TSV, summary reports, HTML/JSON report generation
  - plotting:     Visualization charts (allele heatmaps, indel distributions)
  - config:       Data structures for amplicon configuration and YAML loading
  - models:       Type-safe data models used across all pipeline stages
  - summary:      Summary statistics tables (allele frequency, per-target editing, filtering reasons)
  - denoiser:     UMI/CB denoising via directional adjacency clustering
  - caller:       Allele calling (coarse-grain and exact consensus)
  - threshold:    Statistical read-count threshold computation for UMI/CB filtering
  - metrics:      Diversity and heterogeneity metrics (Shannon entropy, effective alleles)
  - logging_config: Logging configuration framework
"""
from crisviper.config import CutsiteRegion, AmpliconConfig
from crisviper.models import (
    QueryRecord, AlignmentResult, AlignmentStats,
    PipelineConfig, PipelineResult, PipelineStats,
    MutationEvent, MutationType,
)
from crisviper.alignment import (
    affine_gap_alignment, calculate_alignment_stats,
    count_gap_blocks, affine_gap_alignment_position_aware,
)
from crisviper.lineage import (
    lineage_tracer_align, build_gradient_profiles,
    build_homology_penalty_profile, get_amplicon_structure,
)
from crisviper.mutation import (
    extract_mutations, classify_mutation_type,
    build_mutation_summary,
    annotate_mutation, annotate_mutations,
    classify_bp_event,
)
from crisviper.pipeline import Pipeline, align_single
from crisviper.io import (
    fastq_to_dataframe, fastq_to_fasta, fastq_to_fasta_from_rows,
    merge_paired_end, save_tsv,
    read_reference_fasta, read_queries_tsv, read_queries_fasta,
    parse_10x_fastq, parse_10x_provenance,
    filter_sc_cbs_and_umis, parse_indrops_provenance,
)
from crisviper.denoiser import (
    directional_adjacency_top_down_denoiser,
)
from crisviper.threshold import (
    compute_threshold,
)
from crisviper.caller import (
    call_alleles_coarse_grain, call_alleles_exact,
    CalledAllele, _event_structure,
)
from crisviper.reporting import (
    save_alignment_results, generate_report, save_text_report,
)
from crisviper.summary import save_summary_tables
from crisviper.metrics import (
    effective_alleles, diversity_index, alleles_per_cell,
    singletons_per_cell, carlin_potential,
)
from crisviper.logging_config import get_logger, setup_logging
