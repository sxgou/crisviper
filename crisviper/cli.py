"""crisviper CLI — Lineage tracing sequence analysis command-line tool.

Sequence alignment tool based on affine gap penalty algorithm,
supporting FASTQ conversion and parallel batch alignment.

Sub-commands:
  convert    Convert FASTQ files to TSV or FASTA format
  align      Parallel batch alignment of query sequences against a reference

Examples:
  $ crisviper convert fastq-to-tsv --fastq reads.fastq.gz --output reads.tsv
  $ crisviper align --reference ref.fasta --queries queries.tsv --output alignments.json
  $ crisviper align --reference ref.fasta --queries queries.tsv --output results/prefix --format all --report html
"""

import argparse
import sys
import os
import time
import multiprocessing as mp
# Prevent fork + NumPy thread conflicts (must set before importing crisviper)
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from typing import List, Dict

from crisviper import (
    PipelineConfig, QueryRecord,
    fastq_to_dataframe, fastq_to_fasta, fastq_to_fasta_from_rows,
    merge_paired_end, save_tsv,
    read_reference_fasta, read_queries_tsv, read_queries_fasta,
    Pipeline,
    save_alignment_results, save_summary_tables, generate_report,
    get_logger, setup_logging,
)
from crisviper.config import AmpliconConfig, load_yaml_config, cutsites_from_list

log = get_logger(__name__)

__version__ = "1.2.0"


def _queries_to_records(queries: List[Dict]) -> List[QueryRecord]:
    """Convert legacy dict-format queries to QueryRecord objects."""
    return [QueryRecord(**q) for q in queries]


def _log_timing(logger, label: str, t_start: float) -> float:
    """Log elapsed time since t_start (seconds), return current time."""
    elapsed = time.perf_counter() - t_start
    logger.info("  ⏱ %s: %.1fs", label, elapsed)
    return time.perf_counter()


# ── CLI Arg → PipelineConfig field mapping (used for YAML + CLI merging) ──
_CLI_TO_CONFIG = {
    'match_score': 'match_score',
    'mismatch_penalty': 'mismatch_penalty',
    'gap_open': 'gap_open',
    'gap_extend': 'gap_extend',
    'min_scale': 'min_scale',
    'max_scale': 'max_scale',
    'cutsite_edge_scale': 'cutsite_edge_scale',
    'gradient_radius': 'gradient_radius',
    'sub_window': 'sub_window',
    'mismatch_density_threshold': 'mismatch_density_threshold',
    'gap_exit_strength': 'gap_exit_strength',
    'short_match_window': 'short_match_window',
    'short_match_discount': 'short_match_discount',
    'dense_mismatch_window': 'dense_mismatch_window',
    'dense_mismatch_penalty': 'dense_mismatch_penalty',
    'homology_window': 'homology_window',
    'homology_penalty': 'homology_penalty',
    'isolated_base_penalty': 'isolated_base_penalty',
    'primer5_len': 'primer5_len',
    'primer3_len': 'primer3_len',
    'primer5_threshold': 'primer5_threshold',
    'primer3_threshold': 'primer3_threshold',
    'min_reads_sub': 'min_reads_sub',
    'min_reads_indel': 'min_reads_indel',
    'keep_sub_indel_window': 'keep_sub_indel_window',
    'threads': 'threads',
    'chunk_size': 'chunk_size',
    'cutsites': 'cutsites_path',
    'report': 'report_format',
    'allele_top_n': 'allele_top_n',
    'allele_window_start': 'allele_window_start',
    'allele_window_end': 'allele_window_end',
}

# Boolean CLI flags (must detect explicit presence via sys.argv)
_BOOL_CLI_FLAGS = {
    'lineage': 'lineage_mode',
    'correct_bg_sub': 'correct_bg_sub',
}


def _explicit_cli_flags() -> set:
    """Detect which CLI flags were explicitly set in sys.argv.

    Returns a set of flag names (with hyphens converted to underscores).
    Handles both --flag and --no-flag forms, as well as --flag=value.
    For --no-* flags, both the full name (no_xxx) and the base name (xxx)
    are added so that _BOOL_CLI_FLAGS matching succeeds.
    """
    explicit = set()
    for arg in sys.argv[1:]:
        # Handle --flag and --no-flag
        if arg.startswith('--'):
            parts = arg.split('=', 1)  # Handle --flag=value form
            name = parts[0].lstrip('-').replace('-', '_')
            if name:
                explicit.add(name)
                # For --no-* flags, also add the base name (strip no_ prefix)
                if name.startswith('no_'):
                    explicit.add(name[3:])
    return explicit


def _build_pipeline_config(args, yaml_data: dict) -> PipelineConfig:
    """Build PipelineConfig by merging YAML config + CLI arguments.

    Priority (highest to lowest):
    1. Explicit CLI arguments
    2. YAML config file pipeline: section values
    3. PipelineConfig class defaults
    """

    # Step 1: Start with YAML pipeline section as base
    yaml_pipeline = yaml_data.get('pipeline', {}) if yaml_data else {}
    config_dict = dict(yaml_pipeline)

    # Step 2: CLI explicit flags override YAML values
    explicit = _explicit_cli_flags()
    for cli_dest, cfg_field in _CLI_TO_CONFIG.items():
        val = getattr(args, cli_dest, None)
        if val is not None:
            config_dict[cfg_field] = val

    # Step 3: Boolean flags — only override when explicitly in argv
    for cli_dest, cfg_field in _BOOL_CLI_FLAGS.items():
        if cli_dest in explicit:
            config_dict[cfg_field] = getattr(args, cli_dest)

    # Step 4: YAML amplicon config → build AmpliconConfig
    if yaml_data and yaml_data.get('amplicon'):
        config_dict['amplicon_config'] = AmpliconConfig.from_dict(yaml_data['amplicon'])

    # Step 5: YAML explicit cutsites list
    if yaml_data and yaml_data.get('cutsites'):
        config_dict['explicit_cutsites'] = cutsites_from_list(yaml_data['cutsites'])

    return PipelineConfig(**config_dict)


def main():
    try:
        mp.set_start_method('fork')
    except RuntimeError:
        pass
    setup_logging(quiet=False)
    parser = argparse.ArgumentParser(
        description="crisviper — Lineage tracing sequence analysis CLI tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show debug information')
    parser.add_argument('-q', '--quiet', action='store_true', help='Show warnings and errors only')

    subparsers = parser.add_subparsers(dest="command", required=True, help="Sub-command")

    # Shared parent parser for paired-end merge parameters (DRY)
    merge_parent = argparse.ArgumentParser(add_help=False)
    merge_parent.add_argument("--min-overlap", type=int, default=10,
                              help="Minimum overlap length for paired-end merging (bp, default: 10)")
    merge_parent.add_argument("--max-mismatch-rate", type=int, default=20,
                              help="Max mismatch rate in overlap region for merging (%%, default: 20)")
    merge_parent.add_argument("--max-mismatch-diff", type=int, default=5,
                              help="Max absolute mismatches in overlap region (default: 5)")
    merge_parent.add_argument("--require-qual", type=int, default=15,
                              help="Minimum base quality score (Phred) for merging (default: 15)")

    # ── Sub-command: convert ─────────────────────────────────────
    convert_parser = subparsers.add_parser("convert", help="Convert FASTQ to TSV/FASTA format")
    convert_subparsers = convert_parser.add_subparsers(dest="convert_command", required=True, help="Conversion sub-command")

    # fastq-to-tsv
    tsv_parser = convert_subparsers.add_parser("fastq-to-tsv", parents=[merge_parent],
                                               help="Convert FASTQ to TSV format")
    tsv_parser.add_argument("--fastq", default=None, help="Input FASTQ file (single-end, supports .gz)")
    tsv_parser.add_argument("--fastq1", default=None, help="Paired-end R1 FASTQ file (use with --fastq2)")
    tsv_parser.add_argument("--fastq2", default=None, help="Paired-end R2 FASTQ file (use with --fastq1)")
    tsv_parser.add_argument("--output", required=True, help="Output TSV file path")
    tsv_parser.add_argument("--sample-name", default="sample", help="Sample name")
    tsv_parser.add_argument("--min-reads", type=int, default=1,
                            help="Minimum read count threshold (default: 1, no filtering)")

    # fastq-to-fasta
    fasta_parser = convert_subparsers.add_parser("fastq-to-fasta", parents=[merge_parent],
                                                 help="Convert FASTQ to FASTA format")
    fasta_parser.add_argument("--fastq", default=None, help="Input FASTQ file (single-end, supports .gz)")
    fasta_parser.add_argument("--fastq1", default=None, help="Paired-end R1 FASTQ file (use with --fastq2)")
    fasta_parser.add_argument("--fastq2", default=None, help="Paired-end R2 FASTQ file (use with --fastq1)")
    fasta_parser.add_argument("--output", required=True, help="Output FASTA file path")
    fasta_parser.add_argument("--sample-name", default="sample", help="Sample name")

    # ── Sub-command: align (parallel batch alignment) ─────────────
    align_parser = subparsers.add_parser("align", parents=[merge_parent],
                                         help="Parallel batch sequence alignment")
    align_parser.add_argument("--reference", required=True, help="Reference sequence FASTA file")
    align_parser.add_argument("--queries", default=None,
                              help="Query sequence file (TSV, FASTA, or FASTQ). For paired-end input, "
                                   "use --fastq1 and --fastq2 instead.")
    align_parser.add_argument("--fastq1", default=None,
                              help="R1 paired-end FASTQ file (requires --fastq2). Overrides --queries.")
    align_parser.add_argument("--fastq2", default=None,
                              help="R2 paired-end FASTQ file (requires --fastq1).")
    align_parser.add_argument("--output", required=True, help="Output file path (used as prefix with --format all)")
    align_parser.add_argument("--config", default=None,
                              help="YAML config file for target/amplicon structure and pipeline parameters (see crisviper_config.yaml)")
    align_parser.add_argument("--format", choices=["json", "tsv", "all"], default="json",
                              help="Output format: json (default), tsv, all")
    align_parser.add_argument("--sample-name", default="sample",
                              help="Sample name for read labeling (used when input is FASTQ)")

    # Alignment parameters (default=None means: use YAML config first, then PipelineConfig defaults)
    align_parser.add_argument("--match-score", type=float, default=None,
                              help="Score rewarded per aligned identical base in the Smith-Waterman DP matrix (default: 2.0). "
                                   "Purpose: Controls alignment sensitivity — higher values favor longer contiguous matches, "
                                   "lower values allow more gaps/mismatches. Principle: In affine gap alignment, the match "
                                   "score is added at every position where ref[i]==query[j]; the accumulated score determines "
                                   "the optimal path through the DP matrix. Use case: Keep 2.0 for standard Illumina CRISPR "
                                   "amplicon data. Reduce to 1.0-1.5 for noisy long-read (ONT/PacBio) data where substitutions "
                                   "are expected. Avoid values below 0.5 as alignment may fail to find legitimate matches.")
    align_parser.add_argument("--mismatch-penalty", type=float, default=None,
                              help="Penalty subtracted when ref[i] != query[j] in the DP matrix (default: -3.0). "
                                   "Purpose: Discourages accepting mismatched base pairs, controlling the trade-off "
                                   "between substitutions and indels. Principle: A more negative value forces the DP "
                                   "path to open gaps rather than accept mismatches; a less negative value tolerates "
                                   "substitutions, which may lead to fragmented alignments. Use case: -3.0 to -5.0 "
                                   "for high-quality Illumina data (low substitution error). -1.0 to -2.0 for "
                                   "error-prone long-read data or cross-species alignment. Highly conserved target: "
                                   "use -5.0 to -8.0 to aggressively penalize mismatches.")
    align_parser.add_argument("--gap-open", type=float, default=None,
                              help="Penalty for opening a new gap in the alignment (default: -2.0). "
                                   "Purpose: Controls the cost of initiating an insertion or deletion event. "
                                   "Principle: In affine gap alignment, gap_open is added once at the start of a "
                                   "gap, separate from gap_extend which accumulates per-base. A more negative value "
                                   "suppresses all indels; a less negative value allows more indels. Use case: -2.0 "
                                   "to -5.0 for CRISPR amplicons where indels are the signal of interest. For "
                                   "lineage-tracer mode where cutsite indels are expected, use -1.0 to -2.0. "
                                   "For highly conserved regions where indels are rare, use -6.0 to -10.0.")
    align_parser.add_argument("--gap-extend", type=float, default=None,
                              help="Penalty for each additional base in an existing gap (default: -0.1). "
                                   "Purpose: Controls the cost of extending an indel once opened — less negative "
                                   "than gap_open makes long indels relatively cheaper than multiple short indels. "
                                   "Principle: gap_extend << |gap_open| is the defining feature of affine gap "
                                   "alignment: it models biological reality where a single long indel is more "
                                   "likely than several short ones. Use case: -0.05 to -0.2 for most CRISPR data. "
                                   "Use -0.5 to -1.0 with short_match_window to aggressively split long indels "
                                   "into smaller events (useful for mutation counting). Keep gap_extend close to "
                                   "0 (e.g., -0.01) when only gap_open penalty is desired.")
    # Lineage tracer mode
    align_parser.add_argument("--lineage", action="store_true",
                              help="Enable lineage-tracer alignment mode: structure-aware dynamic penalties "
                                   "calibrated for CRISPR-edited target regions (default: disabled). Purpose: "
                                   "This mode applies position-aware penalty scaling centered on known cutsites, "
                                   "making the alignment far more sensitive to small indels near cutsites while "
                                   "remaining conservative in flanking regions. Principle: The DP scoring matrix "
                                   "is modulated by a smoothstep gradient profile — at cutsite centers penalties "
                                   "are reduced (min_scale) to encourage detecting genuine editing events, while "
                                   "primer regions receive maximum penalty protection (max_scale) to suppress "
                                   "artifacts. Use case: ALWAYS enable for lineage-tracing CRISPR amplicon data "
                                   "where cutsite positions are known or auto-detectable. Keep disabled for "
                                   "standard (non-CRISPR) alignment or unknown target structure.")
    align_parser.add_argument("--min-scale", type=float, default=None,
                              help="Minimum penalty scaling factor applied at the exact cutsite center (default: 1.0, "
                                   "range 0.1~1.0). Purpose: Controls how permissive the alignment is at cutsite "
                                   "centers — lower values make indels/substitutions easier to detect at CRISPR "
                                   "target sites. Principle: The smoothstep gradient interpolates between min_scale "
                                   "(at cutsite center) and max_scale (in conserved regions). At min_scale=0.3, "
                                   "all penalty magnitudes (gap_open, mismatch, etc.) are multiplied by 0.3, "
                                   "making it 70%% easier to open a gap at the cutsite. Use case: 0.1-0.3 for "
                                   "high-sensitivity lineage tracing (detect rare editing events). 0.5-0.8 for "
                                   "balanced detection. 1.0 effectively disables the gradient (constant scaling). "
                                   "Must be < max_scale for the gradient to have any effect.")
    align_parser.add_argument("--max-scale", type=float, default=None,
                              help="Maximum penalty scaling factor applied in conserved flanking and primer regions "
                                   "(default: 6.0, range 1.0~10.0). Purpose: Suppresses false-positive indels in "
                                   "regions that should remain intact — primers, barcodes, and conserved backbone. "
                                   "Principle: At max_scale=6.0, all penalty magnitudes are multiplied by 6x, "
                                   "making it 6× harder to open gaps in these protected regions versus unscaled "
                                   "alignment. Use case: 4.0-8.0 for standard CRISPR amplicons. Higher values "
                                   "(8.0-10.0) for targets with critical primer or barcode regions where ANY "
                                   "indel would indicate an artifact. Must be > min_scale.")
    align_parser.add_argument("--cutsite-edge-scale", type=float, default=None,
                              help="Penalty scaling factor at the boundary between cutsite window and conserved "
                                   "region (default: 2.0). Purpose: Controls the sharpness of the transition zone "
                                   "at cutsite edges — a smooth ramp prevents edge artifacts where small indels "
                                   "at boundary positions can be misclassified. Principle: The smoothstep profile "
                                   "uses three anchor points: min_scale (center), cutsite_edge_scale (boundary), "
                                   "and max_scale (outside). edge_scale defines the gradient value exactly at the "
                                   "cutsite boundary ±gradient_radius from center. Use case: 1.5-2.5 for typical "
                                   "CRISPR cutsites. Higher values (3.0-4.0) for sharper transitions (useful when "
                                   "cutsites are close together and windows might overlap). Lower values (1.0-1.5) "
                                   "for gradual transitions in dense cutsite arrays.")
    align_parser.add_argument("--gradient-radius", type=float, default=None,
                              help="Half-radius of the smoothstep penalty gradient in base pairs (default: auto). "
                                   "Purpose: Defines how far from each cutsite center the penalty modulation extends "
                                   "— larger radii widen the permissive zone around cutsites. Principle: The gradient "
                                   "transitions from min_scale (at center) through cutsite_edge_scale (at radius bp) "
                                   "to max_scale (at 2× radius). For multi-cutsite amplicons, radius auto-computes "
                                   "to half the inter-cutsite distance so adjacent gradients merge smoothly. For "
                                   "single cutsites, defaults to 30 bp. Use case: Auto (recommended) for most "
                                   "targets. Manual override to 20-40 bp for single-cutsite targets. 10-15 bp for "
                                   "dense cutsite arrays (<50 bp spacing). 50-80 bp for targets with large "
                                   "editable regions.")
    align_parser.add_argument("--sub-window", type=int, default=None,
                              help="Cutsite-adjacent retention window in base pairs for background substitution "
                                   "correction and mutation annotation (default: 3). Purpose: When background "
                                   "substitution correction is enabled (--correct-bg-sub), substitutions within "
                                   "this distance from any indel are preserved as genuine — others may be filtered "
                                   "as background noise. Also controls the window used to annotate whether a "
                                   "mutation falls 'inside' a cutsite region. Principle: Real Cas9-induced DSB "
                                   "repair often produces compound events (e.g., 1 bp substitution adjacent to a "
                                   "2 bp deletion). The sub_window keeps these linked events intact. Use case: "
                                   "3 (default) for standard CRISPR amplicons. 0-1 to aggressively separate "
                                   "indels from substitutions. 5-8 for targets with frequent compound indel+sub "
                                   "events. Only applies when correct_bg_sub is enabled.")
    align_parser.add_argument("--mismatch-density-threshold", type=float, default=None,
                              help="Density threshold for detecting dense mismatch regions in the DP alignment "
                                   "(default: 0.34, range 0~1). Purpose: When a local region exceeds this mismatch "
                                   "fraction, the DP switches to insertion-favoring mode (lineage mode only), "
                                   "treating the region as a likely Cas9-off-target or sequencing error cluster "
                                   "rather than real biological variation. Principle: The algorithm slides a window "
                                   "across the reference; if the fraction of mismatches > threshold, the insertion "
                                   "path is preferred over substitution. Use case: 0.34 (default) for typical "
                                   "CRISPR amplicon data. 0.2-0.3 for high-stringency (fewer INDEL calls but "
                                   "higher confidence). 0.4-0.5 for low-stringency (more INDEL calls, may include "
                                   "noise). Only applies in lineage or gradient mode.")
    align_parser.add_argument("--cutsites", default=None,
                              help="Path to a JSON file specifying cutsite positions for gradient penalty modulation "
                                   "(optional). Purpose: When provided, the gradient penalty profile centers on "
                                   "these explicitly defined positions rather than using auto-detection from the "
                                    "reference sequence. Principle: JSON format: [{\"name\": \"T1\", \"start\": 42, \"end\": 47}, ...] "
                                    "where start/end are the cutsite boundaries (0-indexed) on the reference. "
                                   "Unless cutsites are provided via JSON, --config YAML, or auto-detection, "
                                   "gradient mode will have no effect. Use case: Explicit cutsite config for "
                                   "multi-target amplicons, non-standard Cas enzymes, or validation against "
                                   "known cutsite positions. Overrides cutsite auto-detection.")

    align_parser.add_argument("--gap-exit-strength", type=float, default=None,
                              help="Gap exit suppression: penalty applied when the DP transitions from a gap back "
                                   "to match state (<=0, 0=disabled, default: 0). Purpose: Suppresses spurious "
                                   "single-base matches inside long indels by penalizing the gap→match transition. "
                                   "This prevents gap regions from being fragmented by isolated matches. Principle: "
                                   "The strength value is scaled by the position-aware gradient (same profile as "
                                   "other penalties), so suppression is strongest at cutsite centers. E.g., "
                                   "gap_exit_strength=-3.5 applied at a cutsite center with min_scale=0.3 produces "
                                   "an effective penalty of -3.5×0.3≈-1.05 at center, but -3.5×6.0≈-21.0 in "
                                   "conserved regions. Use case: 0 (disabled) for standard alignment. -2.0 to -5.0 "
                                   "for lineage mode to clean up gap endpoint fragmentation. -5.0 to -10.0 for "
                                   "aggressive consolidation (may merge distinct nearby indels). Recommended "
                                   "starting point: -3.5.")
    align_parser.add_argument("--short-match-window", type=int, default=None,
                              help="Short match region threshold in bp (default: 0, disabled; recommended: 3-5). "
                                   "Purpose: Detects and penalizes suspiciously short match runs flanked by gaps — "
                                   "a common pattern where an isolated 1-3 bp match breaks a long indel into "
                                   "separate DEL and INS events. Merging them into a single INDEL gives a more "
                                   "biologically accurate representation. Principle: When match runs shorter than "
                                   "this threshold are found between gap regions, their match score is discounted "
                                   "(see --short-match-discount), making the DP prefer a single continuous gap "
                                   "over match+gap+match. Use case: 0 (disabled) for standard alignment. 3-5 for "
                                   "lineage mode to consolidate fragmented indels. Higher values (6-10) for "
                                   "aggressive consolidation (may merge truly separate events).")
    align_parser.add_argument("--short-match-discount", type=float, default=None,
                              help="Discount factor applied to match_score for short match runs (0~1, 1.0=no "
                                   "discount, recommended: 0.5). Purpose: Reduces the match reward for short "
                                   "isolated match segments inside gaps, making the DP less likely to fragment "
                                   "indels. Principle: When a match run is shorter than short_match_window, the "
                                   "match score at each position is multiplied by this factor before being added "
                                   "to the DP matrix. At 0.5, a match_score=2.0 becomes effectively 1.0, "
                                   "reducing the incentive to split a gap. Use case: 0.3-0.5 for lineage mode "
                                   "with gap consolidation enabled. 1.0 (no discount) for standard alignment. "
                                   "Only effective when short_match_window > 0.")
    align_parser.add_argument("--dense-mismatch-window", type=int, default=None,
                              help="Sliding window size (bp) for dense mismatch detection (default: 6). "
                                   "Purpose: Defines the local window over which mismatch density is computed. "
                                   "When mismatches within this window exceed the density threshold, the DP "
                                   "prefers insertion paths over substitution accumulation. Principle: A wider "
                                   "window smooths out local noise but may miss short high-density regions; a "
                                   "narrower window is more sensitive to local clusters but may trigger "
                                   "false-positive insertion paths from random mismatches. Use case: 6 (default) "
                                   "for standard amplicons. 10-15 for long-read data with higher error rates "
                                   "(smoother detection). 3-4 for short amplicons (<200 bp). Only applies when "
                                   "dense_mismatch_penalty is enabled (< 0).")
    align_parser.add_argument("--dense-mismatch-penalty", type=float, default=None,
                              help="Extra penalty subtracted per mismatch in high-density mismatch regions "
                                   "(<=0, 0=disabled, recommended: -2.0). Purpose: When local mismatch density "
                                   "exceeds the threshold, this extra penalty is applied to each mismatch within "
                                   "the window, making the DP strongly prefer opening an insertion rather than "
                                   "accepting many clustered mismatches. This models Cas9 off-target or error "
                                   "clusters as indels rather than multiple substitutions. Principle: The penalty "
                                   "is added on top of mismatch_penalty — so at -2.0, total mismatch cost becomes "
                                   "-5.0 (with default -3.0 mismatch_penalty) in dense regions. Use case: 0 "
                                   "(disabled) for standard alignment. -1.0 to -3.0 for lineage mode to absorb "
                                   "dense substitution clusters into indels. -3.0 to -5.0 for aggressive "
                                   "consolidation.")
    align_parser.add_argument("--homology-window", type=int, default=None,
                              help="Window size (bp) for homology/repeat region detection (default: 8). "
                                   "Purpose: Determines the local context size used to identify repetitive "
                                   "or self-homologous reference regions. When the surrounding sequence is "
                                   "similar to the current position, homology_penalty is applied. This prevents "
                                   "the DP from aligning mismatched bases to repetitive regions by exploiting "
                                   "sequence similarity. Principle: The algorithm compares each reference "
                                   "position's ±window/2 context to detect local self-similarity. A larger "
                                   "window detects longer-range homology but may miss short repeats. Use case: "
                                   "8 (default) for typical amplicons. 12-16 for targets with long repeats "
                                   "(e.g., microsatellites). 4-6 for short targets. Only applies when "
                                   "homology_penalty < 0.")
    align_parser.add_argument("--homology-penalty", type=float, default=None,
                              help="Penalty subtracted from match_score at positions within homologous/repetitive "
                                   "reference regions (<=0, 0=disabled, recommended: -1.0). Purpose: Reduces the "
                                   "match score at repetitive positions, preventing the DP from falsely aligning "
                                   "a query to the wrong copy of a repeat. This is critical for multi-cutsite "
                                   "amplicons where cutsite-flanking sequences may share homology. Principle: At "
                                   "a position detected as homologous, match_score is reduced by this amount, "
                                   "making a gap more competitive versus a weak match. Use case: 0 (disabled) "
                                   "for single-cutsite amplicons. -0.5 to -2.0 for multi-cutsite targets with "
                                   "similar flanking sequences. -2.0 to -5.0 for highly repetitive targets.")
    align_parser.add_argument("--isolated-base-penalty", type=float, default=None,
                              help="Extra penalty for isolated single-base matches adjacent to gap endpoints "
                                   "(<=0, 0=disabled, recommended: -2.0). Purpose: Absorbs isolated single-base "
                                   "matches into adjacent gap regions by penalizing standalone matches that are "
                                   "flanked by gaps on both sides. These isolated matches are often alignment "
                                   "artifacts rather than real biological insertions. Principle: When a match "
                                   "run of exactly 1 bp is found between two gap regions, this penalty is "
                                   "applied to the match score, making the DP prefer extending the gap. "
                                   "Works synergistically with gap_exit_strength: gap_exit_strength penalizes "
                                   "the gap→M transition; isolated_base_penalty penalizes the actual match "
                                   "if it was too short. Use case: 0 (disabled) for standard alignment. -1.0 "
                                   "to -3.0 for lineage mode to clean up isolated bases. -5.0 for aggressive "
                                   "gap consolidation.")

    # Parallelism parameters
    align_parser.add_argument("--threads", "-t", type=int, default=None,
                              help="Number of parallel worker processes for batch alignment (default: 1, single-threaded). "
                                   "Purpose: Enables multi-CPU parallel processing by splitting query sequences into "
                                   "chunks and processing them across worker processes. Each worker independently "
                                   "aligns its chunk against the reference. Principle: Uses Python multiprocessing "
                                   "with fork start method — each worker loads its own copy of the reference and "
                                   "scoring matrices. Speedup is near-linear for large batch sizes, limited by "
                                   "memory bandwidth. Use case: Set to the number of physical CPU cores (not "
                                   "hyperthreads) for batch sizes > 500 sequences. For < 100 sequences, "
                                   "single-threaded is faster due to fork overhead. Max 12 recommended (system "
                                   "stability limit). Set to 1 for debugging (easier stack traces).")
    align_parser.add_argument("--chunk-size", type=int, default=None,
                              help="Number of query sequences assigned to each worker process per chunk "
                                   "(PipelineConfig default: 500; set to 0 for auto-compute as total_queries / (threads × 3)). Purpose: Controls "
                                   "the granularity of parallel work distribution. Smaller chunks improve load "
                                   "balancing (workers finish at similar times) but increase inter-process "
                                   "communication overhead. Principle: Each worker processes one chunk at a time; "
                                   "the main process distributes chunks via a multiprocessing queue. With ~3× "
                                   "as many chunks as workers, load imbalance from variable-length sequences "
                                   "is minimized. Use case: 500 (default) for batch sizes up to 5000 sequences. 50-200 for "
                                   "short sequences (< 200 bp). 200-500 for long sequences (> 500 bp). "
                                   "Set explicitly when auto-computed chunk counts are suboptimal.")

    # Primer parameters
    align_parser.add_argument("--primer5-len", type=int, default=None,
                              help="Length in bp of the 5' (forward) primer region (default: 23). "
                                   "Purpose: Defines the expected primer binding region at the 5' end of the "
                                   "amplicon. This region is excluded from mutation detection and is protected "
                                   "by max_scale gradient penalty. Principle: The primer5_len bases at the "
                                   "reference 5' end are treated as the primer annealing region during "
                                   "alignment. They must match above the primer5_threshold for a successful "
                                   "alignment. Use case: 23 for standard CRISPR amplicons. Adjust to match "
                                   "your experimental primer design. Must match the actual primer length "
                                   "used in library preparation.")
    align_parser.add_argument("--primer3-len", type=int, default=None,
                              help="Length in bp of the 3' (reverse) primer region (default: 33). "
                                   "Purpose: Defines the expected primer binding region at the 3' end of the "
                                   "amplicon. Like primer5_len, this region is excluded from mutation detection "
                                   "and gradient-modulated. Principle: The primer3_len bases at the reference "
                                   "3' end are checked against primer3_threshold for alignment success. "
                                   "Use case: 33 for standard CRISPR amplicons. Match your experimental "
                                   "reverse primer design. Adjust for different primer lengths.")
    align_parser.add_argument("--primer5-threshold", type=int, default=None,
                              help="Minimum number of matching bases required in the 5' primer region for "
                                   "alignment to succeed (default: 19). Purpose: Quality check — sequences "
                                   "with too many mismatches in the primer region likely represent off-target "
                                   "amplification or primer-dimer artifacts and should be rejected. Principle: "
                                   "If the 5' primer matches fewer than this threshold, the alignment is "
                                   "marked as failed (anchor failure). Use case: 19 for standard 23 bp primer "
                                   "(≈80%% match required). Lower threshold (14-16) for high-diversity targets "
                                   "or degenerate primers. Higher threshold (21-22) for stringent filtering.")
    align_parser.add_argument("--primer3-threshold", type=int, default=None,
                              help="Minimum number of matching bases required in the 3' primer region for "
                                   "alignment to succeed (default: 29). Purpose: Same as primer5_threshold but "
                                   "for the 3' end — serves as the second anchor quality check. Principle: "
                                   "Both primer5 and primer3 must pass their respective thresholds for a "
                                   "successful alignment. Use case: 29 for standard 33 bp primer (≈88%% "
                                   "match required). Lower for degenerate primers. Higher for stringent "
                                   "quality filtering.")

    # Allele filtering parameters
    align_parser.add_argument("--min-reads", type=int, default=1,
                              help="Input-side minimum read count threshold (default: 1, no filtering, used for pre-filtering)")
    align_parser.add_argument("--min-reads-sub", type=int, default=None,
                              help="Minimum read count (inclusive, >=threshold passes) for pure-substitution alleles "
                                   "(default: 5). Purpose: Filters out low-abundance substitution-only alleles "
                                   "that are likely sequencing errors or PCR artifacts rather than genuine "
                                   "biological variation. Principle: Pure substitution alleles (no indels) at "
                                   "low read counts are more likely PCR/sequencing errors than true mutations. "
                                   "A higher threshold removes more noise but may discard rare real variants. "
                                   "Use case: 5 for standard Illumina data (good base quality). 10-20 for "
                                   "noisy data or when high specificity is needed. 0-2 for maximum sensitivity "
                                   "(e.g., rare clone detection). Only applies when the allele has only "
                                   "substitutions (no indels).")
    align_parser.add_argument("--min-reads-indel", type=int, default=None,
                              help="Minimum read count (inclusive, >=threshold passes) for alleles containing indels (default: 0, "
                                   "no filter). Purpose: Indel alleles are the primary signal of Cas9 editing "
                                   "activity and are typically trusted even at low read counts. However, in "
                                   "noisy datasets, a low threshold can help filter spurious indel artifacts. "
                                   "Principle: Unlike substitutions, indels in CRISPR amplicons are expected "
                                   "to be genuine editing outcomes, so the default allows all indel alleles "
                                   "through. Use case: 0 (no filter) for standard lineage tracing. 2-5 for "
                                   "noisy long-read data where indels may include alignment errors. Keep at "
                                   "0 when detecting rare editing events.")

    # Background substitution correction parameters
    align_parser.add_argument("--correct-bg-sub", action=argparse.BooleanOptionalAction,
                              default=True, help="Enable background substitution correction (default: enabled). "
                                   "Purpose: Automatically identifies and filters background-level substitution "
                                   "noise (e.g., from PCR errors or sequencing mistakes) while preserving "
                                   "genuine substitutions near indels. Principle: Substitutions are classified "
                                   "as 'background' when they occur in isolation (no nearby indel) and at "
                                   "frequencies consistent with the global error profile. Substitutions within "
                                   "--keep-sub-indel-window of an indel are always preserved as they may be "
                                   "part of a compound editing event. Use case: Enabled (default) for standard "
                                   "CRISPR amplicon data where PCR errors are expected. Disable when analyzing "
                                   "samples with known high mutation rates or when every substitution must "
                                   "be reported regardless of context.")
    align_parser.add_argument("--keep-sub-indel-window", type=int, default=None,
                              help="Retention window in bp around each indel for preserving substitutions during "
                                   "background correction (default: 3). Purpose: Substitutions within this distance "
                                   "from an indel are preserved as potentially meaningful (compound editing "
                                   "events), while isolated substitutions further away may be filtered as "
                                   "background noise. Principle: Cas9-mediated repair often produces complex "
                                   "alleles with both indels and nearby substitutions. This window protects "
                                   "those linked events from background correction. Use case: 3 (default) for "
                                   "standard CRISPR editing. 5-10 for targets with known high compound event "
                                   "rates. 0 to disable protection (filter all isolated substitutions). "
                                   "Only applies when --correct-bg-sub is enabled.")

    # Report parameters
    align_parser.add_argument("--report", choices=["json", "html"], default=None,
                              help="Generate an analysis report after alignment completes (optional). "
                                   "Purpose: Produces a comprehensive visual or structured report containing "
                                   "QC metrics (read counts, alignment rates), mutation type distributions, "
                                   "indel length histograms, allele frequency heatmaps, and per-sequence "
                                   "mutation summaries. Principle: The report is generated from the aggregated "
                                   "alignment results and provides an at-a-glance overview of the experiment "
                                   "quality and editing outcomes. Use case: 'html' for a self-contained browser-"
                                   "viewable report with interactive visualizations. 'json' for programmatic "
                                   "access to aggregated report data. Omit to skip report generation (faster "
                                   "for batch processing pipelines). Requires --report-output or defaults "
                                   "to output path + '_report' suffix.")
    align_parser.add_argument("--report-output", default=None,
                              help="Output file path for the analysis report (optional). When omitted, the report "
                                   "path is auto-inferred from --output by stripping the extension and appending "
                                   "'_report'. Purpose: Specifies where the generated report file(s) will be "
                                   "written. For HTML reports, a .html extension is appended. For JSON reports, "
                                   "a .json extension is used. Use case: Provide a custom path to organize "
                                   "reports in a specific directory or naming convention. Auto-inferred path "
                                   "works well for most workflows.")
    align_parser.add_argument("--allele-window-start", type=int, default=None,
                              help="Start position (0-indexed) for the allele display window in reports "
                                   "(default: 0 = beginning of reference). Purpose: Restricts the allele "
                                   "heatmap and sequence display to a specific region of the reference, "
                                   "focusing attention on the target area. Principle: Full-length alleles "
                                   "include primer regions which are typically invariant. Truncating to "
                                   "the internal region produces cleaner visualizations. Use case: Set "
                                   "to primer5_len (e.g., 23) to exclude the 5' primer from allele displays. "
                                   "Use with --allele-window-end to zoom into a specific cutsite region. "
                                   "Leave at default (0) to show the full amplicon.")
    align_parser.add_argument("--allele-window-end", type=int, default=None,
                               help="End position (inclusive, 0-indexed) for the allele display window in reports "
                                   "(default: end of reference). Purpose: Defines the end boundary of the "
                                   "display window — see --allele-window-start for details. Use case: Set "
                                   "to len(ref) - primer3_len to exclude the 3' primer region. Set to a "
                                   "cutsite position + padding to zoom into a specific target.")
    align_parser.add_argument("--allele-top-n", type=int, default=None,
                              help="Maximum number of most abundant alleles to display in the allele heatmap "
                                   "report (default: 50). Purpose: Limits the heatmap to the top N most "
                                   "frequent alleles, avoiding visual clutter from rare alleles that may "
                                   "be noise. Principle: Alleles are sorted by read count (descending), "
                                   "and only the top N are included in the heatmap visualization. Use case: "
                                   "50 for standard reports with good coverage. 100-200 for comprehensive "
                                   "views of complex editing outcomes. 10-20 for quick overviews or when "
                                   "only dominant clones are of interest.")
    align_parser.add_argument("--read-to-allele", action="store_true",
                              help="Output read-to-allele mapping table to output folder. Requires FASTQ input (--queries "
                                   "or --fastq1/--fastq2). Only supported with FASTQ input; TSV and FASTA inputs "
                                   "will ignore this option with a warning.")

    args = parser.parse_args()
    setup_logging(verbose=args.verbose, quiet=args.quiet)

    # ═══════════════════════════════════════════════════════════════
    # Execute command
    # ═══════════════════════════════════════════════════════════════
    if args.command == "convert":
        if args.convert_command == "fastq-to-tsv":
            # Paired-end mode
            if args.fastq1 and args.fastq2:
                rows = merge_paired_end(
                    args.fastq1, args.fastq2,
                    min_overlap=args.min_overlap,
                    max_mismatch_rate=args.max_mismatch_rate,
                    max_mismatch_diff=args.max_mismatch_diff,
                    require_qual=args.require_qual,
                    sample_name=args.sample_name,
                )
            elif args.fastq:
                rows = fastq_to_dataframe(args.fastq, args.sample_name)
            else:
                log.error("Specify --fastq (single-end) or --fastq1 + --fastq2 (paired-end)")
                sys.exit(1)

            if args.min_reads > 1:
                before = len(rows)
                rows = [r for r in rows if r['readCount'] >= args.min_reads]
                log.info("Allele filtering (min_reads>=%d): %d -> %d", args.min_reads, before, len(rows))
            save_tsv(rows, args.output)

        elif args.convert_command == "fastq-to-fasta":
            if args.fastq1 and args.fastq2:
                rows = merge_paired_end(
                    args.fastq1, args.fastq2,
                    min_overlap=args.min_overlap,
                    max_mismatch_rate=args.max_mismatch_rate,
                    max_mismatch_diff=args.max_mismatch_diff,
                    require_qual=args.require_qual,
                    sample_name=args.sample_name,
                )
                fastq_to_fasta_from_rows(rows, args.output)
            elif args.fastq:
                fastq_to_fasta(args.fastq, args.output, args.sample_name)
            else:
                log.error("Specify --fastq (single-end) or --fastq1 + --fastq2 (paired-end)")
                sys.exit(1)

    elif args.command == "align":
        t0 = time.perf_counter()

        # Read reference sequence
        log.info("Reading reference sequence: %s", args.reference)
        ref_seq = read_reference_fasta(args.reference)
        log.info("Reference sequence length: %d bp", len(ref_seq))
        t = _log_timing(log, "Read reference", t0)

        # Read query sequences
        read_to_allele_active = args.read_to_allele
        if args.fastq1 and args.fastq2:
            # Paired-end FASTQ input
            log.info("Processing paired-end FASTQ: %s + %s", args.fastq1, args.fastq2)
            queries = merge_paired_end(
                args.fastq1, args.fastq2,
                min_overlap=args.min_overlap,
                max_mismatch_rate=args.max_mismatch_rate,
                max_mismatch_diff=args.max_mismatch_diff,
                require_qual=args.require_qual,
                sample_name=args.sample_name,
                keep_read_names=read_to_allele_active,
            )
        elif args.fastq1 or args.fastq2:
            log.error("--fastq1 and --fastq2 must be used together")
            sys.exit(1)
        elif args.queries:
            queries_path = args.queries
            if queries_path.endswith('.tsv'):
                queries = read_queries_tsv(queries_path)
                if read_to_allele_active:
                    log.warning("--read-to-allele is not supported with TSV input. "
                                "Use FASTQ input (single or paired-end) to generate read-to-allele mapping.")
                    read_to_allele_active = False
            elif queries_path.endswith(('.fasta', '.fa', '.fas')):
                queries = read_queries_fasta(queries_path)
                if read_to_allele_active:
                    log.warning("--read-to-allele is not supported with FASTA input. "
                                "Use FASTQ input (single or paired-end) to generate read-to-allele mapping.")
                    read_to_allele_active = False
            elif queries_path.endswith(('.fq', '.fastq', '.fq.gz', '.fastq.gz')):
                queries = fastq_to_dataframe(queries_path, args.sample_name,
                                              keep_read_names=read_to_allele_active)
            else:
                log.error("Unsupported query file format: %s", queries_path)
                sys.exit(1)
        else:
            log.error("Specify --queries (single-end) or --fastq1 + --fastq2 (paired-end) for input")
            sys.exit(1)
        t = _log_timing(log, "Read queries", t)

        # Convert to type-safe QueryRecord objects
        query_records = _queries_to_records(queries)

        # Allele quality pre-filtering (before conversion)
        if args.min_reads > 1:
            before = len(query_records)
            query_records = [q for q in query_records if q.readCount >= args.min_reads]
            log.info("Allele filtering (min_reads>=%d): %d -> %d",
                     args.min_reads, before, len(query_records))

        # ── Load YAML configuration ──
        yaml_data = {}
        if args.config:
            log.info("Loading config: %s", args.config)
            yaml_data = load_yaml_config(args.config)

        # ── Build PipelineConfig (YAML + CLI merge) ──
        config = _build_pipeline_config(args, yaml_data)

        # ── Compute chunk size (backward-compatible logic) ──
        total_queries = len(query_records)
        if total_queries == 0:
            log.warning("No query sequences to align (empty input)")
        if not config.chunk_size and total_queries > 0:
            target_chunks = max(config.threads * 3, 12)
            config.chunk_size = max(100, total_queries // target_chunks)
        elif not config.chunk_size:
            config.chunk_size = 100
        log.info("  Chunk size: %d seqs/chunk (%d chunks total)", config.chunk_size,
                 (total_queries + config.chunk_size - 1) // config.chunk_size if total_queries else 0)

        # Display configuration info
        radius_info = "auto" if config.gradient_radius is None else f"{config.gradient_radius}bp"
        log.info("  Gradient penalties: min_scale=%s, max_scale=%s, edge_scale=%s, radius=%s",
                 config.min_scale, config.max_scale, config.cutsite_edge_scale, radius_info)
        if config.lineage_mode:
            log.info("  sub_window: cutsite+/-%s bp, mismatch_density_threshold: %s",
                     config.sub_window, config.mismatch_density_threshold)

        if config.threads > 1:
            log.info("Starting parallel batch alignment (%d threads)...", config.threads)
        else:
            log.info("Starting batch alignment (single-threaded)...")
        # ── Run pipeline ──
        pipeline = Pipeline(config=config, ref_seq=ref_seq)
        pipeline_result = pipeline.run(query_records)
        t = _log_timing(log, "Alignment (parallel)", t)

        # ── Convert to compatible output format ──
        output_results = [r.to_dict() for r in pipeline_result.results]
        t = _log_timing(log, "Results to dict", t)

        # ── Save results ──
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        save_alignment_results(output_results, args.output, args.format)
        t = _log_timing(log, "Save TSV results", t)

        # ── Save summary tables ──
        output_dir = os.path.dirname(os.path.abspath(args.output))
        summary_data = save_summary_tables(output_results, output_dir, ref_seq=ref_seq,
                                            cutsites=_get_display_cutsites(ref_seq, config.amplicon_config),
                                            read_to_allele_path=os.path.join(output_dir, "read_to_allele.tsv") if read_to_allele_active else None)
        t = _log_timing(log, "Save summary tables", t)

        # ── Generate analysis report ──
        if args.report:
            report_output = args.report_output or _default_report_path(args.output)
            log.info("Generating %s format analysis report...", args.report.upper())
            generate_report(output_results, report_output, args.report,
                             ref_length=len(ref_seq),
                             ref_seq=ref_seq,
                             cutsites=_get_display_cutsites(ref_seq, config.amplicon_config),
                             allele_window_start=args.allele_window_start,
                             allele_window_end=args.allele_window_end,
                             allele_top_n=args.allele_top_n,
                             version=__version__,
                             summary_data=summary_data,
                             target_region_left=config.amplicon_config.target_region_left if config.amplicon_config else 13,
                             target_region_right=config.amplicon_config.target_region_right if config.amplicon_config else 7)
            t = _log_timing(log, "Generate HTML report", t)

        total = time.perf_counter() - t0
        log.info("  Total time: %.1fs (%.1f min)", total, total / 60)

    else:
        parser.print_help()
        sys.exit(1)


def _default_report_path(output_path: str) -> str:
    """Infer default report path from the --output argument."""
    base = output_path
    for ext in [".json", ".tsv"]:
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    return f"{base}_report"


def _get_display_cutsites(ref_seq: str, amplicon_config=None):
    """Get cutsite information for heatmap annotations."""
    try:
        from crisviper import get_amplicon_structure
        from crisviper.config import AmpliconConfig
        if amplicon_config is None:
            amplicon_config = AmpliconConfig.carlin_standard()
        cs = get_amplicon_structure(ref_seq, config=amplicon_config)
        if cs:
            log.info("Auto-detected %d cutsite regions (for heatmap annotation)", len(cs))
        return cs
    except Exception:
        return None


if __name__ == "__main__":
    main()
