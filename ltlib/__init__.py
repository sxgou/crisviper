"""
ltlib — Lineage Tracer 库
=====================
核心比对算法、谱系示踪管线、后处理矫正和 I/O 函数。

模块划分（按管道流程）：
  - io:           数据转换（FASTQ→TSV/FASTA）
  - alignment:    核心比对算法
  - lineage:      谱系示踪比对模式
  - corrections:  比对后矫正
  - mutation:     突变识别
  - pipeline:     管道编排
  - reporting:    输出和报告生成
  - plotting:     可视化图表
  - config:       数据结构定义
  - models:       类型安全的数据模型
  - logging_config: 日志框架
"""
from ltlib.config import CutsiteRegion, AmpliconConfig
from ltlib.models import (
    QueryRecord, AlignmentResult, AlignmentStats,
    PipelineConfig, PipelineResult, PipelineStats,
    MutationEvent, MutationType,
)
from ltlib.alignment import (
    affine_gap_alignment, calculate_alignment_stats,
    count_gap_blocks, affine_gap_alignment_position_aware,
)
from ltlib.lineage import (
    lineage_tracer_align, build_gap_penalty_profile,
    get_amplicon_structure,
)
from ltlib.corrections import (
    convert_dense_mismatch_to_indel,
    filter_point_mutations,
    correct_repetitive_misalignment,
    correct_target_misalignments,
    remove_isolated_matches,
)
from ltlib.mutation import (
    extract_mutations, classify_mutation_type,
    build_mutation_summary, format_mutations_for_display,
    annotate_mutation, annotate_mutations,
    identify_sequence_events, identify_cas9_events, classify_bp_event,
)
from ltlib.pipeline import Pipeline, align_single, check_primer_anchoring
from ltlib.io import (
    fastq_to_dataframe, fastq_to_fasta, save_tsv,
    read_reference_fasta, read_queries_tsv, read_queries_fasta,
    parse_10x_fastq, parse_10x_provenance,
    filter_sc_cbs_and_umis, parse_indrops_provenance,
)
from ltlib.denoiser import (
    directional_adjacency_top_down_denoiser,
)
from ltlib.threshold import (
    compute_threshold,
)
from ltlib.caller import (
    call_alleles_coarse_grain, call_alleles_exact,
    CalledAllele, _event_structure,
)
from ltlib.reporting import (
    save_alignment_results, generate_report, save_text_report,
)
from ltlib.metrics import (
    effective_alleles, diversity_index, alleles_per_cell,
    singletons_per_cell, carlin_potential,
)
from ltlib.logging_config import get_logger, setup_logging
