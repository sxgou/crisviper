"""
crisviper — CrisViper 库
=====================
核心比对算法、谱系示踪管线、后处理矫正和 I/O 函数。

模块划分（按管道流程）：
  - io:           数据转换（FASTQ→TSV/FASTA）
  - alignment:    核心比对算法
  - lineage:      谱系示踪比对模式
  - mutation:     突变识别
  - pipeline:     管道编排
  - reporting:    输出和报告生成
  - plotting:     可视化图表
  - config:       数据结构定义
  - models:       类型安全的数据模型
  - logging_config: 日志框架
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
    build_mutation_summary, format_mutations_for_display,
    annotate_mutation, annotate_mutations,
    classify_bp_event,
)
from crisviper.pipeline import Pipeline, align_single, check_primer_anchoring
from crisviper.io import (
    fastq_to_dataframe, fastq_to_fasta, save_tsv,
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
