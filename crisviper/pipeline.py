"""crisviper/pipeline.py — 谱系示踪分析 管道编排模块

将完整分析流程组织为可配置的 Pipeline，包含：
  1. 数据转换（DataLoader）
  2. 序列比对（Aligner）
  3. 突变识别（MutationDetector）
  4. Allele过滤（AlleleFilter）
  5. 数据统计（StatsCollector）
  6. 输出结果（ResultSaver）

每个步骤可独立测试、替换、跳过。
"""

import os
# 防止 fork + NumPy 线程冲突（在 import numpy 之前必须设置）
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool
from functools import partial
from typing import List, Dict, Optional, Tuple

from crisviper.logging_config import get_logger
from crisviper.models import (
    QueryRecord, AlignmentResult, AlignmentStats,
    PipelineConfig, PipelineResult, PipelineStats,
    MutationEvent, MutationType,
)
from crisviper.config import CutsiteRegion, AmpliconConfig
from crisviper.mutation import extract_mutations, build_mutation_summary, annotate_mutations, _build_ref_pos_map_full
from crisviper.alignment import (
    affine_gap_alignment,
    affine_gap_alignment_position_aware,
    calculate_alignment_stats,
)
from crisviper.lineage import (
    lineage_tracer_align,
    build_gradient_profiles,
    build_homology_penalty_profile,
    get_amplicon_structure,
)
from crisviper.denoiser import directional_adjacency_top_down_denoiser
from crisviper.caller import call_alleles_coarse_grain, call_alleles_exact, CalledAllele

log = get_logger(__name__)


# ═══════════════════════════════════════════════════════════════
# 步骤1: 双端锚定（Primer检测）
# ═══════════════════════════════════════════════════════════════

def _check_primer_quality(
    ar: str,
    aq: str,
    ref_seq: str,
    p5: int,
    p3: int,
    p5_threshold: int,
    p3_threshold: int,
) -> Tuple[bool, int, bool, int]:
    """从全长比对结果中检查引物区比对质量

    在全局比对完成后调用，检查 aligned query 在引物区域的匹配度。
    取代 check_primer_anchoring 的先验检查。

    返回: (p5通过, p5匹配数, p3通过, p3匹配数)
    """
    pos_map, total = _build_ref_pos_map_full(ar)
    if total < p5 + p3:
        return False, 0, False, 0

    # 5' 引物区: ref 位置 0..p5-1
    p5_end = next((i for i, p in enumerate(pos_map) if p == p5 - 1), None)
    if p5_end is None:
        return False, 0, False, 0
    p5_match = sum(1 for i in range(p5_end + 1)
                   if ar[i] != '-' and aq[i] != '-' and ar[i] == aq[i])

    # 3' 引物区: ref 位置 total-p3..total-1
    p3_start = next((i for i, p in enumerate(pos_map) if p == total - p3), None)
    if p3_start is None:
        return False, p5_match, False, 0
    p3_match = sum(1 for i in range(p3_start, len(ar))
                   if ar[i] != '-' and aq[i] != '-' and ar[i] == aq[i])

    return p5_match >= p5_threshold, p5_match, p3_match >= p3_threshold, p3_match


def _extract_internal_region(
    ar: str,
    aq: str,
    r_seq: str,
    p5: int,
    p3: int,
) -> Tuple[Optional[str], Optional[str], str]:
    """从全长比对中提取引物之间的内部区域

    返回: (aligned_internal_ref, aligned_internal_query, ungapped_internal_ref)
    """
    pos_map, total = _build_ref_pos_map_full(ar)

    int_start = next((i for i, p in enumerate(pos_map) if p == p5), None)
    int_end = next((i for i, p in enumerate(pos_map) if p == total - p3 - 1), None)
    if int_start is None or int_end is None:
        return None, None, ""

    ar_int = ar[int_start:int_end + 1]
    aq_int = aq[int_start:int_end + 1]
    int_r = r_seq[p5:total - p3] if p3 > 0 else r_seq[p5:]

    return ar_int, aq_int, int_r


# ═══════════════════════════════════════════════════════════════
# 步骤2: Allele置信度过滤
# ═══════════════════════════════════════════════════════════════

def check_allele_confidence(
    stats: AlignmentStats,
    read_count: int,
    min_reads_sub: int = 5,
    min_reads_indel: int = 0,
) -> Tuple[bool, str]:
    """检查Allele置信度

    规则:
      - 仅点突变（无indel）需 readCount > min_reads_sub
      - 有indel的需 readCount > min_reads_indel
      - 无突变的野生型直接通过

    返回: (通过与否, 原因)
    """
    if not stats.has_indel and stats.mismatches > 0:
        # 仅点突变
        if read_count <= min_reads_sub:
            return False, f"假阳性: 仅点突变但readCount不足({read_count}<={min_reads_sub})"
    elif stats.has_indel:
        # 有indel
        if read_count <= min_reads_indel:
            return False, f"假阳性: 有indel但readCount不足({read_count}<={min_reads_indel})"
    return True, ""


# ═══════════════════════════════════════════════════════════════
# 步骤3: 背景点突变矫正
# ═══════════════════════════════════════════════════════════════

def correct_background_substitutions(
    ar_int: str,
    aq_int: str,
    cutsites: Optional[List[CutsiteRegion]],
    sub_window: int = 3,
    keep_sub_indel_window: int = 3,
    lineage_mode: bool = True,
) -> str:
    """校正背景点突变：将cutsite窗口和indel邻近区外的substitution改回reference

    在合并allele之前调用，消除测序背景点突变导致的等位基因碎片化。

    保留规则（满足任一即保留）：
      1. 谱系模式：cutsite ± sub_window 内的 substitution
      2. 标准模式：cutsite 自身范围内的 substitution
      3. 任意模式：indel 任一端点 ± keep_sub_indel_window 内的 substitution

    Args:
        ar_int: 内部区域比对的参考序列（含gap）
        aq_int: 内部区域比对的查询序列（含gap）
        cutsites: cutsite区域列表（已转换到内部区域坐标）
        sub_window: cutsite邻近保留窗口(bp)
        keep_sub_indel_window: indel邻近保留窗口(bp)
        lineage_mode: 是否谱系示踪模式

    Returns:
        矫正后的 aq_int
    """
    if cutsites is None or not ar_int or len(ar_int) != len(aq_int):
        return aq_int

    pos_map, total = _build_ref_pos_map_full(ar_int)

    # Step 1: 收集 indel 涉及的 ref 坐标区间（扩展 flank）
    indel_keep_regions = []
    i = 0
    alen = len(ar_int)
    while i < alen:
        # Deletion: ref has base, query is gap
        if ar_int[i] != '-' and aq_int[i] == '-':
            start = pos_map[i]
            while i < alen and ar_int[i] != '-' and aq_int[i] == '-':
                i += 1
            end = pos_map[i - 1]
            rs = max(0, start - keep_sub_indel_window)
            re = min(total - 1, end + keep_sub_indel_window)
            indel_keep_regions.append((rs, re))
        # Insertion: ref is gap, query has base
        elif ar_int[i] == '-' and aq_int[i] != '-':
            left_pos = pos_map[i - 1] if i > 0 and pos_map[i - 1] >= 0 else 0
            while i < alen and ar_int[i] == '-' and aq_int[i] != '-':
                i += 1
            rs = max(0, left_pos - keep_sub_indel_window)
            re = min(total - 1, left_pos + keep_sub_indel_window)
            indel_keep_regions.append((rs, re))
        else:
            i += 1

    # Step 2: 构建保留区域集合（区间合并）
    keep_intervals = []  # (start, end) inclusive

    # 2a: cutsite 区域
    for cs in cutsites:
        if lineage_mode:
            ks = max(0, cs.start - sub_window)
            ke = min(total - 1, cs.end + sub_window)
        else:
            ks = max(0, cs.start)
            ke = min(total - 1, cs.end)
        keep_intervals.append((ks, ke))

    # 2b: indel flanking
    keep_intervals.extend(indel_keep_regions)

    if not keep_intervals:
        return aq_int

    # 合并重叠区间
    keep_intervals.sort()
    merged = [keep_intervals[0]]
    for s, e in keep_intervals[1:]:
        if s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))

    # Step 3: 遍历 substitution 列，判断是否保留
    aq_list = list(aq_int)
    for i in range(alen):
        if ar_int[i] != '-' and aq_int[i] != '-' and ar_int[i] != aq_int[i]:
            ref_pos = pos_map[i]
            if ref_pos < 0:
                continue
            in_keep = False
            for ks, ke in merged:
                if ks <= ref_pos <= ke:
                    in_keep = True
                    break
            if not in_keep:
                aq_list[i] = ar_int[i]

    return ''.join(aq_list)


# ═══════════════════════════════════════════════════════════════
# 步骤4: 单序列比对（完整的 align_single 逻辑）
# ═══════════════════════════════════════════════════════════════


def align_single(
    query: QueryRecord,
    ref_seq: str,
    config: PipelineConfig,
    cutsites: Optional[List[CutsiteRegion]] = None,
) -> AlignmentResult:
    """单条序列的完整比对管道

    管道步骤:
      1. 全长全局比对（ref 332bp vs query）
      2. 引物区比对质量检查（取代双端锚定先验检查）
      3. 内部区域提取（trim 掉引物列）
      4. 背景点突变矫正（合并allele前，消除测序背景）
      5. 使用 WT 引物组装全长输出
      6. 内部区域重计分
      7. 突变识别
      8. Allele置信度过滤
    """
    query_seq = query.seq
    q_seq = query_seq
    r_seq = ref_seq
    p5 = config.primer5_len
    p3 = config.primer3_len

    # ── 步骤1: 全长全局比对 ──
    if config.lineage_mode and cutsites:
        score, ar, aq, raw_stats = _align_full_lineage(
            r_seq, q_seq, cutsites, config,
        )
    else:
        score, ar, aq, raw_stats = _align_full_standard(
            r_seq, q_seq, config, cutsites=cutsites,
        )

    if not ar:
        return AlignmentResult.error_result(query, "全序列比对失败")

    # ── 步骤2: 引物区比对质量检查（从比对结果判断）──
    p5_ok, p5_match, p3_ok, p3_match = _check_primer_quality(
        ar, aq, r_seq, p5, p3,
        config.primer5_threshold, config.primer3_threshold,
    )
    if not (p5_ok and p3_ok):
        n5 = config.primer5_threshold
        n3 = config.primer3_threshold
        return AlignmentResult.error_result(
            query,
            f"引物锚定失败: Primer5({p5_match}/{p5}<{n5}) "
            f"Primer3({p3_match}/{p3}<{n3})",
        )

    # ── 步骤3: 提取内部区域 ──
    ar_int, aq_int, int_r = _extract_internal_region(ar, aq, r_seq, p5, p3)
    if ar_int is None:
        return AlignmentResult.error_result(query, "无法从全长比对中提取内部区域")

    # ── 步骤4: 背景点突变矫正（合并allele前） ──
    if config.correct_bg_sub and cutsites is not None:
        int_r_len = len(r_seq) - p5 - p3
        internal_cutsites = _adjust_cutsites_to_internal(cutsites, p5, int_r_len)
        aq_int = correct_background_substitutions(
            ar_int, aq_int,
            cutsites=internal_cutsites,
            sub_window=config.sub_window,
            keep_sub_indel_window=config.keep_sub_indel_window,
            lineage_mode=config.lineage_mode,
        )

    # ── 步骤5: WT 引物组装 ──
    aligned_ref, aligned_query = _assemble_full_length(
        ar_int, aq_int, r_seq, p5, p3,
    )

    # ── 步骤6: 内部区域重计分 ──
    stats_dict = calculate_alignment_stats(ar_int, aq_int)
    score = (stats_dict['matches'] * 2 +
             stats_dict['mismatches'] * (-3) +
             stats_dict['gaps_in_ref'] * (-2) +
             stats_dict['gaps_in_query'] * (-2))
    stats_dict['score'] = score
    stats = AlignmentStats.from_dict(stats_dict)

    # ── 步骤7: 突变识别（仅内部区域）──
    if cutsites is not None:
        int_r_len = len(r_seq) - p5 - p3
        internal_cutsites = _adjust_cutsites_to_internal(cutsites, p5, int_r_len)
    else:
        internal_cutsites = None
    mutations = extract_mutations(
        ar_int, aq_int,
        cutsites=internal_cutsites,
        sub_window=config.sub_window,
    )

    # 转换 ref_pos 从内部坐标 → 全长参考坐标
    for m in mutations:
        m.ref_pos += config.primer5_len

    # ── 步骤8: Allele置信度过滤 ──
    passed, reason = check_allele_confidence(
        stats, query.readCount,
        min_reads_sub=config.min_reads_sub,
        min_reads_indel=config.min_reads_indel,
    )
    if not passed:
        return AlignmentResult.error_result(query, reason)

    return AlignmentResult(
        query=query,
        success=True,
        score=score,
        aligned_ref=aligned_ref,
        aligned_query=aligned_query,
        stats=stats,
        mutations=mutations,
        mode="lineage" if config.lineage_mode else "standard",
    )


# ═══════════════════════════════════════════════════════════════
# 谱系示踪比对模式
# ═══════════════════════════════════════════════════════════════

def _align_full_lineage(
    r_seq: str,
    q_seq: str,
    cutsites: List[CutsiteRegion],
    config: PipelineConfig,
) -> Tuple[float, str, str, Dict]:
    """全长谱系示踪比对 ref vs query（332bp 全局对齐）

    直接调用 lineage_tracer_align 做全长比对，cutsite 已在全长坐标中。
    无需 CGCCG 前缀特殊处理 — DP 在 332bp 上下文中自然处理。
    """
    if not q_seq:
        return 0.0, "", "", {"matches": 0, "mismatches": 0, "gaps_in_ref": 0,
                             "gaps_in_query": 0, "gap_blocks_ref": [],
                             "gap_blocks_query": [], "alignment_length": 0,
                             "similarity": 0.0, "identity": 0.0, "score": 0.0}

    score, ar, aq, raw_stats = lineage_tracer_align(
        r_seq, q_seq, cutsites,
        match_score=config.match_score,
        mismatch_penalty=config.mismatch_penalty,
        base_gap_open=config.gap_open,
        base_gap_extend=config.gap_extend,
        min_scale=config.min_scale,
        max_scale=config.max_scale,
        cutsite_edge_scale=config.cutsite_edge_scale,
        gradient_radius=config.gradient_radius,
        mismatch_density_threshold=config.mismatch_density_threshold,
        sub_window=config.sub_window,
        gap_exit_bonus=0.0,
        base_gap_exit=config.gap_exit_strength,
        short_match_window=config.short_match_window,
        short_match_discount=config.short_match_discount,
        dense_mismatch_window=config.dense_mismatch_window,
        dense_mismatch_penalty=config.dense_mismatch_penalty,
        homology_window=config.homology_window,
        homology_penalty=config.homology_penalty,
        isolated_base_penalty=config.isolated_base_penalty,
    )
    return score, ar, aq, raw_stats


# ═══════════════════════════════════════════════════════════════
# 标准 Gotoh 比对模式
# ═══════════════════════════════════════════════════════════════

def _align_full_standard(
    r_seq: str,
    q_seq: str,
    config: PipelineConfig,
    cutsites: Optional[List[CutsiteRegion]] = None,
) -> Tuple[float, str, str, Dict]:
    """全长标准 Gotoh 比对 ref vs query

    当 cutsites 可用且 config.gradient_mode 启用时，使用平滑梯度
    位置感知的 gap/mismatch 惩罚；否则使用全局固定惩罚。
    """
    if not q_seq:
        return 0.0, "", "", {"matches": 0, "mismatches": 0, "gaps_in_ref": 0,
                             "gaps_in_query": 0, "gap_blocks_ref": [],
                             "gap_blocks_query": [], "alignment_length": 0,
                             "similarity": 0.0, "identity": 0.0, "score": 0.0}

    if cutsites and config.gradient_mode:
        gap_open_p, gap_extend_p, mismatch_p, _ = build_gradient_profiles(
            ref_length=len(r_seq), cutsites=cutsites,
            base_gap_open=config.gap_open, base_gap_extend=config.gap_extend,
            mismatch_penalty=config.mismatch_penalty,
            min_scale=config.min_scale, max_scale=config.max_scale,
            cutsite_edge_scale=config.cutsite_edge_scale,
            gradient_radius=config.gradient_radius,
        )
        homology_p = build_homology_penalty_profile(
            r_seq, homology_window=config.homology_window,
            homology_penalty=config.homology_penalty,
        )
        score, ar, aq, raw_stats = affine_gap_alignment_position_aware(
            r_seq, q_seq, gap_open_p, gap_extend_p,
            match_score=config.match_score,
            mismatch_penalty=config.mismatch_penalty,
            mismatch_penalty_profile=mismatch_p,
            gap_exit_bonus=config.gap_exit_strength,
            short_match_window=config.short_match_window,
            short_match_discount=config.short_match_discount,
            dense_mismatch_window=config.dense_mismatch_window,
            dense_mismatch_threshold=config.mismatch_density_threshold,
            dense_mismatch_penalty=config.dense_mismatch_penalty,
            homology_profile=homology_p,
            isolated_base_penalty=config.isolated_base_penalty,
        )
    else:
        score, ar, aq, raw_stats = affine_gap_alignment(
            r_seq, q_seq,
            match_score=config.match_score,
            mismatch_penalty=config.mismatch_penalty,
            gap_open=config.gap_open,
            gap_extend=config.gap_extend,
            gap_exit_bonus=config.gap_exit_strength,
            short_match_window=config.short_match_window,
            short_match_discount=config.short_match_discount,
        )
    return score, ar, aq, raw_stats


# ═══════════════════════════════════════════════════════════════
# 全长比对组装
# ═══════════════════════════════════════════════════════════════

def _assemble_full_length(
    ar_int: str,
    aq_int: str,
    r_seq: str,
    p5: int,
    p3: int,
) -> Tuple[str, str]:
    """使用 WT 引物序列组装全长比对

    两端直接使用参考序列的引物（WT），避免 query 引物区突变
    导致等位基因碎片化。
    """
    aligned_ref = r_seq[:p5] + ar_int + r_seq[-p3:]
    aligned_query = r_seq[:p5] + aq_int + r_seq[-p3:]
    return aligned_ref, aligned_query


# ═══════════════════════════════════════════════════════════════
# Cutsite坐标调整
# ═══════════════════════════════════════════════════════════════

def _adjust_cutsites_to_internal(
    cutsites: List[CutsiteRegion],
    primer5_len: int,
    int_r_len: int,
) -> List[CutsiteRegion]:
    """将cutsite坐标从全长转换到内部区域"""
    result = []
    for cs in cutsites:
        ns = cs.start - primer5_len
        ne = cs.end - primer5_len
        if ns < 0 and ne < 0:
            continue
        if ne >= int_r_len:
            ne = int_r_len - 1
        if ns < 0:
            ns = 0
        if ns <= ne:
            result.append(CutsiteRegion(name=cs.name, start=ns, end=ne))
    return result


# ═══════════════════════════════════════════════════════════════
# 分块处理（用于并行）
# ═══════════════════════════════════════════════════════════════

def _process_chunk(
    queries_chunk: List[QueryRecord],
    ref_seq: str,
    config: PipelineConfig,
    cutsites: Optional[List[CutsiteRegion]] = None,
) -> List[AlignmentResult]:
    """处理一批序列（供并行调用）"""
    return [
        align_single(q, ref_seq, config, cutsites)
        for q in queries_chunk
    ]


# ═══════════════════════════════════════════════════════════════
# Pipeline 主类
# ═══════════════════════════════════════════════════════════════

class Pipeline:
    """谱系示踪分析管道

    完整工作流:
      1. 全长全局比对         (affine_gap_alignment / lineage_tracer_align)
      2. 引物区比对质量检查     (_check_primer_quality)
      3. 内部区域提取          (_extract_internal_region)

      5. WT 引物组装          (_assemble_full_length)
      6. 内部区域重计分         (calculate_alignment_stats)
      7. 突变识别             (extract_mutations)
      8. Allele置信度过滤      (check_allele_confidence)

    每一步都可以独立测试和替换。
    """

    def __init__(
        self,
        config: Optional[PipelineConfig] = None,
        ref_seq: str = "",
    ):
        self.config = config or PipelineConfig()
        self.ref_seq = ref_seq
        self.cutsites: Optional[List[CutsiteRegion]] = None

    def load_cutsites(self, ref_seq: Optional[str] = None) -> None:
        """加载/检测cutsite位置"""
        seq = ref_seq or self.ref_seq
        if not seq:
            return

        if self.config.cutsites_path:
            # 从JSON文件加载
            import json
            try:
                with open(self.config.cutsites_path) as f:
                    cs_data = json.load(f)
                raw = cs_data.get("cutsites", cs_data)
                self.cutsites = [
                    CutsiteRegion(name=c.get("name", f"Target{i+1}"),
                                  start=c["start"], end=c["end"])
                    for i, c in enumerate(raw)
                ]
                log.info("从配置文件加载 %d 个cutsite区域", len(self.cutsites))
            except Exception as e:
                log.error("读取cutsite配置文件失败 - %s", e)
                sys.exit(1)
        elif self.config.auto_detect_cutsites:
            self.cutsites = get_amplicon_structure(seq)
            if self.cutsites:
                log.info("自动检测到 %d 个cutsite区域", len(self.cutsites))
            else:
                log.warning("无法自动推断cutsite位置")

    def run(
        self,
        queries: List[QueryRecord],
        ref_seq: Optional[str] = None,
    ) -> PipelineResult:
        """运行完整分析管道

        参数:
            queries: 查询序列列表（QueryRecord对象）
            ref_seq: 参考序列（可覆盖构造函数中设置的）

        返回:
            PipelineResult 包含所有比对结果和统计数据
        """
        if ref_seq:
            self.ref_seq = ref_seq
        if not self.ref_seq:
            raise ValueError("参考序列未设置")

        # 自动加载cutsite
        if self.config.lineage_mode and self.cutsites is None:
            self.load_cutsites()

        total = len(queries)
        if total == 0:
            return PipelineResult(
                results=[], config=self.config,
                stats=PipelineStats(),
                ref_length=len(self.ref_seq),
            )

        # ── 并行比对 ──
        mode_label = "谱系示踪比对" if self.config.lineage_mode else "标准比对"

        threads = self.config.threads or 1
        threads = min(threads, total)

        log.info("  使用 %d 个并行进程处理 %d 条序列 (%s)...",
                 threads, total, mode_label)

        # 批次划分
        chunk_size = self.config.chunk_size
        chunks = [queries[i:i + chunk_size] for i in range(0, total, chunk_size)]
        log.info("  共 %d 个批次 (每批最多 %d 条)...", len(chunks), chunk_size)

        # 并行执行
        results = []

        if threads > 1:
            try:
                with ProcessPoolExecutor(max_workers=threads) as executor:
                    chunk_func = partial(
                        _process_chunk,
                        ref_seq=self.ref_seq,
                        config=self.config,
                        cutsites=self.cutsites,
                    )
                    futures = {executor.submit(chunk_func, ch): ch for ch in chunks}
                    for future in as_completed(futures):
                        try:
                            results.extend(future.result())
                        except Exception as e:
                            log.error("  批次处理失败 (跳过): %s", e)
            except BrokenProcessPool:
                log.warning("ProcessPoolExecutor 崩溃 (BrokenProcessPool)，回退到单线程处理")
        if not results:
            # ProcessPoolExecutor 失败或 threads <= 1
            log.info("  使用单线程处理 %d 条序列...", total)
            for chunk in chunks:
                results.extend(_process_chunk(chunk, self.ref_seq, self.config, self.cutsites))

        # ── 等位基因调用（可选） ──
        called_alleles = []
        successful_results = [r for r in results if r.success]
        if self.config.call_alleles_enabled and successful_results:
            method = call_alleles_coarse_grain if self.config.call_alleles_mode == "coarse" else call_alleles_exact
            called_alleles = method(successful_results, dominant_frac=self.config.dominant_frac)
            log.info("等位基因调用 (%s): %d 个等位基因",
                     self.config.call_alleles_mode, len(called_alleles))

        # ── 构建统计数据 ──
        return self._build_pipeline_result(results, called_alleles)

    def _build_pipeline_result(self, results: List[AlignmentResult],
                               called_alleles: list = None) -> PipelineResult:
        """从比对结果构建管道统计"""
        stats = PipelineStats()
        stats.total_queries = len(results)

        # 分类丢弃原因
        n_anchor = sum(1 for r in results if not r.success and "锚定" in r.error)
        n_noise = sum(1 for r in results if not r.success and "假阳性" in r.error)

        # 只保留成功结果
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        stats.successful = len(successful)
        stats.failed = len(failed)
        stats.total_reads = sum(r.query.readCount for r in successful)
        stats.n_anchor_failed = n_anchor
        stats.n_noise_filtered = n_noise

        # 突变计数
        mutated = [r for r in successful if r.stats and r.stats.has_mutation]
        unmutated = [r for r in successful if r.stats and not r.stats.has_mutation]
        stats.mutated_sequences = len(mutated)
        stats.unmutated_sequences = len(unmutated)
        stats.mutated_reads = sum(r.query.readCount for r in mutated)

        # 突变汇总
        summary = build_mutation_summary(successful)

        # 日志
        if stats.failed > 0:
            parts = []
            if n_anchor:
                parts.append(f"引物锚定失败 {n_anchor}")
            if n_noise:
                parts.append(f"假阳性allele {n_noise}")
            other = stats.failed - n_anchor - n_noise
            if other:
                parts.append(f"其他 {other}")
            log.warning("已丢弃 %d 条：%s", stats.failed, "，".join(parts))
        log.info("批量比对完成: %d 条有效结果", stats.successful)

        return PipelineResult(
            results=results,
            config=self.config,
            stats=stats,
            ref_length=len(self.ref_seq),
            mutation_type_counts=summary["type_counts"],
            total_mismatches=summary["total_mismatches"],
            insertion_lengths=summary["ins_lengths"],
            deletion_lengths=summary["del_lengths"],
            called_alleles=called_alleles or [],
        )
