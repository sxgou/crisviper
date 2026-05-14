"""ltlib/pipeline.py — 谱系示踪分析 管道编排模块

将完整分析流程组织为可配置的 Pipeline，包含：
  1. 数据转换（DataLoader）
  2. 序列比对（Aligner）
  3. 突变识别（MutationDetector）
  4. 序列矫正（CorrectionPipeline）
  5. Allele过滤（AlleleFilter）
  6. 数据统计（StatsCollector）
  7. 输出结果（ResultSaver）

每个步骤可独立测试、替换、跳过。
"""

import sys
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
from typing import List, Dict, Optional, Tuple, Callable

from ltlib.logging_config import get_logger
from ltlib.models import (
    QueryRecord, AlignmentResult, AlignmentStats,
    PipelineConfig, PipelineResult, PipelineStats,
    MutationEvent, MutationType,
)
from ltlib.config import CutsiteRegion, AmpliconConfig
from ltlib.mutation import extract_mutations, build_mutation_summary, annotate_mutations, _build_ref_pos_map_full
from ltlib.corrections import (
    convert_dense_mismatch_to_indel,
    filter_point_mutations,
    correct_repetitive_misalignment,
    correct_target_misalignments,
    remove_isolated_matches,
)
from ltlib.alignment import (
    affine_gap_alignment,
    affine_gap_alignment_position_aware,
    calculate_alignment_stats,
)
from ltlib.lineage import (
    lineage_tracer_align,
    build_gap_penalty_profile,
    get_amplicon_structure,
)
from ltlib.denoiser import directional_adjacency_top_down_denoiser
from ltlib.caller import call_alleles_coarse_grain, call_alleles_exact, CalledAllele

log = get_logger(__name__)


def _find_anchor_start(ar_bases: str, int_r: str, n_anchor: int = 20) -> int:
    """在 ref 中查找比对 ref 碱基的最佳起始位置"""
    if not ar_bases:
        return 0
    n_anchor = min(n_anchor, len(ar_bases))
    anchor = ar_bases[:n_anchor]
    best_pos, best_cnt = 0, -1
    for i in range(len(int_r) - n_anchor + 1):
        cnt = sum(1 for k in range(n_anchor) if int_r[i + k] == anchor[k])
        if cnt > best_cnt:
            best_cnt = cnt
            best_pos = i
    return best_pos


# ═══════════════════════════════════════════════════════════════
# 步骤1: 双端锚定（Primer检测）
# ═══════════════════════════════════════════════════════════════

def check_primer_anchoring(
    query_seq: str,
    ref_seq: str,
    primer5_len: int = 23,
    primer3_len: int = 33,
    primer5_threshold: int = 19,
    primer3_threshold: int = 29,
) -> Tuple[bool, int, int]:
    """检查 query 是否能双端锚定到 ref (deprecated — 推荐使用 _check_primer_quality)

    返回: (通过与否, p5匹配数, p3匹配数)
    """
    if len(query_seq) < primer5_len + primer3_len:
        return False, 0, 0

    p5_matches = sum(1 for a, b in zip(query_seq[:primer5_len], ref_seq[:primer5_len]) if a == b)
    q3_check = query_seq[-primer3_len:] if len(query_seq) >= primer3_len else ""
    p3_matches = sum(1 for a, b in zip(q3_check, ref_seq[-primer3_len:]) if a == b)

    return (p5_matches >= primer5_threshold and p3_matches >= primer3_threshold,
            p5_matches, p3_matches)


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


def apply_corrections_pipeline(
    aligned_ref: str,
    aligned_query: str,
    ref_seq: str,
    query_seq: str,
    correction_list: List[str],
    cutsites: Optional[List[CutsiteRegion]] = None,
    mismatch_density_threshold: float = 0.34,
    mutation_window: int = 3,
    config: Optional[AmpliconConfig] = None,
) -> Tuple[str, str]:
    """应用矫正函数管线

    correction_list 中的函数名按顺序执行。
    每个函数对齐后的输出作为下一个的输入。

    参数:
        config: 可选的 AmpliconConfig。提供时启用动态坐标检测。

    返回: (矫正后的 aligned_ref, aligned_query)
    """
    ar, aq = aligned_ref, aligned_query

    for fn_name in correction_list:
        try:
            if fn_name == "convert_dense_mismatch_to_indel":
                ar, aq, _ = convert_dense_mismatch_to_indel(
                    ar, aq,
                    threshold=mismatch_density_threshold,
                )
            elif fn_name == "correct_repetitive_misalignment":
                ar, aq, _ = correct_repetitive_misalignment(ar, aq, ref_seq)
            elif fn_name == "correct_target_misalignments":
                if config is not None:
                    cs_start = len(config.prefix) + 1 * config.period + config.cutsite_offset - 1
                    cs_end = len(config.prefix) + (config.n_targets - 1) * config.period + config.cutsite_offset - 1
                else:
                    cs_start = cs_end = None
                ar, aq, _ = correct_target_misalignments(ar, aq, ref_seq, cs_start, cs_end)
            elif fn_name == "remove_isolated_matches":
                ar, aq, _ = remove_isolated_matches(ar, aq)
            elif fn_name == "filter_point_mutations":
                if cutsites:
                    ar, aq, _ = filter_point_mutations(
                        ar, aq, cutsites, window=mutation_window,
                    )
            else:
                log.warning("未知的矫正函数: %s，跳过", fn_name)
        except Exception as e:
            log.warning("矫正函数 %s 执行失败: %s，继续尝试后续矫正", fn_name, e)

    return ar, aq


# ═══════════════════════════════════════════════════════════════
# 步骤2: Allele置信度过滤
# ═══════════════════════════════════════════════════════════════

def check_allele_confidence(
    stats: AlignmentStats,
    read_count: int,
    min_reads_snv: int = 10,
    min_reads_indel: int = 3,
) -> Tuple[bool, str]:
    """检查Allele置信度

    规则:
      - 仅点突变（无indel）需 readCount >= min_reads_snv
      - 有indel的需 readCount >= min_reads_indel
      - 无突变的野生型直接通过

    返回: (通过与否, 原因)
    """
    if not stats.has_indel and stats.mismatches > 0:
        # 仅点突变
        if read_count < min_reads_snv:
            return False, f"假阳性: 仅点突变但readCount不足({read_count}<{min_reads_snv})"
    elif stats.has_indel:
        # 有indel
        if read_count < min_reads_indel:
            return False, f"假阳性: 有indel但readCount不足({read_count}<{min_reads_indel})"
    return True, ""


# ═══════════════════════════════════════════════════════════════
# 步骤3: 单序列比对（完整的 align_single 逻辑）
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
      4. 矫正管线
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
            r_seq, q_seq, config,
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
            f"引物比对质量不足: Primer5({p5_match}/{p5}<{n5}) "
            f"Primer3({p3_match}/{p3}<{n3})",
        )

    # ── 步骤3: 提取内部区域 ──
    ar_int, aq_int, int_r = _extract_internal_region(ar, aq, r_seq, p5, p3)
    if ar_int is None:
        return AlignmentResult.error_result(query, "无法从全长比对中提取内部区域")

    # ── 步骤4: 矫正管线 ──
    corrections_to_apply = config.corrections[:]
    if config.lineage_mode and cutsites:
        corrections_to_apply = [c for c in corrections_to_apply
                                if c not in ("convert_dense_mismatch_to_indel", "filter_point_mutations")]
    ar_int_corrected, aq_int_corrected = apply_corrections_pipeline(
        ar_int, aq_int, int_r, int_r,
        correction_list=corrections_to_apply,
        cutsites=cutsites,
        mismatch_density_threshold=config.mismatch_density_threshold,
        mutation_window=config.mutation_window,
    )
    if ar_int_corrected != ar_int or aq_int_corrected != aq_int:
        ar_int, aq_int = ar_int_corrected, aq_int_corrected

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
        mutation_window=config.mutation_window,
    )

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
        cutsite_gap_scale=config.cutsite_gap_scale,
        flank_gap_scale=config.flank_gap_scale,
        far_gap_scale=config.far_gap_scale,
        flank_width=config.flank_width,
        mismatch_density_threshold=config.mismatch_density_threshold,
        mutation_window=config.mutation_window,
    )
    return score, ar, aq, raw_stats


# ═══════════════════════════════════════════════════════════════
# 标准 Gotoh 比对模式
# ═══════════════════════════════════════════════════════════════

def _align_full_standard(
    r_seq: str,
    q_seq: str,
    config: PipelineConfig,
) -> Tuple[float, str, str, Dict]:
    """全长标准 Gotoh 比对 ref vs query"""
    if not q_seq:
        return 0.0, "", "", {"matches": 0, "mismatches": 0, "gaps_in_ref": 0,
                             "gaps_in_query": 0, "gap_blocks_ref": [],
                             "gap_blocks_query": [], "alignment_length": 0,
                             "similarity": 0.0, "identity": 0.0, "score": 0.0}

    score, ar, aq, raw_stats = affine_gap_alignment(
        r_seq, q_seq,
        match_score=config.match_score,
        mismatch_penalty=config.mismatch_penalty,
        gap_open=config.gap_open,
        gap_extend=config.gap_extend,
    )
    return score, ar, aq, raw_stats


# ═══════════════════════════════════════════════════════════════
# 片段化比对修复
# ═══════════════════════════════════════════════════════════════

def _fix_fragmented_alignment(
    ar_int: str,
    aq_int: str,
    int_r: str,
    int_q: str,
    config: PipelineConfig,
) -> Tuple[str, str, float]:
    """修复因重复元件导致的片段化或位置偏移比对"""
    # ── 从 aq_int 中提取所有非 gap 片段 ──
    blocks = []
    in_block = False
    for i, c in enumerate(aq_int):
        if c != '-':
            if not in_block:
                block_start = i
                in_block = True
        else:
            if in_block:
                blocks.append((block_start, i, aq_int[block_start:i]))
                in_block = False
    if in_block:
        blocks.append((block_start, len(aq_int), aq_int[block_start:]))

    if not blocks:
        return ar_int, aq_int, 0.0

    # ── 判断是否需要矫正 ──
    first_nongap = blocks[0][0]
    last_block_end = blocks[-1][1]

    shifted = first_nongap > 20
    scattered = len(blocks) >= 2 and last_block_end < len(int_r) - 5
    needs_anchoring = shifted or scattered or len(blocks) > 2

    if not (needs_anchoring and len(int_q) <= len(int_r) - 10):
        if len(blocks) > 2:
            # 轻度片段化：基于参考位置重新放置每个片段
            return _place_blocks_at_ref_positions(ar_int, blocks)
        return ar_int, aq_int, 0.0

    # ── 计算每个片段在 int_r 上的起止坐标（strip gap） ──
    block_ref_positions = []
    for start, end, seq in blocks:
        r_start = sum(1 for c in ar_int[:start] if c != '-')
        r_end = sum(1 for c in ar_int[:end] if c != '-')
        block_ref_positions.append((r_start, r_end, start, seq))

    # ── 安全守卫：检测片段间是否存在真实缺失 ──
    #   间距 > 3 且间隔区域非简单重复 → 视为真实 indel，跳过修复
    has_real_indel = False
    for i in range(len(block_ref_positions) - 1):
        _, cur_end, _, _ = block_ref_positions[i]
        next_start, _, _, _ = block_ref_positions[i + 1]
        gap = next_start - cur_end
        if gap > 3:
            gap_seq = int_r[cur_end:next_start]
            if len(gap_seq) > 3 and len(set(gap_seq)) > 2:
                has_real_indel = True
                break

    if has_real_indel:
        return ar_int, aq_int, 0.0

    # ── 基于参考位置构建新的 aq_int ──
    #   每个片段放在 ar_int 中对应参考位置处，而非暴力拼接
    new_aq = ['-'] * len(ar_int)
    for r_start, _r_end, _orig_start, seq in block_ref_positions:
        insert_pos = 0
        non_gap_cnt = 0
        while insert_pos < len(ar_int) and non_gap_cnt < r_start:
            if ar_int[insert_pos] != '-':
                non_gap_cnt += 1
            insert_pos += 1
        if non_gap_cnt == r_start:
            limit = min(len(seq), len(new_aq) - insert_pos)
            for i in range(limit):
                new_aq[insert_pos + i] = seq[i]

    new_aq_str = ''.join(new_aq)

    # ── 基于真实比对位置计算分数 ──
    matches = sum(1 for i in range(len(new_aq_str))
                  if new_aq_str[i] != '-' and new_aq_str[i] == ar_int[i])
    mismatches = sum(1 for i in range(len(new_aq_str))
                     if new_aq_str[i] != '-' and new_aq_str[i] != ar_int[i])
    score = (matches * config.match_score +
             mismatches * config.mismatch_penalty)
    return ar_int, new_aq_str, score


def _place_blocks_at_ref_positions(
    ar_int: str,
    blocks: list,
) -> Tuple[str, str, float]:
    """将 aq_int 的各个非 gap 片段按其参考坐标重新放置。

    每个片段原本在 ar_int 中有对应的参考区域，此函数保持每个片段
    与其参考区域的对应关系，而非简单拼接。
    """
    new_aq = ['-'] * len(ar_int)
    for start, end, seq in blocks:
        # 计算此片段覆盖的参考区域起止（去除 ar_int 中的 gap）
        r_start = sum(1 for c in ar_int[:start] if c != '-')
        # 找到 ar_int 中第 r_start 个非 gap 字符的起始位置
        insert_pos = 0
        cnt = 0
        while insert_pos < len(ar_int) and cnt < r_start:
            if ar_int[insert_pos] != '-':
                cnt += 1
            insert_pos += 1
        if cnt == r_start:
            limit = min(len(seq), len(new_aq) - insert_pos)
            for i in range(limit):
                new_aq[insert_pos + i] = seq[i]
    return ar_int, ''.join(new_aq), 0.0


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
      4. 序列矫正             (apply_corrections_pipeline)
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
        log.info("  使用 %d 个并行进程处理 %d 条序列 (%s)...",
                 self.config.threads or mp.cpu_count(), total, mode_label)

        # 批次划分
        chunk_size = self.config.chunk_size
        chunks = [queries[i:i + chunk_size] for i in range(0, total, chunk_size)]
        log.info("  共 %d 个批次 (每批最多 %d 条)...", len(chunks), chunk_size)

        # 并行执行
        threads = self.config.threads or mp.cpu_count()
        threads = min(threads, total)
        results = []

        if threads > 1:
            with ProcessPoolExecutor(max_workers=threads) as executor:
                chunk_func = partial(
                    _process_chunk,
                    ref_seq=self.ref_seq,
                    config=self.config,
                    cutsites=self.cutsites,
                )
                futures = {executor.submit(chunk_func, ch): ch for ch in chunks}
                for future in as_completed(futures):
                    results.extend(future.result())
        else:
            # 单线程
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
