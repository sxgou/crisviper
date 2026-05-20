"""crisviper/models.py — 谱系示踪分析 数据模型

类型安全的数据模型，贯穿整个管道的数据流通。
所有管道阶段之间的数据交换通过这些模型进行。
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Literal
from enum import Enum


# ═══════════════════════════════════════════════════════════════
# 突变类型
# ═══════════════════════════════════════════════════════════════

class MutationType(Enum):
    """突变类型枚举"""
    SUBSTITUTION = "substitution"       # 点突变（替换）
    DELETION = "deletion"               # 删除（query 中有 gap）
    INSERTION = "insertion"             # 插入（ref 中有 gap）
    INDEL = "indel"                     # 复合突变（插入+删除相邻）


# ═══════════════════════════════════════════════════════════════
# 单突变事件
# ═══════════════════════════════════════════════════════════════

@dataclass
class MutationEvent:
    """单个突变事件

    从比对结果中提取的独立突变事件。
    一个比对可能包含多个 MutationEvent。
    """
    type: MutationType                   # 突变类型
    ref_pos: int                         # 参考序列上的起始位置（0-indexed）
    ref_base: str = ""                   # 参考序列上的碱基（substitution/del）
    query_base: str = ""                 # 查询序列上的碱基（substitution/ins）
    length: int = 1                      # 突变长度（bp）
    in_cutsite_window: bool = False      # 是否在 cutsite 窗口中
    raw_ref_segment: str = ""            # 原始比对中 ref 片段
    raw_query_segment: str = ""          # 原始比对中 query 片段
    score: float = 0.0                   # 置信度分数

    def to_dict(self) -> Dict:
        return {
            "type": self.type.value if hasattr(self.type, 'value') else str(self.type),
            "ref_pos": self.ref_pos,
            "ref_base": self.ref_base,
            "query_base": self.query_base,
            "length": self.length,
            "in_cutsite_window": self.in_cutsite_window,
            "score": self.score,
        }


# ═══════════════════════════════════════════════════════════════
# 输入记录
# ═══════════════════════════════════════════════════════════════

@dataclass
class QueryRecord:
    """单条查询序列记录"""
    readName: str
    cellBC: str = "unknown"
    UMI: str = "unknown"
    readCount: int = 1
    seq: str = ""


# ═══════════════════════════════════════════════════════════════
# 比对统计信息
# ═══════════════════════════════════════════════════════════════

@dataclass
class AlignmentStats:
    """比对统计信息"""
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
        """从字典创建（兼容旧代码）"""
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
        """转换为字典（兼容旧代码）"""
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
        """是否有插入或删除"""
        return self.gaps_in_ref > 0 or self.gaps_in_query > 0

    @property
    def has_mutation(self) -> bool:
        """是否有任何突变（点突变或indel）"""
        return self.mismatches > 0 or self.has_indel


# ═══════════════════════════════════════════════════════════════
# 比对结果
# ═══════════════════════════════════════════════════════════════

@dataclass
class AlignmentResult:
    """单条序列的比对结果"""
    query: QueryRecord           # 原始查询序列
    success: bool = True         # 是否成功
    score: float = 0.0           # 比对分数
    aligned_ref: str = ""        # 比对后的参考序列
    aligned_query: str = ""      # 比对后的查询序列
    stats: Optional[AlignmentStats] = None  # 比对统计
    error: str = ""              # 错误消息（失败时）
    mutations: List[MutationEvent] = field(default_factory=list)  # 突变事件列表
    mode: str = "standard"       # 比对模式 (standard / lineage)

    @classmethod
    def error_result(cls, query: QueryRecord, error_msg: str) -> "AlignmentResult":
        """创建一个错误结果"""
        return cls(query=query, success=False, error=error_msg)

    def to_dict(self) -> Dict:
        """转换为字典（兼容旧输出格式）"""
        base = {
            "readName": self.query.readName,
            "cellBC": self.query.cellBC,
            "UMI": self.query.UMI,
            "readCount": self.query.readCount,
        }
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
# 管道配置
# ═══════════════════════════════════════════════════════════════

@dataclass
class PipelineConfig:
    """管道配置 — 所有默认值集中管理

    覆盖此配置即可调整管道的所有行为，
    无需修改任何函数代码。
    """
    # ── 比对参数 ──
    match_score: float = 2.0
    mismatch_penalty: float = -3.0
    gap_open: float = -2.0
    gap_extend: float = -0.1

    # ── 谱系示踪模式 ──
    lineage_mode: bool = False
    # 梯度惩罚: 以 cutsite 为中心平滑变化，替代离散区域倍率
    gradient_mode: bool = True      # 标准模式下也启用位置感知（需要cutsite信息）
    min_scale: float = 1.0          # 切割点处最低惩罚倍率
    max_scale: float = 6.0          # 保守区最高惩罚倍率
    cutsite_edge_scale: float = 2.0  # Cutsite 边界惩罚倍率
    gradient_radius: Optional[float] = None  # 梯度半径: None=自动(多cutsite)/30bp(单)
    mismatch_density_threshold: float = 0.34
    sub_window: int = 3

    # ── Gap退出惩罚 ──
    # gap_exit_strength: Gap exit抑制强度（≤0，0=关闭）
    # 按cutsite位置自动计算梯度惩罚，在cutsite中心产生最强抑制
    # 默认0.0=关闭，推荐-3.5在cutsite中心产生约-21.0峰值惩罚
    gap_exit_strength: float = 0.0

    # ── 短匹配区域折扣 ──
    # short_match_window: 短匹配区域阈值（bp），0=关闭
    # short_match_discount: 短匹配区域match_score折扣系数（0~1），1.0=不打折
    # 例如 window=3, discount=0.5 时，≤3bp的匹配区域match_score减半
    short_match_window: int = 0
    short_match_discount: float = 1.0

    # ── 密集错配区域惩罚 ──
    # dense_mismatch_window: 密集错配检测窗口大小（bp）
    # dense_mismatch_penalty: 密集错配区域额外惩罚（≤0，0=关闭）
    # 负值使DP在密集错配区域倾向选择insertion路径
    dense_mismatch_window: int = 6
    dense_mismatch_penalty: float = 0.0

    # ── 同源区域重复惩罚（跨靶点homology保护） ──
    # homology_window: 同源性检测窗口大小（bp）
    # homology_penalty: ≤0，同源区域match_score减分，0=关闭
    # 负值使DP在参考序列的重复区域不倾向匹配，减少跨靶点错配
    homology_window: int = 8
    homology_penalty: float = 0.0

    # ── 孤立碱基端点和并惩罚（吸收孤立匹配到gap端点） ──
    # 孤立匹配：前后被gap包围的单个碱基匹配
    # 负值使DP倾向将孤立匹配吸收到gap中，减少碎片化比对
    # 与 gap_exit_strength 协同：gap_exit_strength 惩罚gap→M的过渡，
    # isolated_base_penalty 惩罚过渡后只有1bp匹配的场景
    isolated_base_penalty: float = 0.0

    # ── 引物参数 ──
    primer5_len: int = 23
    primer3_len: int = 33
    primer5_threshold: int = 19
    primer3_threshold: int = 29

    # ── Allele过滤（exclusive阈值，>threshold通过） ──
    min_reads_sub: int = 5       # 纯点突变最小readCount（默认>5通过）
    min_reads_indel: int = 0     # 含indel最小readCount（0=不过滤）

    # ── 背景点突变矫正 ──
    correct_bg_sub: bool = True           # 启用背景点突变矫正
    keep_sub_indel_window: int = 3        # 矫正时indel邻近保留窗口(bp)

    # ── 多线程 ──
    threads: int = 1
    chunk_size: int = 500

    # ── 报告 ──
    report_format: Optional[str] = None   # json / html
    allele_top_n: int = 50
    allele_window_start: int = 0
    allele_window_end: Optional[int] = None

    # ── Cutsite配置 ──
    cutsites_path: Optional[str] = None    # JSON配置文件路径
    auto_detect_cutsites: bool = True

    # ── 降噪与等位基因调用 ──
    denoise_enabled: bool = False          # 是否启用UMI/CB降噪
    call_alleles_enabled: bool = False     # 是否启用等位基因调用
    call_alleles_mode: str = "coarse"      # "coarse" 或 "exact"
    dominant_frac: float = 0.5             # 显性等位基因阈值


# ═══════════════════════════════════════════════════════════════
# 管道统计结果
# ═══════════════════════════════════════════════════════════════

@dataclass
class PipelineStats:
    """全管道统计数据（用于报告生成）"""
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
        """编辑效率（百分比）"""
        if self.successful == 0:
            return 0.0
        return self.mutated_sequences / self.successful * 100


@dataclass
class PipelineResult:
    """完整管道运行的最终结果"""
    results: List[AlignmentResult]         # 所有比对结果
    config: PipelineConfig                 # 使用的配置
    stats: PipelineStats                   # 管道统计数据
    ref_length: int = 0                    # 参考序列长度
    mutation_type_counts: Dict = field(default_factory=dict)  # 突变类型计数
    total_mismatches: int = 0              # 点突变总计
    insertion_lengths: List[int] = field(default_factory=list)
    deletion_lengths: List[int] = field(default_factory=list)
    called_alleles: List = field(default_factory=list)  # CalledAllele列表（可选）

    def get_successful(self) -> List[AlignmentResult]:
        """获取所有成功比对的结果"""
        return [r for r in self.results if r.success]

    def get_failed(self) -> List[AlignmentResult]:
        """获取所有失败比对的结果"""
        return [r for r in self.results if not r.success]

    def get_mutated(self) -> List[AlignmentResult]:
        """获取所有包含突变的结果"""
        return [r for r in self.get_successful()
                if r.stats and r.stats.has_mutation]

    def get_unmutated(self) -> List[AlignmentResult]:
        """获取所有无突变的结果"""
        return [r for r in self.get_successful()
                if r.stats and not r.stats.has_mutation]


# ═══════════════════════════════════════════════════════════════
# 突变统计（用于报告）
# ═══════════════════════════════════════════════════════════════

@dataclass
class MutationTypeCounts:
    """各类突变类型的序列数和Reads数"""
    only_insertion: int = 0
    only_deletion: int = 0
    only_substitution: int = 0
    insertion_and_deletion: int = 0
    insertion_and_substitution: int = 0
    deletion_and_substitution: int = 0
    all_three: int = 0

    # 对应reads版本
    only_insertion_reads: int = 0
    only_deletion_reads: int = 0
    only_substitution_reads: int = 0
    insertion_and_deletion_reads: int = 0
    insertion_and_substitution_reads: int = 0
    deletion_and_substitution_reads: int = 0
    all_three_reads: int = 0
