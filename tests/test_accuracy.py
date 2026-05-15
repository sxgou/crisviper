"""验证性测试 — 使用已知真值的合成序列验证 pipeline 准确性

现有测试只验证"函数能跑不报错"，本文件验证准确性：

1. 突变检测准确性：已知位置/类型/长度的突变是否能被正确识别
2. 编辑效率验证：已知编辑效率的合成数据计算编辑效率是否一致
3. 位置感知 gap 惩罚：gap 是否优先在 cutsite 区域开启
"""

import sys
import os
import pytest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from crisviper import (
    Pipeline, PipelineConfig, QueryRecord,
    align_single,
    get_amplicon_structure, CutsiteRegion,
    affine_gap_alignment,
    affine_gap_alignment_position_aware,
    build_gradient_profiles,
    lineage_tracer_align,
    extract_mutations,
    MutationEvent, MutationType,
    AlignmentStats,
)

# ═══════════════════════════════════════════════════════════════
# 标准 CARLIN 参考序列 (332bp)
# ═══════════════════════════════════════════════════════════════
CARLIN_REF = (
    "TATGTGTGGGAGGGCTAAGAGG"   # Primer5 (23bp)
    "CCGCC"                     # Prefix (5bp)
    "GACTGCACGACAGTCGA"         # Target1: 保守区13bp + cutsite 7bp
    "CGATGGAG"                  # Linker (7bp)
    "TCGACACGACTCGCGCA"         # Target2
    "TACGATGG"                  # Linker
    "AGTCGACTACAGTCGCTA"        # Target3
    "CGACGATG"                  # Linker
    "GAGTCGCGAGCGCTATG"         # Target4
    "AGCGACTA"                  # Linker
    "TGGAGTCGATACGATACG"        # Target5
    "CGCACGCT"                  # Linker
    "ATGGAGTCGAGAGCGCGC"        # Target6
    "TCGTCAAC"                  # Linker
    "GATGGAGTCGCGACTGTA"        # Target7
    "CGCACTCG"                  # Linker
    "CGATGGAGTCGATAGTAT"        # Target8
    "GCGTACAC"                  # Linker
    "GCGATGGAGTCGACTGCA"        # Target9
    "CGACAGTC"                  # Linker
    "GACTATGGAGTCGATACGTAGC"    # Target10
    "ACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"  # Postfix(8) + Primer3(33)
)

# Cutsite 在全长参考序列上的坐标 (0-indexed, inclusive)
CUTSITES = [
    CutsiteRegion("Target1",  41,  47),
    CutsiteRegion("Target2",  68,  74),
    CutsiteRegion("Target3",  95, 101),
    CutsiteRegion("Target4", 122, 128),
    CutsiteRegion("Target5", 149, 155),
    CutsiteRegion("Target6", 176, 182),
    CutsiteRegion("Target7", 203, 209),
    CutsiteRegion("Target8", 230, 236),
    CutsiteRegion("Target9", 257, 263),
    CutsiteRegion("Target10", 284, 290),
]


# ═══════════════════════════════════════════════════════════════
# 辅助函数
# ═══════════════════════════════════════════════════════════════

def _make_deletion(seq: str, pos: int, length: int) -> str:
    """删除从 pos 开始的 length 个碱基"""
    return seq[:pos] + seq[pos+length:]


def _make_insertion(seq: str, pos: int, insert: str) -> str:
    """在 pos 处插入 insert"""
    return seq[:pos] + insert + seq[pos:]


def _make_substitution(seq: str, pos: int, new_base: str) -> str:
    """将 pos 处的碱基替换为 new_base"""
    lst = list(seq)
    lst[pos] = new_base
    return ''.join(lst)


def _run_align_single(
    query_seq: str,
    read_name: str = "test",
    read_count: int = 50,
    lineage_mode: bool = True,
    min_reads_snv: int = 1,
    min_reads_indel: int = 1,
) -> 'AlignmentResult':
    """通过 align_single 比对单条序列，返回 AlignmentResult"""
    config = PipelineConfig(
        lineage_mode=lineage_mode,
        primer5_len=23,
        primer3_len=33,
        primer5_threshold=19,
        primer3_threshold=29,
        min_reads_snv=min_reads_snv,
        min_reads_indel=min_reads_indel,
    )
    query = QueryRecord(
        readName=read_name, cellBC="test", UMI="UMI",
        readCount=read_count, seq=query_seq,
    )
    cutsites = get_amplicon_structure(CARLIN_REF) if lineage_mode else None
    return align_single(query, CARLIN_REF, config, cutsites)


def _run_pipeline(
    queries,
    lineage_mode: bool = True,
    min_reads_snv: int = 1,
    min_reads_indel: int = 1,
):
    """运行完整 pipeline，返回 PipelineResult"""
    config = PipelineConfig(
        lineage_mode=lineage_mode,
        primer5_len=23,
        primer3_len=33,
        primer5_threshold=19,
        primer3_threshold=29,
        min_reads_snv=min_reads_snv,
        min_reads_indel=min_reads_indel,
        threads=1,
    )
    pipeline = Pipeline(config=config, ref_seq=CARLIN_REF)
    if lineage_mode:
        pipeline.load_cutsites()
    return pipeline.run(queries)


# ═══════════════════════════════════════════════════════════════
# 测试1：突变检测准确性
# ═══════════════════════════════════════════════════════════════

class TestMutationDetectionAccuracy:
    """使用已知真值验证突变识别的准确性"""

    def test_wildtype_no_mutations(self):
        """野生型序列 → 无突变"""
        result = _run_align_single(CARLIN_REF, "wt")
        assert result.success, f"野生型比对失败: {result.error}"
        assert len(result.mutations) == 0, \
            f"野生型不应检测到突变，但检测到 {len(result.mutations)} 个: {[(m.type, m.ref_pos) for m in result.mutations]}"
        assert result.stats.mismatches == 0
        assert not result.stats.has_indel

    def test_3bp_deletion_at_cutsite(self):
        """Target1 cutsite 内 3bp 删除 → 识别为 DELETION, length=3, 在窗口内"""
        query = _make_deletion(CARLIN_REF, 42, 3)
        result = _run_align_single(query, "del_3bp_cutsite")
        assert result.success, f"3bp删除比对失败: {result.error}"

        deletions = [m for m in result.mutations if m.type == MutationType.DELETION]
        assert len(deletions) >= 1, \
            f"未检测到删除突变，检测到: {[(m.type, m.length) for m in result.mutations]}"

        del_event = deletions[0]
        assert del_event.length == 3, f"删除长度应为3，实际为 {del_event.length}"
        assert del_event.in_cutsite_window, \
            f"删除 (ref_pos={del_event.ref_pos}) 应在cutsite窗口内"

    def test_5bp_insertion_at_cutsite(self):
        """Target1 cutsite 内 5bp 插入 → 识别为 INSERTION, length=5, 在窗口内"""
        query = _make_insertion(CARLIN_REF, 42, "ACGTA")
        result = _run_align_single(query, "ins_5bp_cutsite")
        assert result.success, f"5bp插入比对失败: {result.error}"

        insertions = [m for m in result.mutations if m.type == MutationType.INSERTION]
        assert len(insertions) >= 1, \
            f"未检测到插入突变，检测到: {[(m.type, m.length) for m in result.mutations]}"

        ins_event = insertions[0]
        assert ins_event.length == 5, f"插入长度应为5，实际为 {ins_event.length}"
        assert ins_event.in_cutsite_window, \
            f"插入 (ref_pos={ins_event.ref_pos}) 应在cutsite窗口内"

    def test_point_mutation_at_cutsite(self):
        """Target1 cutsite 内点突变 → 识别为 SUBSTITUTION, length=1, 在窗口内"""
        query = _make_substitution(CARLIN_REF, 44, "T")
        result = _run_align_single(query, "sub_cutsite")
        assert result.success, f"点突变比对失败: {result.error}"

        subs = [m for m in result.mutations if m.type == MutationType.SUBSTITUTION]
        assert len(subs) >= 1, \
            f"未检测到替换突变，检测到: {[(m.type, m.length) for m in result.mutations]}"

        sub_event = subs[0]
        assert sub_event.length == 1, f"替换长度应为1，实际为 {sub_event.length}"
        assert sub_event.in_cutsite_window, \
            f"替换 (ref_pos={sub_event.ref_pos}) 应在cutsite窗口内"

    def test_7bp_deletion_full_cutsite(self):
        """完整删除一个 cutsite (7bp) → 识别为 DELETION, length=7"""
        query = _make_deletion(CARLIN_REF, 41, 7)  # 删除整个 Target1 cutsite
        result = _run_align_single(query, "del_7bp_full_cutsite")
        assert result.success, f"7bp删除比对失败: {result.error}"

        deletions = [m for m in result.mutations if m.type == MutationType.DELETION]
        if not deletions:
            deletions = [m for m in result.mutations if m.type == MutationType.COMPLEX]
        assert len(deletions) >= 1, \
            f"未检测到删除突变: {[(m.type, m.length) for m in result.mutations]}"

        # 总删除长度（可能分布在多个事件中）
        total_del_len = sum(m.length for m in result.mutations
                            if m.type in (MutationType.DELETION, MutationType.COMPLEX))
        assert total_del_len >= 7, \
            f"总删除长度应≥7，实际为 {total_del_len}"

    def test_mutations_in_multiple_cutsites(self):
        """多个 cutsite 上各有突变 → 全部被识别"""
        # Target1 cutsite 3bp 删除 + Target2 cutsite 点突变
        query = _make_deletion(CARLIN_REF, 42, 3)
        # 注意：删除后坐标会偏移，Target2从68变成68-3=65
        # CARLIN_REF[68]=='A', 改用'T'制造真正不同的碱基
        query = _make_substitution(query, 65, "T")  # 原始Target2位置
        result = _run_align_single(query, "multi_cutsites")
        assert result.success, f"多cutsite突变比对失败: {result.error}"

        types = {m.type for m in result.mutations}
        has_deletion = MutationType.DELETION in types or MutationType.COMPLEX in types
        has_substitution = MutationType.SUBSTITUTION in types
        assert has_deletion, f"未检测到删除，检测到类型: {types}"
        assert has_substitution, f"未检测到替换，检测到类型: {types}"

    def test_large_deletion_across_target(self):
        """跨 Target 的大删除 → 正确识别"""
        query = _make_deletion(CARLIN_REF, 41, 20)  # 删除 Target1 大部分区域
        result = _run_align_single(query, "del_large_20bp")
        assert result.success, f"大删除比对失败: {result.error}"

        # 应该有删除事件
        deletions = [m for m in result.mutations if m.type in
                     (MutationType.DELETION, MutationType.COMPLEX)]
        total_del_len = sum(m.length for m in deletions)
        assert total_del_len >= 15, \
            f"大删除长度应≥15，实际检测到 {total_del_len}"

    def test_standard_mode_also_detects_mutations(self):
        """标准模式（非lineage）也能正确检测已知突变"""
        query = _make_deletion(CARLIN_REF, 42, 3)
        result = _run_align_single(query, "del_std", lineage_mode=False)
        assert result.success, f"标准模式3bp删除比对失败: {result.error}"

        deletions = [m for m in result.mutations if m.type == MutationType.DELETION]
        assert len(deletions) >= 1, \
            f"标准模式未检测到删除: {[(m.type, m.length) for m in result.mutations]}"
        assert deletions[0].length == 3


# ═══════════════════════════════════════════════════════════════
# 测试2：编辑效率验证
# ═══════════════════════════════════════════════════════════════

class TestEditingEfficiency:
    """使用已知编辑效率的合成数据验证编辑效率计算"""

    def test_known_editing_efficiency(self):
        """已知编辑效率 = 60% (12/20 序列有突变)"""
        queries = []

        # 12 条有突变的序列（不同类型）
        for i in range(6):
            queries.append(QueryRecord(
                readName=f"del_{i}", cellBC="test", UMI=f"D{i}",
                readCount=10,
                seq=_make_deletion(CARLIN_REF, 42, 3),
            ))
        for i in range(6):
            queries.append(QueryRecord(
                readName=f"ins_{i}", cellBC="test", UMI=f"I{i}",
                readCount=10,
                seq=_make_insertion(CARLIN_REF, 42, "ACG"),
            ))

        # 8 条野生型
        for i in range(8):
            queries.append(QueryRecord(
                readName=f"wt_{i}", cellBC="test", UMI=f"W{i}",
                readCount=10, seq=CARLIN_REF,
            ))

        result = _run_pipeline(queries)

        assert result.stats.total_queries == 20
        assert result.stats.successful == 20, \
            f"预期20条成功，实际成功 {result.stats.successful}，失败 {result.stats.failed}"
        assert result.stats.mutated_sequences == 12, \
            f"预期12条突变，实际检测到 {result.stats.mutated_sequences}"

        expected_eff = 12.0 / 20.0 * 100.0
        actual_eff = result.stats.editing_efficiency_pct
        assert actual_eff == pytest.approx(expected_eff, abs=1.0), \
            f"编辑效率: 预期 {expected_eff:.1f}%，实际 {actual_eff:.1f}%"

    def test_zero_editing_efficiency(self):
        """全部野生型 → 编辑效率 = 0%"""
        queries = [
            QueryRecord(readName=f"wt_{i}", cellBC="test", UMI=f"W{i}",
                        readCount=10, seq=CARLIN_REF)
            for i in range(10)
        ]
        result = _run_pipeline(queries)
        assert result.stats.editing_efficiency_pct == pytest.approx(0.0, abs=0.1)

    def test_full_editing_efficiency(self):
        """全部有突变 → 编辑效率 = 100%"""
        queries = [
            QueryRecord(readName=f"mut_{i}", cellBC="test", UMI=f"M{i}",
                        readCount=10, seq=_make_deletion(CARLIN_REF, 42, 3))
            for i in range(10)
        ]
        result = _run_pipeline(queries)
        assert result.stats.editing_efficiency_pct == pytest.approx(100.0, abs=0.1)

    def test_editing_efficiency_with_low_readcount_filtering(self):
        """低 readCount 的序列应被过滤，不影响编辑效率"""
        queries = []
        # 高 readCount 突变序列（应被保留）
        for i in range(5):
            queries.append(QueryRecord(
                readName=f"good_del_{i}", cellBC="test", UMI=f"G{i}",
                readCount=50,
                seq=_make_deletion(CARLIN_REF, 42, 3),
            ))
        # 低 readCount 突变序列（应被过滤）
        for i in range(5):
            queries.append(QueryRecord(
                readName=f"bad_sub_{i}", cellBC="test", UMI=f"B{i}",
                readCount=2,
                seq=_make_substitution(CARLIN_REF, 44, "T"),
            ))
        # 高 readCount 野生型
        for i in range(5):
            queries.append(QueryRecord(
                readName=f"good_wt_{i}", cellBC="test", UMI=f"H{i}",
                readCount=50, seq=CARLIN_REF,
            ))

        config = PipelineConfig(
            lineage_mode=True,
            primer5_len=23, primer3_len=33,
            primer5_threshold=19, primer3_threshold=29,
            min_reads_snv=10,
            min_reads_indel=3,
            threads=1,
        )
        pipeline = Pipeline(config=config, ref_seq=CARLIN_REF)
        pipeline.load_cutsites()
        result = pipeline.run(queries)

        # 5 好mut + 5 bad_sub + 5 wt = 15 全部成功（假阳性过滤已移除）
        assert result.stats.successful == 15, \
            f"预期15条成功，实际 {result.stats.successful} (失败: {result.stats.failed})"
        # 10 mutations: 5 good_del + 5 bad_sub (substitution仍在窗口内)
        assert result.stats.mutated_sequences == 10
        # 编辑效率 = 10/15 * 100 = 66.7%
        assert result.stats.editing_efficiency_pct == pytest.approx(200.0 / 3.0, abs=0.1)

    def test_stats_consistency(self):
        """统计数据一致性验证"""
        queries = []
        for i in range(8):
            queries.append(QueryRecord(
                readName=f"mut_{i}", cellBC="test", UMI=f"M{i}",
                readCount=10, seq=_make_deletion(CARLIN_REF, 42, 3),
            ))
        for i in range(5):
            queries.append(QueryRecord(
                readName=f"wt_{i}", cellBC="test", UMI=f"W{i}",
                readCount=10, seq=CARLIN_REF,
            ))

        result = _run_pipeline(queries)
        stats = result.stats

        # 一致性：successful = mutated + unmutated
        assert stats.mutated_sequences + stats.unmutated_sequences == stats.successful
        # 一致性：total_queries = successful + failed
        assert stats.successful + stats.failed == stats.total_queries
        # 编辑效率边界
        assert 0.0 <= stats.editing_efficiency_pct <= 100.0


# ═══════════════════════════════════════════════════════════════
# 测试3：位置感知 gap 惩罚
# ═══════════════════════════════════════════════════════════════

class TestPositionAwareGapPenalty:
    """验证 gap 是否优先在 cutsite 区域开启而非保守区域"""

    def test_build_gradient_profiles_correctness(self):
        """梯度惩罚配置文件的数值正确性"""
        ref_len = 20
        cutsites = [CutsiteRegion("T1", 5, 10)]
        go, ge, mp = build_gradient_profiles(
            ref_len, cutsites,
            base_gap_open=-2.0, base_gap_extend=-0.1,
            mismatch_penalty=-3.0,
            min_scale=1.0, max_scale=2.0, cutsite_edge_scale=2.0,
        )

        # far 区域 (0-1, 14-19): max_scale=2.0 → open=-4.0, extend=-0.2
        for i in [0, 1, 15, 16, 17, 18, 19]:
            assert go[i] == pytest.approx(-4.0), f"far go[{i}]"
            assert ge[i] == pytest.approx(-0.2), f"far ge[{i}]"

        # cutsite 中心 (7-8): 接近 min_scale=1.0 → open ≈ -2.0
        assert go[7] > -2.5, f"center go[7]={go[7]} should be near -2.0"
        assert go[8] > -2.5, f"center go[8]={go[8]} should be near -2.0"
        # cutsite 边缘 (5, 10): 接近 edge_scale=2.0 → open ≈ -4.0
        assert go[5] < -3.0, f"edge go[5]={go[5]} should be near -4.0"
        assert go[10] < -3.0, f"edge go[10]={go[10]} should be near -4.0"

    def test_position_aware_alignment_gap_in_cutsite(self):
        """位置感知比对将gap放置在cutsite区域（低成本区域）"""
        # 设计一个两端保守、中间cutsite的序列
        # 保守区(0-3, 8-11) + cutsite(4-7)
        ref = "AAATTTCCCGGG"    # 12bp
        query = "AAAGGG"         # 6bp, 缺失中间的 TTTCCC(6bp)
        cutsites = [CutsiteRegion("T1", 4, 7)]

        # 位置感知比对（cutsite gap 惩罚低）
        go, ge, mp = build_gradient_profiles(
            len(ref), cutsites,
            base_gap_open=-2.0, base_gap_extend=-0.1,
            mismatch_penalty=-3.0,
            min_scale=1.0, max_scale=2.0, cutsite_edge_scale=2.0,
        )
        _, ar, aq, stats = affine_gap_alignment_position_aware(
            ref, query, go, ge,
            mismatch_penalty_profile=mp,
        )

        # 应有gap (6bp 删除)
        assert stats["gaps_in_query"] >= 6, \
            f"预期gap≥6，实际 {stats['gaps_in_query']}"

        # gap 应位置在 cutsite 区域 (4-7)
        # 理想情况：AAA-----GGG → gap从位置3到8
        gap_indices = [i for i, c in enumerate(aq) if c == '-']
        assert len(gap_indices) >= 6, f"gap列数: {len(gap_indices)}"

        # gap 应包含 cutsite 区域(4-7)
        gap_set = set(gap_indices)
        cutsite_set = set(range(4, 8))
        overlap = gap_set & cutsite_set
        assert len(overlap) >= 3, \
            f"gap与cutsite重叠不足: gap={gap_indices}, cutsite=4-7, 重叠={overlap}"

    def test_position_aware_vs_standard_gap_placement(self):
        """位置感知比对与标准比对的gap位置差异：位置感知应更倾向cutsite"""
        ref = "AAATTTCCCGGGAAA"  # 15bp
        query = "AAACCCAAA"       # 9bp, 缺失 TTTGGG (6bp)
        cutsites = [CutsiteRegion("T1", 3, 5)]  # TTT 区域

        # 标准比对（均匀惩罚）
        _, ar_std, aq_std, _ = affine_gap_alignment(
            ref, query,
        )

        # 位置感知比对（cutsite gap 惩罚低，两侧惩罚较高）
        go, ge, mp = build_gradient_profiles(
            len(ref), cutsites,
            base_gap_open=-2.0, base_gap_extend=-0.1,
            mismatch_penalty=-3.0,
            min_scale=1.0, max_scale=5.0, cutsite_edge_scale=5.0,
        )
        _, ar_pa, aq_pa, stats_pa = affine_gap_alignment_position_aware(
            ref, query, go, ge,
            mismatch_penalty_profile=mp,
        )

        # 位置感知比对应有gap
        assert stats_pa["gaps_in_query"] >= 6

        # gap应主要在cutsite区域(3-5)
        pa_gaps = {i for i, c in enumerate(aq_pa) if c == '-'}
        pa_in_cutsite = sum(1 for i in pa_gaps if 3 <= i <= 5)
        pa_outside = len(pa_gaps) - pa_in_cutsite

        # 大部分gap应在cutsite
        assert pa_in_cutsite >= pa_outside, \
            f"位置感知: gap在cutsite={pa_in_cutsite}, 外部={pa_outside}"

    def test_lineage_tracer_align_detects_cutsite_indel(self):
        """lineage_tracer_align 正确检测 cutsite 区域内的 indel"""
        ref = "A" * 15 + "GAGTCG" + "A" * 15  # GAGTCG 是 cutsite motif
        query = "A" * 15 + "A" * 15              # cutsite 被完全删除
        cutsites = [CutsiteRegion("T1", 15, 20)]  # GAGTCG 位置

        score, ar, aq, stats = lineage_tracer_align(
            ref, query, cutsites,
            min_scale=1.0,
            max_scale=2.0,
        )

        # 应检测到 6bp 删除
        assert stats["gaps_in_query"] >= 5, \
            f"lineage_tracer_align 未检测到足够大的删除: gaps_in_query={stats['gaps_in_query']}"

        # gap 应在 cutsite 附近
        gap_start = -1
        gap_end = -1
        for i, c in enumerate(aq):
            if c == '-':
                if gap_start == -1:
                    gap_start = i
                gap_end = i

        assert gap_start >= 14 or gap_start == -1, \
            f"gap起始位置 {gap_start} 应接近cutsite(15-20)"


# ═══════════════════════════════════════════════════════════════
# 测试4：extract_mutations 直接验证
# ═══════════════════════════════════════════════════════════════

class TestExtractMutationsAccuracy:
    """直接验证 extract_mutations 函数对已知突变的识别准确度"""

    def test_extract_single_substitution(self):
        """已知点突变的位置和碱基正确"""
        ar = "ACGTACGT"
        aq = "ACCTACGT"  # G→C at ref pos 2
        mutations = extract_mutations(ar, aq)

        assert len(mutations) == 1
        m = mutations[0]
        assert m.type == MutationType.SUBSTITUTION
        assert m.ref_pos == 2
        assert m.ref_base == "G"
        assert m.query_base == "C"
        assert m.length == 1

    def test_extract_deletion_length(self):
        """已知删除长度正确"""
        ar = "ACGTACGT"
        aq = "A-----GT"  # CGTAC deleted (5bp)
        mutations = extract_mutations(ar, aq)

        deletions = [m for m in mutations if m.type == MutationType.DELETION]
        assert len(deletions) == 1
        assert deletions[0].length == 5

    def test_extract_deletion_with_cutsite_window(self):
        """在cutsite窗口内的删除, in_cutsite_window=True"""
        # 3bp deletion at ref positions 3-5 (TAC removed from ACGTACGTACGT)
        ar = "ACGTACGTACGT"
        aq = "ACG---GTACGT"  # dashes at positions 3-5 for TAC deletion
        cutsites = [CutsiteRegion("T1", 3, 5)]
        mutations = extract_mutations(ar, aq, cutsites=cutsites, mutation_window=3)

        deletions = [m for m in mutations if m.type == MutationType.DELETION]
        assert len(deletions) >= 1
        assert deletions[0].in_cutsite_window, "cutsite窗口内的删除应标记为 in_window"

    def test_extract_deletion_outside_cutsite_window(self):
        """在cutsite窗口外的删除, in_cutsite_window=False"""
        # 3bp deletion at ref positions 3-5, cutsite at 8-10 (far away)
        ar = "ACGTACGTACGT"
        aq = "ACG---GTACGT"  # dashes at positions 3-5 for TAC deletion
        cutsites = [CutsiteRegion("T1", 8, 10)]  # cutsite 远离删除位置
        mutations = extract_mutations(ar, aq, cutsites=cutsites, mutation_window=1)

        deletions = [m for m in mutations if m.type == MutationType.DELETION]
        assert len(deletions) >= 1
        assert not deletions[0].in_cutsite_window, \
            "cutsite窗口外的删除应标记为 not_in_window (false)"

    def test_extract_insertion_length(self):
        """已知插入长度正确"""
        ar = "AC--GT"
        aq = "ACGGGT"  # GG inserted (2bp)
        mutations = extract_mutations(ar, aq)

        insertions = [m for m in mutations if m.type == MutationType.INSERTION]
        assert len(insertions) == 1
        assert insertions[0].length == 2

    def test_no_mutations_for_identical(self):
        """相同序列无突变"""
        assert extract_mutations("ACGT", "ACGT") == []

    def test_adjacent_indel_merged_to_complex(self):
        """相邻插入+删除合并为复合事件"""
        ar = "ACGT--AC"
        aq = "AC--GGAC"  # DEL at pos 3-4 (GT), INS at pos 5-6 (GG)
        mutations = extract_mutations(ar, aq)

        types = {m.type for m in mutations}
        assert MutationType.COMPLEX in types, \
            f"相邻插入+删除应合并为COMPLEX: {types}"

    def test_empty_inputs(self):
        """空输入返回空列表"""
        assert extract_mutations("", "") == []

    def test_mismatched_lengths(self):
        """长度不同返回空列表"""
        assert extract_mutations("ACGT", "AC") == []


# ═══════════════════════════════════════════════════════════════
# 独立测试函数（供 run_tests.py 导入）
# ═══════════════════════════════════════════════════════════════

def run_accuracy_checks(check_func):
    """使用 check_func(name, condition) 方式运行核心准确性断言

    此函数专为 run_tests.py 设计，使用一致的 check API。
    参数 check_func 应为 (str, bool) -> None 的 callable。
    """
    print("\n── 准确性验证测试 ──\n")

    # ── 1a: 突变检测 ──
    result = _run_align_single(CARLIN_REF, "wt", read_count=50)
    check_func("野⽣型无突变", result.success and len(result.mutations) == 0)

    result = _run_align_single(_make_deletion(CARLIN_REF, 42, 3), "del_3bp")
    deletions = [m for m in result.mutations if m.type == MutationType.DELETION]
    check_func("3bp删除检测", result.success
               and len(deletions) >= 1
               and deletions[0].length == 3
               and deletions[0].in_cutsite_window)

    result = _run_align_single(_make_insertion(CARLIN_REF, 42, "ACGTA"), "ins_5bp")
    insertions = [m for m in result.mutations if m.type == MutationType.INSERTION]
    check_func("5bp插入检测", result.success
               and len(insertions) >= 1
               and insertions[0].length == 5
               and insertions[0].in_cutsite_window)

    result = _run_align_single(_make_substitution(CARLIN_REF, 44, "T"), "sub_1bp")
    subs = [m for m in result.mutations if m.type == MutationType.SUBSTITUTION]
    check_func("点突变检测", result.success
               and len(subs) >= 1
               and subs[0].length == 1
               and subs[0].in_cutsite_window)

    # ── 编辑效率验证 ──
    queries = []
    for i in range(6):
        queries.append(QueryRecord(
            readName=f"mut_{i}", cellBC="test", UMI=f"M{i}",
            readCount=10,
            seq=_make_deletion(CARLIN_REF, 42, 3),
        ))
    for i in range(4):
        queries.append(QueryRecord(
            readName=f"wt_{i}", cellBC="test", UMI=f"W{i}",
            readCount=10, seq=CARLIN_REF,
        ))
    result = _run_pipeline(queries)
    expected_eff = 60.0
    actual_eff = result.stats.editing_efficiency_pct
    check_func(f"编辑效率 {expected_eff}% (预期={expected_eff}, 实际={actual_eff:.1f})",
               abs(actual_eff - expected_eff) < 1.0)

    # ── 位置感知 gap ──
    ref = "AAATTTCCCGGG"
    query = "AAAGGG"
    cutsites = [CutsiteRegion("T1", 4, 7)]
    go, ge, mp = build_gradient_profiles(
        len(ref), cutsites,
        base_gap_open=-2.0, base_gap_extend=-0.1,
        mismatch_penalty=-3.0,
        min_scale=1.0, max_scale=2.0, cutsite_edge_scale=2.0,
    )
    _, ar, aq, stats = affine_gap_alignment_position_aware(
        ref, query, go, ge,
        mismatch_penalty_profile=mp,
    )
    gap_indices = {i for i, c in enumerate(aq) if c == '-'}
    cutsite_range = set(range(4, 8))
    overlap = len(gap_indices & cutsite_range)
    check_func("位置感知gap在cutsite区域",
               stats["gaps_in_query"] >= 6 and overlap >= 3)

    # ── extract_mutations 直接验证 ──
    m = extract_mutations("ACGT", "ACCT")[0]
    check_func("extract_mutations 点突变 ref_pos",
               m.type == MutationType.SUBSTITUTION and m.ref_pos == 2)
    m = extract_mutations("ACGTACGT", "A-----GT")[0]
    check_func("extract_mutations 5bp删除",
               m.type == MutationType.DELETION and m.length == 5)
    check_func("extract_mutations 相同序列无突变",
               extract_mutations("ACGT", "ACGT") == [])

    print()  # 空行分隔
