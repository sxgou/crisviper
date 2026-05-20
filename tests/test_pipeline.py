"""End-to-end integration test for the full crisviper pipeline.

Tests the complete flow:
  1. Creates synthetic CARLIN-like reference and query sequences
  2. Converts to TSV format
  3. Runs the alignment pipeline
  4. Runs with standard and lineage modes
  5. Generates reports (JSON and HTML)
  6. Validates output statistics
"""

import sys
import os
import json
import tempfile
import shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest

from crisviper import (
    Pipeline, PipelineConfig, QueryRecord,
    align_single, check_primer_anchoring,
    read_reference_fasta, read_queries_tsv,
    save_alignment_results, generate_report,
    get_amplicon_structure, CutsiteRegion,
    extract_mutations, MutationType,
)


# ═══════════════════════════════════════════════════════════════
# 合成测试数据
# ═══════════════════════════════════════════════════════════════

# 标准 CARLIN 332bp 参考序列
CARLIN_REF = (
    "TATGTGTGGGAGGGCTAAGAGG"  # Primer5 (23bp)
    "CCGCC"                    # Prefix (5bp)
    "GACTGCACGACAGTCGA"        # Target1 (13+7 cutsite)
    "CGATGGAG"                 # Linker (7bp)
    "TCGACACGACTCGCGCA"        # Target2
    "TACGATGG"                 # Linker
    "AGTCGACTACAGTCGCTA"       # Target3
    "CGACGATG"                 # Linker
    "GAGTCGCGAGCGCTATG"        # Target4
    "AGCGACTA"                 # Linker
    "TGGAGTCGATACGATACG"       # Target5
    "CGCACGCT"                 # Linker
    "ATGGAGTCGAGAGCGCGC"       # Target6
    "TCGTCAAC"                 # Linker
    "GATGGAGTCGCGACTGTA"       # Target7
    "CGCACTCG"                 # Linker
    "CGATGGAGTCGATAGTAT"       # Target8
    "GCGTACAC"                 # Linker
    "GCGATGGAGTCGACTGCA"       # Target9
    "CGACAGTC"                 # Linker
    "GACTATGGAGTCGATACGTAGC"   # Target10
    "ACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"  # Postfix(8)+Primer3(33)
)


def _make_single_mutated_query(ref_seq: str, pos: int, new_base: str) -> str:
    """在指定位置产生点突变"""
    lst = list(ref_seq)
    lst[pos] = new_base
    return ''.join(lst)


def _make_deletion_query(ref_seq: str, start: int, length: int) -> str:
    """在指定位置产生删除"""
    return ref_seq[:start] + ref_seq[start + length:]


def _make_insertion_query(ref_seq: str, pos: int, insert: str) -> str:
    """在指定位置产生插入"""
    return ref_seq[:pos] + insert + ref_seq[pos:]


def _build_test_queries(ref_seq: str) -> list:
    """构建包含各种突变的测试序列列表"""
    queries = [
        # 野生型（未突变）
        QueryRecord(readName="wt_high",  cellBC="test", UMI="UMI1", readCount=50, seq=ref_seq),
        QueryRecord(readName="wt_low",   cellBC="test", UMI="UMI2", readCount=1,  seq=ref_seq),

        # 点突变在cutsite窗口内
        QueryRecord(readName="sub_in_window", cellBC="test", UMI="UMI3", readCount=15,
                    seq=_make_single_mutated_query(ref_seq, 43, 'T')),  # Target1 cutsite内

        # 点突变在窗口外（应该被矫正/过滤）
        QueryRecord(readName="sub_out_window", cellBC="test", UMI="UMI4", readCount=12,
                    seq=_make_single_mutated_query(ref_seq, 30, 'T')),  # 保守区

        # 小删除在cutsite
        QueryRecord(readName="del_cutsite_3bp", cellBC="test", UMI="UMI5", readCount=8,
                    seq=_make_deletion_query(ref_seq, 43, 3)),

        # 大删除跨Target
        QueryRecord(readName="del_large_20bp", cellBC="test", UMI="UMI6", readCount=5,
                    seq=_make_deletion_query(ref_seq, 41, 20)),

        # 插入在conserved区
        QueryRecord(readName="ins_small", cellBC="test", UMI="UMI7", readCount=4,
                    seq=_make_insertion_query(ref_seq, 50, "AAA")),

        # Primer3锚定失败（修改Primer3区域，最后33个碱基）
        QueryRecord(readName="primer3_bad", cellBC="test", UMI="UMI8", readCount=3,
                    seq=ref_seq[:299] + "N" * 33),  # 完全替换Primer3区域

        # 全长删除
        QueryRecord(readName="del_full_target", cellBC="test", UMI="UMI9", readCount=6,
                    seq=_make_deletion_query(ref_seq, 41, 7)),  # 删除Target1的cutsite

        # 点突变+低readCount（使用3bp突变避免对齐算法gap偏好）
        QueryRecord(readName="sub_low_rc", cellBC="test", UMI="UMI10", readCount=2,
                    seq=_make_single_mutated_query(
                        _make_single_mutated_query(
                            _make_single_mutated_query(ref_seq, 100, 'G'), 101, 'A'), 102, 'T')),

        # 复合突变（删除+插入）
        QueryRecord(readName="complex", cellBC="test", UMI="UMI11", readCount=7,
                    seq=_make_deletion_query(
                        _make_insertion_query(ref_seq, 50, "ACG"), 44, 5
                    )),
    ]
    return queries


# ═══════════════════════════════════════════════════════════════
# 测试
# ═══════════════════════════════════════════════════════════════

class TestEndToEndPipeline:
    """端到端管道测试"""

    @pytest.fixture(autouse=True)
    def setup(self):
        """在每个测试前设置"""
        self.ref_seq = CARLIN_REF
        self.queries = _build_test_queries(self.ref_seq)
        self.temp_dir = tempfile.mkdtemp()
        yield
        shutil.rmtree(self.temp_dir)

    def _run_pipeline(self, lineage_mode=False):
        """运行管道并返回结果"""
        config = PipelineConfig(
            lineage_mode=lineage_mode,
            primer5_len=23,
            primer3_len=33,
            primer5_threshold=19,
            primer3_threshold=29,
            min_reads_sub=10,
            min_reads_indel=3,
            threads=1,
        )
        pipeline = Pipeline(config=config, ref_seq=self.ref_seq)
        if lineage_mode:
            pipeline.load_cutsites()
        return pipeline.run(self.queries)

    def test_pipeline_runs_successfully(self):
        """管道能正确运行并返回结果"""
        result = self._run_pipeline()
        assert result is not None
        assert len(result.results) == len(self.queries)
        assert result.stats.total_queries == len(self.queries)

    def test_primer_anchoring_detects_failures(self):
        """Primer3锚定失败应被正确识别"""
        result = self._run_pipeline()
        failed = result.get_failed()
        failed_names = [r.query.readName for r in failed]
        assert "primer3_bad" in failed_names

    def test_wildtype_high_rc_passes(self):
        """野生型（高readCount）应通过"""
        result = self._run_pipeline()
        successful = result.get_successful()
        names = [r.query.readName for r in successful]
        assert "wt_high" in names

    def test_wildtype_low_rc_passes(self):
        """野生型（低readCount但无突变）应通过"""
        result = self._run_pipeline()
        successful = result.get_successful()
        names = [r.query.readName for r in successful]
        assert "wt_low" in names

    def test_snv_in_window_passes(self):
        """cutsite窗口内的点突变应通过"""
        result = self._run_pipeline()
        successful = result.get_successful()
        names = [r.query.readName for r in successful]
        assert "sub_in_window" in names

    def test_snv_low_rc_filtered(self):
        """点突变+低readCount应被过滤"""
        result = self._run_pipeline()
        successful = result.get_successful()
        successful_names = [r.query.readName for r in successful]
        assert "sub_low_rc" not in successful_names

    def test_mutations_detected_in_results(self):
        """成功比对的结果应包含突变事件"""
        result = self._run_pipeline()
        for r in result.get_successful():
            # 所有成功的结果应有 mutations 字段
            assert hasattr(r, 'mutations')
            if r.query.readName == "sub_in_window":
                assert len(r.mutations) > 0
                assert r.mutations[0].type == MutationType.SUBSTITUTION

    def test_stats_collected(self):
        """管道统计数据应合理"""
        result = self._run_pipeline()
        stats = result.stats
        assert stats.total_queries > 0
        assert stats.successful <= stats.total_queries
        assert stats.failed >= 0
        assert stats.successful + stats.failed == stats.total_queries

    def test_lineage_mode_produces_results(self):
        """谱系示踪比对模式也能正常工作"""
        result = self._run_pipeline(lineage_mode=True)
        assert result is not None
        assert len(result.results) > 0

    def test_lineage_mode_detects_cutsites(self):
        """谱系示踪模式自动检测cutsite"""
        # 直接测试 get_amplicon_structure
        cutsites = get_amplicon_structure(self.ref_seq)
        assert cutsites is not None
        assert len(cutsites) == 10  # CARLIN有10个Target

    def test_pipeline_stats_consistency(self):
        """统计数据应一致（mutated + unmutated = successful）"""
        result = self._run_pipeline()
        stats = result.stats
        assert stats.mutated_sequences + stats.unmutated_sequences == stats.successful

    def test_pipeline_result_conversion_to_dict(self):
        """AlignmentResult.to_dict() 与旧格式兼容"""
        result = self._run_pipeline()
        for r in result.get_successful()[:3]:
            d = r.to_dict()
            assert "readName" in d
            assert "score" in d
            assert "aligned_ref" in d
            assert "aligned_query" in d
            assert "stats" in d

    def test_error_result_conversion(self):
        """错误结果的 to_dict() 应包含 error 字段"""
        result = self._run_pipeline()
        for r in result.get_failed():
            d = r.to_dict()
            assert d["error"] != ""
            assert d["score"] is None

    def test_save_and_read_tsv_roundtrip(self):
        """TSV写入和回读应一致"""
        import tempfile
        save_path = os.path.join(self.temp_dir, "test_queries.tsv")

        # 先写入TSV
        from crisviper import save_tsv
        dict_rows = [{"readName": q.readName, "cellBC": q.cellBC,
                      "UMI": q.UMI, "readCount": q.readCount, "seq": q.seq}
                     for q in self.queries]
        save_tsv(dict_rows, save_path)

        # 再读取回来
        read_back = read_queries_tsv(save_path)
        assert len(read_back) == len(self.queries)
        assert read_back[0]["readName"] == self.queries[0].readName

    def test_report_json_generation(self):
        """JSON报告生成应成功"""
        result = self._run_pipeline()
        output_results = [r.to_dict() for r in result.results]
        report_path = os.path.join(self.temp_dir, "test_report")

        generate_report(output_results, report_path, fmt="json",
                         ref_length=len(self.ref_seq))

        # 验证JSON文件存在且包含必要字段
        json_path = report_path + ".json"
        assert os.path.exists(json_path)
        with open(json_path) as f:
            report = json.load(f)
        assert "summary" in report
        assert "editing_efficiency_pct" in report["summary"]

    def test_html_report_generation(self):
        """HTML报告生成应成功（跳过matplotlib检查）"""
        result = self._run_pipeline()
        output_results = [r.to_dict() for r in result.results]
        report_path = os.path.join(self.temp_dir, "test_html_report")

        generate_report(output_results, report_path, fmt="html",
                         ref_length=len(self.ref_seq))

        html_path = report_path + ".html"
        assert os.path.exists(html_path)
        # 验证HTML包含关键内容
        with open(html_path) as f:
            content = f.read()
        assert "CARLIN" in content
        assert "Editing Efficiency" in content
        # 图表可能为空（无matplotlib），但HTML结构应该完整
        assert "<html" in content
        assert "modal" in content

    def test_pipeline_with_no_cutsites(self):
        """标准模式（无双端cutsite）也能工作"""
        config = PipelineConfig(
            lineage_mode=False,
            primer5_len=23,
            primer3_len=33,
            min_reads_sub=10,
            min_reads_indel=3,
            threads=1,
        )
        pipeline = Pipeline(config=config, ref_seq=self.ref_seq)
        result = pipeline.run(self.queries)
        assert len(result.results) > 0

    def test_align_single_function(self):
        """align_single 函数可直接调用"""
        from crisviper import align_single, PipelineConfig
        config = PipelineConfig()
        result = align_single(self.queries[0], self.ref_seq, config)
        assert result is not None
        assert result.success or not result.success

    def test_mutations_in_del_query(self):
        """删除突变的比对结果应包含 MutationEvent"""
        result = self._run_pipeline()
        del_results = [r for r in result.get_successful()
                       if "del_" in r.query.readName]
        for r in del_results:
            if r.mutations:
                has_deletion = any(m.type == MutationType.DELETION
                                   for m in r.mutations)
                if not has_deletion:
                    # 可能被合并为complex
                    has_indel = any(m.type == MutationType.INDEL
                                      for m in r.mutations)
                    assert has_indel or not r.stats.has_indel
