"""crisviper CLI — 谱系示踪序列分析命令行工具

基于仿射gap惩罚算法的序列比对工具，支持FASTQ转换和并行批量比对。

子命令:
  convert   将FASTQ文件转换为TSV或FASTA格式
  align     并行批量将TSV或FASTA格式中的序列与reference进行比对

示例:
  $ crisviper convert fastq-to-tsv --fastq reads.fastq.gz --output reads.tsv
  $ crisviper align --reference ref.fasta --queries queries.tsv --output alignments.json
  $ crisviper align --reference ref.fasta --queries queries.tsv --output results/prefix --format all --report html
"""

import argparse
import sys
import os
import time
# 防止 fork + NumPy 线程冲突（在 import crisviper 之前设置）
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
    get_amplicon_structure,
    get_logger, setup_logging,
)

log = get_logger(__name__)

__version__ = "1.1.0"


def _queries_to_records(queries: List[Dict]) -> List[QueryRecord]:
    """将旧版 dict 格式转为 QueryRecord"""
    return [QueryRecord(**q) for q in queries]


def _log_timing(logger, label: str, t_start: float) -> float:
    """记录自 t_start 以来的耗时(s)，返回当前时间。"""
    elapsed = time.perf_counter() - t_start
    logger.info("  ⏱ %s: %.1fs", label, elapsed)
    return time.perf_counter()


def main():
    setup_logging(quiet=False)
    parser = argparse.ArgumentParser(
        description="crisviper — 谱系示踪序列分析命令行工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-v', '--verbose', action='store_true', help='显示调试信息')
    parser.add_argument('-q', '--quiet', action='store_true', help='仅显示警告和错误')

    subparsers = parser.add_subparsers(dest="command", required=True, help="命令")

    # ── 子命令: convert ──────────────────────────────────────────
    convert_parser = subparsers.add_parser("convert", help="转换FASTQ文件格式")
    convert_subparsers = convert_parser.add_subparsers(dest="convert_command", required=True, help="转换子命令")

    # fastq-to-tsv
    tsv_parser = convert_subparsers.add_parser("fastq-to-tsv", help="将FASTQ转换为TSV格式")
    tsv_parser.add_argument("--fastq", default=None, help="输入FASTQ文件（单端模式，支持.gz）")
    tsv_parser.add_argument("--fastq1", default=None, help="双端模式R1 FASTQ文件（配合 --fastq2）")
    tsv_parser.add_argument("--fastq2", default=None, help="双端模式R2 FASTQ文件（配合 --fastq1）")
    tsv_parser.add_argument("--output", required=True, help="输出TSV文件路径")
    tsv_parser.add_argument("--sample-name", default="sample", help="样本名称")
    tsv_parser.add_argument("--min-reads", type=int, default=1,
                            help="最小read数阈值（默认：1，即不过滤）")
    # 双端合并参数
    tsv_parser.add_argument("--min-overlap", type=int, default=10,
                            help="双端合并最小overlap长度 (bp)，默认10")
    tsv_parser.add_argument("--max-mismatch-rate", type=int, default=20,
                            help="双端合并overlap错配比例上限 (%%)，默认20")
    tsv_parser.add_argument("--max-mismatch-diff", type=int, default=5,
                            help="双端合并overlap最大错配绝对数，默认5")
    tsv_parser.add_argument("--require-qual", type=int, default=15,
                            help="双端合并碱基质量阈值(phred)，默认15")

    # fastq-to-fasta
    fasta_parser = convert_subparsers.add_parser("fastq-to-fasta", help="将FASTQ转换为FASTA格式")
    fasta_parser.add_argument("--fastq", default=None, help="输入FASTQ文件（单端模式，支持.gz）")
    fasta_parser.add_argument("--fastq1", default=None, help="双端模式R1 FASTQ文件（配合 --fastq2）")
    fasta_parser.add_argument("--fastq2", default=None, help="双端模式R2 FASTQ文件（配合 --fastq1）")
    fasta_parser.add_argument("--output", required=True, help="输出FASTA文件路径")
    fasta_parser.add_argument("--sample-name", default="sample", help="样本名称")
    # 双端合并参数
    fasta_parser.add_argument("--min-overlap", type=int, default=10,
                              help="双端合并最小overlap长度 (bp)，默认10")
    fasta_parser.add_argument("--max-mismatch-rate", type=int, default=20,
                              help="双端合并overlap错配比例上限 (%%)，默认20")
    fasta_parser.add_argument("--max-mismatch-diff", type=int, default=5,
                              help="双端合并overlap最大错配绝对数，默认5")
    fasta_parser.add_argument("--require-qual", type=int, default=15,
                              help="双端合并碱基质量阈值(phred)，默认15")

    # ── 子命令: align（并行批量比对） ──────────────────────────────
    align_parser = subparsers.add_parser("align", help="并行批量序列比对")
    align_parser.add_argument("--reference", required=True, help="reference序列FASTA文件")
    align_parser.add_argument("--queries", required=True, help="查询序列文件（TSV或FASTA格式）")
    align_parser.add_argument("--output", required=True, help="输出文件路径（使用--format all时作为前缀）")
    align_parser.add_argument("--format", choices=["json", "tsv", "all"], default="json",
                              help="输出格式: json（默认）、tsv、all")

    # 比对参数
    align_parser.add_argument("--match-score", type=float, default=2.0)
    align_parser.add_argument("--mismatch-penalty", type=float, default=-3.0)
    align_parser.add_argument("--gap-open", type=float, default=-2.0)
    align_parser.add_argument("--gap-extend", type=float, default=-0.1)
    # 谱系示踪模式
    align_parser.add_argument("--lineage", action="store_true",
                              help="启用谱系示踪比对模式")
    align_parser.add_argument("--min-scale", type=float, default=1.0,
                              help="切割点处最低惩罚倍率（默认：1.0）")
    align_parser.add_argument("--max-scale", type=float, default=6.0,
                              help="保守区最高惩罚倍率（默认：6.0）")
    align_parser.add_argument("--cutsite-edge-scale", type=float, default=2.0,
                              help="Cutsite边界惩罚倍率（默认：2.0）")
    align_parser.add_argument("--gradient-radius", type=float, default=None,
                              help="梯度有效半径 (bp)，省略则自动计算")
    align_parser.add_argument("--sub-window", type=int, default=3,
                              help="cutsite邻近保留窗口(bp)，控制背景矫正和突变标注")
    align_parser.add_argument("--mismatch-density-threshold", type=float, default=0.34,
                              help="密集错配检测密度阈值")
    align_parser.add_argument("--cutsites", default=None,
                              help="cutsite位置配置文件（JSON格式）")

    align_parser.add_argument("--gap-exit-strength", type=float, default=0.0,
                              help="Gap exit抑制强度（≤0，0=关闭，推荐-3.5在cutsite中心产生约-21.0峰值惩罚）")
    align_parser.add_argument("--short-match-window", type=int, default=0,
                              help="短匹配区域阈值（bp，0=关闭，推荐3~5）")
    align_parser.add_argument("--short-match-discount", type=float, default=1.0,
                              help="短匹配区域match_score折扣系数（0~1，1.0=不打折，推荐0.5）")
    align_parser.add_argument("--dense-mismatch-window", type=int, default=6,
                              help="密集错配检测窗口大小（bp）")
    align_parser.add_argument("--dense-mismatch-penalty", type=float, default=0.0,
                              help="密集错配区域额外惩罚（≤0，0=关闭，推荐-2.0）")
    align_parser.add_argument("--homology-window", type=int, default=8,
                              help="同源区域检测窗口大小（bp）")
    align_parser.add_argument("--homology-penalty", type=float, default=0.0,
                              help="同源区域match_score惩罚（≤0，0=关闭，推荐-1.0）")
    align_parser.add_argument("--isolated-base-penalty", type=float, default=0.0,
                              help="孤立碱基匹配额外惩罚（≤0，0=关闭，推荐-2.0），吸收孤立匹配到gap端点")

    # 并行参数
    align_parser.add_argument("--threads", "-t", type=int, default=1)
    align_parser.add_argument("--chunk-size", type=int, default=None,
                              help="每批次序列数（默认：自动计算，约 total/(threads×3)）")

    # 引物参数
    align_parser.add_argument("--primer5-len", type=int, default=23)
    align_parser.add_argument("--primer3-len", type=int, default=33)
    align_parser.add_argument("--primer5-threshold", type=int, default=19)
    align_parser.add_argument("--primer3-threshold", type=int, default=29)

    # Allele过滤参数
    align_parser.add_argument("--min-reads", type=int, default=1,
                              help="输入侧最小read数阈值（默认：1，即不过滤，仅用于预过滤）")
    align_parser.add_argument("--min-reads-sub", type=int, default=5,
                              help="纯点突变allele最小read数阈值（exclusive，>此值通过，默认5）")
    align_parser.add_argument("--min-reads-indel", type=int, default=0,
                              help="含indel的allele最小read数阈值（0=不过滤）")

    # 背景点突变矫正参数
    align_parser.add_argument("--correct-bg-sub", action=argparse.BooleanOptionalAction,
                              default=True, help="启用背景点突变矫正（默认开启）")
    align_parser.add_argument("--keep-sub-indel-window", type=int, default=3,
                              help="背景矫正时indel邻近保留窗口(bp)")

    # 报告参数
    align_parser.add_argument("--report", choices=["json", "html"], default=None)
    align_parser.add_argument("--report-output", default=None)
    align_parser.add_argument("--allele-window-start", type=int, default=0)
    align_parser.add_argument("--allele-window-end", type=int, default=None)
    align_parser.add_argument("--allele-top-n", type=int, default=50)

    args = parser.parse_args()
    setup_logging(verbose=args.verbose, quiet=args.quiet)

    # ═══════════════════════════════════════════════════════════════
    # 执行命令
    # ═══════════════════════════════════════════════════════════════
    if args.command == "convert":
        if args.convert_command == "fastq-to-tsv":
            # 双端模式
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
                log.error("请指定 --fastq（单端）或 --fastq1 + --fastq2（双端）")
                sys.exit(1)

            if args.min_reads > 1:
                before = len(rows)
                rows = [r for r in rows if r['readCount'] >= args.min_reads]
                log.info("Allele过滤 (min_reads>=%d): %d → %d", args.min_reads, before, len(rows))
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
                log.error("请指定 --fastq（单端）或 --fastq1 + --fastq2（双端）")
                sys.exit(1)

    elif args.command == "align":
        t0 = time.perf_counter()

        # 读取reference
        log.info("读取reference序列: %s", args.reference)
        ref_seq = read_reference_fasta(args.reference)
        log.info("reference序列长度: %d bp", len(ref_seq))
        t = _log_timing(log, "读取reference", t0)

        # 读取查询序列
        queries_path = args.queries
        if queries_path.endswith('.tsv'):
            queries = read_queries_tsv(queries_path)
        elif queries_path.endswith(('.fasta', '.fa', '.fas')):
            queries = read_queries_fasta(queries_path)
        else:
            log.error("不支持的查询文件格式: %s", queries_path)
            sys.exit(1)
        t = _log_timing(log, "读取查询序列", t)

        # 转为类型安全的 QueryRecord
        query_records = _queries_to_records(queries)

        # Allele质量过滤（转换前）
        if args.min_reads > 1:
            before = len(query_records)
            query_records = [q for q in query_records if q.readCount >= args.min_reads]
            log.info("Allele过滤 (min_reads>=%d): %d → %d",
                     args.min_reads, before, len(query_records))

        # ── 计算批次大小 ──
        total_queries = len(query_records)
        if args.chunk_size is not None:
            chunk_size = args.chunk_size
        else:
            target_chunks = max(args.threads * 3, 12)
            chunk_size = max(100, total_queries // target_chunks)
        log.info("  批次大小: %d 条/批（共 %d 批）", chunk_size,
                 (total_queries + chunk_size - 1) // chunk_size)

        # ── 构建 PipelineConfig ──
        config = PipelineConfig(
            match_score=args.match_score,
            mismatch_penalty=args.mismatch_penalty,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend,
            lineage_mode=args.lineage,
            min_scale=args.min_scale,
            max_scale=args.max_scale,
            cutsite_edge_scale=args.cutsite_edge_scale,
            gradient_radius=args.gradient_radius,
            mismatch_density_threshold=args.mismatch_density_threshold,
            sub_window=args.sub_window,
            primer5_len=args.primer5_len,
            primer3_len=args.primer3_len,
            primer5_threshold=args.primer5_threshold,
            primer3_threshold=args.primer3_threshold,
            gap_exit_strength=args.gap_exit_strength,
            short_match_window=args.short_match_window,
            short_match_discount=args.short_match_discount,
            dense_mismatch_window=args.dense_mismatch_window,
            dense_mismatch_penalty=args.dense_mismatch_penalty,
            homology_window=args.homology_window,
            homology_penalty=args.homology_penalty,
            isolated_base_penalty=args.isolated_base_penalty,
            min_reads_sub=args.min_reads_sub,
            min_reads_indel=args.min_reads_indel,
            threads=args.threads,
            chunk_size=chunk_size,
            cutsites_path=args.cutsites,
            report_format=args.report,
            allele_top_n=args.allele_top_n,
            allele_window_start=args.allele_window_start,
            allele_window_end=args.allele_window_end,
        )

        # 显示配置信息
        radius_info = "auto" if args.gradient_radius is None else f"{args.gradient_radius}bp"
        log.info("  梯度惩罚参数: min_scale=%s, max_scale=%s, edge_scale=%s, radius=%s",
                 args.min_scale, args.max_scale, args.cutsite_edge_scale, radius_info)
        if args.lineage:
            log.info("  sub窗口: cutsite±%s bp, mismatch密度阈值: %s",
                     args.sub_window, args.mismatch_density_threshold)

        if args.threads > 1:
            log.info("开始并行批量比对（全局，%d 线程）...", args.threads)
        else:
            log.info("开始批量比对（全局，单线程）...")
        # ── 运行管道 ──
        pipeline = Pipeline(config=config, ref_seq=ref_seq)
        pipeline_result = pipeline.run(query_records)
        t = _log_timing(log, "谱系示踪比对（并行）", t)

        # ── 转换为兼容的输出格式 ──
        output_results = [r.to_dict() for r in pipeline_result.results]
        t = _log_timing(log, "结果转dict", t)

        # ── 保存结果 ──
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        save_alignment_results(output_results, args.output, args.format)
        t = _log_timing(log, "保存TSV结果", t)

        # ── 保存总结表格 ──
        output_dir = os.path.dirname(os.path.abspath(args.output))
        save_summary_tables(output_results, output_dir, ref_seq=ref_seq,
                            cutsites=_get_display_cutsites(ref_seq))
        t = _log_timing(log, "保存总结表格", t)

        # ── 生成分析报告 ──
        if args.report:
            report_output = args.report_output or _default_report_path(args.output)
            log.info("生成%s格式分析报告...", args.report.upper())
            generate_report(output_results, report_output, args.report,
                             ref_length=len(ref_seq),
                             ref_seq=ref_seq,
                             cutsites=_get_display_cutsites(ref_seq),
                             allele_window_start=args.allele_window_start,
                             allele_window_end=args.allele_window_end,
                             allele_top_n=args.allele_top_n,
                             version=__version__)
            t = _log_timing(log, "生成HTML报告", t)

        total = time.perf_counter() - t0
        log.info("  ⏱ 总耗时: %.1fs (%.1f min)", total, total / 60)

    else:
        parser.print_help()
        sys.exit(1)


def _default_report_path(output_path: str) -> str:
    """从--output参数推断默认报告路径"""
    base = output_path
    for ext in [".json", ".tsv"]:
        if base.endswith(ext):
            base = base[:-len(ext)]
            break
    return f"{base}_report"


def _get_display_cutsites(ref_seq: str):
    """获取用于热图标注的cutsite信息"""
    try:
        from crisviper import get_amplicon_structure
        cs = get_amplicon_structure(ref_seq)
        if cs:
            log.info("自动检测到 %d 个cutsite区域（用于热图标注）", len(cs))
        return cs
    except Exception:
        return None


if __name__ == "__main__":
    try:
        mp.set_start_method('fork')
    except RuntimeError:
        pass
    main()
