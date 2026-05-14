#!/usr/bin/env python3
"""
CARLIN序列分析命令行工具

基于仿射gap惩罚算法的序列比对工具，支持FASTQ转换和并行批量比对。

子命令:
  convert   将FASTQ文件转换为TSV或FASTA格式
  align     并行批量将TSV或FASTA格式中的序列与reference进行比对

示例:
  # 将FASTQ转换为TSV
  carlin_tool.py convert fastq-to-tsv --fastq reads.fastq.gz --output reads.tsv

  # 并行批量比对TSV文件（自动使用所有CPU核心）
  carlin_tool.py align --reference ref.fasta --queries queries.tsv --output alignments.json

  # 同时输出JSON和TSV，并生成分析报告
  carlin_tool.py align --reference ref.fasta --queries queries.tsv --output results/prefix --format all --report html

  # 指定8个并行进程
  carlin_tool.py align --reference ref.fasta --queries queries.tsv --output alignments.json --threads 8
"""

import argparse
import sys
import os
import multiprocessing as mp
from typing import List, Dict

from ltlib import (
    # 数据模型 — 类型安全
    PipelineConfig, PipelineResult, AlignmentResult, QueryRecord,

    # I/O
    fastq_to_dataframe, fastq_to_fasta, save_tsv,
    read_reference_fasta, read_queries_tsv, read_queries_fasta,

    # 管道编排
    Pipeline,

    # 报告
    save_alignment_results, generate_report,

    # 配置
    get_amplicon_structure, CutsiteRegion,

    # 日志
    get_logger, setup_logging,
)

log = get_logger(__name__)

__version__ = "3.0.0"


def _queries_to_records(queries: List[Dict]) -> List[QueryRecord]:
    """将旧版 dict 格式转为 QueryRecord"""
    return [QueryRecord(**q) for q in queries]


def main():
    setup_logging(quiet=False)
    parser = argparse.ArgumentParser(
        description="CARLIN序列分析命令行工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase verbosity (show debug output)')
    parser.add_argument('-q', '--quiet', action='store_true', help='Suppress info output, show only warnings/errors')

    subparsers = parser.add_subparsers(dest="command", required=True, help="命令")

    # ── 子命令: convert ──────────────────────────────────────────
    convert_parser = subparsers.add_parser("convert", help="转换FASTQ文件格式")
    convert_subparsers = convert_parser.add_subparsers(dest="convert_command", required=True, help="转换子命令")

    # fastq-to-tsv
    tsv_parser = convert_subparsers.add_parser("fastq-to-tsv", help="将FASTQ转换为TSV格式")
    tsv_parser.add_argument("--fastq", required=True, help="输入FASTQ文件（支持.gz压缩）")
    tsv_parser.add_argument("--output", required=True, help="输出TSV文件路径")
    tsv_parser.add_argument("--sample-name", default="sample", help="样本名称")
    tsv_parser.add_argument("--min-reads", type=int, default=1,
                            help="最小read数阈值（默认：1，即不过滤）")

    # fastq-to-fasta
    fasta_parser = convert_subparsers.add_parser("fastq-to-fasta", help="将FASTQ转换为FASTA格式")
    fasta_parser.add_argument("--fastq", required=True, help="输入FASTQ文件（支持.gz压缩）")
    fasta_parser.add_argument("--output", required=True, help="输出FASTA文件路径")
    fasta_parser.add_argument("--sample-name", default="sample", help="样本名称")

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
    align_parser.add_argument("--global", action="store_true", dest="use_global",
                              help="使用全局比对（默认使用半全局比对）")

    # 谱系示踪模式
    align_parser.add_argument("--lineage", action="store_true",
                              help="启用谱系示踪比对模式")
    align_parser.add_argument("--cutsite-scale", type=float, default=1.0)
    align_parser.add_argument("--flank-scale", type=float, default=2.0)
    align_parser.add_argument("--far-scale", type=float, default=6.0)
    align_parser.add_argument("--flank-width", type=int, default=3)
    align_parser.add_argument("--mutation-window", type=int, default=3)
    align_parser.add_argument("--density-threshold", type=float, default=0.34)
    align_parser.add_argument("--cutsites", default=None,
                              help="cutsite位置配置文件（JSON格式）")

    # 并行参数
    align_parser.add_argument("--threads", "-t", type=int, default=None)

    # 引物参数
    align_parser.add_argument("--primer5-len", type=int, default=23)
    align_parser.add_argument("--primer3-len", type=int, default=33)
    align_parser.add_argument("--primer5-threshold", type=int, default=19)
    align_parser.add_argument("--primer3-threshold", type=int, default=29)

    # Allele过滤参数
    align_parser.add_argument("--min-reads", type=int, default=1,
                              help="最小read数阈值（默认：1，即不过滤）")

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
            rows = fastq_to_dataframe(args.fastq, args.sample_name)
            if args.min_reads > 1:
                before = len(rows)
                rows = [r for r in rows if r['readCount'] >= args.min_reads]
                log.info("Allele过滤 (min_reads>=%d): %d → %d", args.min_reads, before, len(rows))
            save_tsv(rows, args.output)

        elif args.convert_command == "fastq-to-fasta":
            fastq_to_fasta(args.fastq, args.output, args.sample_name)

    elif args.command == "align":
        # 读取reference
        log.info("读取reference序列: %s", args.reference)
        ref_seq = read_reference_fasta(args.reference)
        log.info("reference序列长度: %d bp", len(ref_seq))

        # 读取查询序列
        queries_path = args.queries
        if queries_path.endswith('.tsv'):
            queries = read_queries_tsv(queries_path)
        elif queries_path.endswith(('.fasta', '.fa', '.fas')):
            queries = read_queries_fasta(queries_path)
        else:
            log.error("不支持的查询文件格式: %s", queries_path)
            sys.exit(1)

        # 转为类型安全的 QueryRecord
        query_records = _queries_to_records(queries)

        # Allele质量过滤（转换前）
        if args.min_reads > 1:
            before = len(query_records)
            query_records = [q for q in query_records if q.readCount >= args.min_reads]
            log.info("Allele过滤 (min_reads>=%d): %d → %d",
                     args.min_reads, before, len(query_records))

        # ── 构建 PipelineConfig ──
        config = PipelineConfig(
            match_score=args.match_score,
            mismatch_penalty=args.mismatch_penalty,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend,
            semi_global=not args.use_global,
            lineage_mode=args.lineage,
            cutsite_gap_scale=args.cutsite_scale,
            flank_gap_scale=args.flank_scale,
            far_gap_scale=args.far_scale,
            flank_width=args.flank_width,
            mismatch_density_threshold=args.density_threshold,
            mutation_window=args.mutation_window,
            primer5_len=args.primer5_len,
            primer3_len=args.primer3_len,
            primer5_threshold=args.primer5_threshold,
            primer3_threshold=args.primer3_threshold,
            # 智能过滤（min-reads=1时启用智能过滤，>1时用阈值过滤）
            min_reads_snv=args.min_reads if args.min_reads > 1 else 10,
            min_reads_indel=args.min_reads if args.min_reads > 1 else 3,
            threads=args.threads,
            cutsites_path=args.cutsites,
            report_format=args.report,
            allele_top_n=args.allele_top_n,
            allele_window_start=args.allele_window_start,
            allele_window_end=args.allele_window_end,
        )

        # 显示配置信息
        align_type = "全局" if args.use_global else "半全局"
        if args.lineage:
            log.info("  谱系示踪参数: cutsite倍率=%s, 侧翼倍率=%s, 远端倍率=%s",
                     args.cutsite_scale, args.flank_scale, args.far_scale)
            log.info("  突变窗口: cutsite±%s bp, mismatch密度阈值: %s",
                     args.mutation_window, args.density_threshold)

        if args.threads:
            log.info("开始并行批量比对（%s，%d 线程）...", align_type, args.threads)
        else:
            log.info("开始并行批量比对（%s，自动检测到 %d 个CPU核心）...",
                     align_type, mp.cpu_count())

        # ── 运行管道 ──
        pipeline = Pipeline(config=config, ref_seq=ref_seq)
        pipeline_result = pipeline.run(query_records)

        # ── 转换为兼容的输出格式 ──
        output_results = [r.to_dict() for r in pipeline_result.results]

        # ── 保存结果 ──
        save_alignment_results(output_results, args.output, args.format)

        # ── 生成分析报告 ──
        if args.report:
            report_output = args.report_output or _default_report_path(args.output)
            log.info("生成%s格式分析报告...", args.report.upper())
            generate_report(output_results, report_output, args.report,
                             ref_length=len(ref_seq),
                             ref_seq=ref_seq,
                             cutsites=_get_display_cutsites(ref_seq, args.lineage),
                             allele_window_start=args.allele_window_start,
                             allele_window_end=args.allele_window_end,
                             allele_top_n=args.allele_top_n,
                             version=__version__)

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


def _get_display_cutsites(ref_seq: str, lineage_mode: bool):
    """获取用于热图标注的cutsite信息"""
    try:
        from ltlib import get_amplicon_structure
        cs = get_amplicon_structure(ref_seq)
        if cs:
            log.info("自动检测到 %d 个cutsite区域（用于热图标注）", len(cs))
        return cs
    except Exception:
        return None


if __name__ == "__main__":
    # init process
    try:
        mp.set_start_method('fork')
    except RuntimeError:
        pass
    main()
