import sys, os

BASE = "/Users/mac/Documents/软件开发/谱系示踪分析/lineage-tracer"
sys.path.insert(0, BASE)

from ltlib import (
    read_reference_fasta, read_queries_tsv,
    Pipeline, PipelineConfig, QueryRecord,
    save_alignment_results, generate_report, save_text_report,
    fastq_to_dataframe, get_amplicon_structure,
)

PRINT = lambda *a: print(*a, flush=True)


def main():
    DATA = os.path.join(BASE, "example_data")
    OUT = os.path.join(BASE, "results")
    os.makedirs(OUT, exist_ok=True)

    ref_seq = read_reference_fasta(os.path.join(DATA, "reference.fa"))
    PRINT(f"参考序列: {len(ref_seq)}bp")

    # ── TSV 全量 ──
    all_q = read_queries_tsv(os.path.join(DATA, "test_queries.tsv"))
    PRINT(f"\nTSV: {len(all_q)} 条")
    recs = [QueryRecord(**q) for q in all_q]

    config = PipelineConfig(
        primer5_len=23, primer3_len=33,
        primer5_threshold=19, primer3_threshold=29,
        min_reads_snv=10, min_reads_indel=3, threads=4,
        lineage_mode=True,
        auto_detect_cutsites=True,
    )
    p = Pipeline(config=config, ref_seq=ref_seq)
    r = p.run(recs)

    PRINT(f"\n>>> TSV Pipeline:")
    PRINT(f"  总数:{r.stats.total_queries} 成功:{r.stats.successful} 失败:{r.stats.failed}")
    PRINT(f"  突变:{r.stats.mutated_sequences} 未突变:{r.stats.unmutated_sequences}")

    # WT引物验证（新方案：全长输出使用WT引物替换query引物）
    wt_p5 = sum(1 for rr in r.results if rr.success and rr.aligned_query[:23] == rr.aligned_ref[:23])
    wt_p3 = sum(1 for rr in r.results if rr.success and rr.aligned_query[-33:] == rr.aligned_ref[-33:])
    PRINT(f"  WT引物验证: p5={wt_p5}/{r.stats.successful}, p3={wt_p3}/{r.stats.successful}")

    results_dict = [rr.to_dict() for rr in r.results]
    save_alignment_results(results_dict, os.path.join(OUT, "results.json"), "json")

    # ── FASTQ 全量 ──
    PRINT(f"\nFASTQ: converting...")
    rows = fastq_to_dataframe(os.path.join(DATA, "test.fastq.gz"), "ESPC2")
    PRINT(f"  {len(rows)} 条")
    r2 = p.run([QueryRecord(**r) for r in rows])
    PRINT(f">>> FASTQ Pipeline:")
    PRINT(f"  总数:{r2.stats.total_queries} 成功:{r2.stats.successful} 失败:{r2.stats.failed}")

    # ── 报告生成（全量数据）──
    PRINT(f"\n生成报告（{len(results_dict)} 条全量数据）...")
    cutsites = get_amplicon_structure(ref_seq)
    PRINT(f"  检测到 {len(cutsites)} 个Target cutsite区域")
    generate_report(results_dict, os.path.join(OUT, "report"), "html",
                    ref_seq=ref_seq, ref_length=len(ref_seq),
                    cutsites=cutsites)
    save_text_report(results_dict, OUT, ref_seq=ref_seq, cutsites=cutsites)

    PRINT(f"\n{'='*50}")
    PRINT(f"输出文件:")
    for fn in sorted(os.listdir(OUT)):
        fp = os.path.join(OUT, fn)
        PRINT(f"  {fn}: {os.path.getsize(fp):,} bytes")
    PRINT(f"{'='*50}")
    PRINT(f"全流程测试完成 ✓")


if __name__ == "__main__":
    main()
