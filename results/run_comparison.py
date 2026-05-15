#!/usr/bin/env python3
"""对比测试：矫正开/关 两组比对，输出比对结果 + HTML报告"""
import sys, os, time, json, shutil

# OMP threading protection (must be before any import)
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

import multiprocessing as mp
mp.set_start_method('fork', force=True)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ltlib.models import PipelineConfig, QueryRecord
from ltlib.pipeline import Pipeline
from ltlib.lineage import get_amplicon_structure
from ltlib import read_reference_fasta, read_queries_tsv
from ltlib.reporting import save_alignment_results, generate_report
from ltlib.logging_config import get_logger, setup_logging

setup_logging(verbose=True)
log = get_logger(__name__)

# ── 参数 ──
TOP_N = None          # None = 全部序列
THREADS = 12
CHUNK_SIZE = 100       # 全量测试用稍大batch提升吞吐

ref_path = "example_data/reference.fa"
queries_path = "example_data/test_queries.tsv"
out_root = "results"

# ── 读取数据 ──
log.info("读取reference...")
ref_seq = read_reference_fasta(ref_path)
log.info("Ref length: %d bp", len(ref_seq))

log.info("读取queries...")
queries = read_queries_tsv(queries_path)
log.info("Total queries: %d", len(queries))

queries.sort(key=lambda q: q['readCount'], reverse=True)
records = [QueryRecord(**q) for q in queries] if TOP_N is None else [QueryRecord(**q) for q in queries[:TOP_N]]
log.info("总序列: %d 条 (readCount range: %d ~ %d)", len(records),
         records[-1].readCount, records[0].readCount)

cutsites = get_amplicon_structure(ref_seq)
log.info("检测到 %d 个cutsite区域", len(cutsites))

# ══════════════════════════════════════════════════════════════════
# 通用 DP 特征参数（两轮共享）
# ══════════════════════════════════════════════════════════════════
DP_FEATURES = dict(
    match_score=2.0, mismatch_penalty=-3.0, gap_open=-2.0, gap_extend=-0.1,
    lineage_mode=True,
    gap_exit_bonus=-1.0,
    short_match_window=3, short_match_discount=0.5,
    dense_mismatch_window=6, dense_mismatch_penalty=-2.0,
    homology_window=8, homology_penalty=-1.0,
    isolated_base_penalty=-2.0,
    cutsite_gap_scale=1.0, flank_gap_scale=2.0, far_gap_scale=6.0,
    flank_width=3, mismatch_density_threshold=0.34, mutation_window=3,
)
COMMON_PARAMS = dict(
    threads=THREADS,
    chunk_size=CHUNK_SIZE,
    auto_detect_cutsites=True,
    # Allele过滤
    min_reads_snv=10,
    min_reads_indel=3,
    # 报告相关
    call_alleles_enabled=False,
)


# ══════════════════════════════════════════════════════════════════
def run_one(label: str, corrections_config: dict, subdir: str):
    """运行一轮比对，保存结果 + HTML报告"""
    log.info("=" * 60)
    log.info("  运行: %s", label)
    log.info("=" * 60)

    config = PipelineConfig(**DP_FEATURES, **COMMON_PARAMS, **corrections_config)

    os.makedirs(os.path.join(out_root, subdir), exist_ok=True)
    prefix = os.path.join(out_root, subdir, "result")

    pipeline = Pipeline(config=config, ref_seq=ref_seq)
    pipeline.cutsites = cutsites

    log.info("开始比对...")
    t0 = time.time()
    pipeline_result = pipeline.run(records)
    elapsed = time.time() - t0

    results = [r.to_dict() for r in pipeline_result.results]
    successful = [r for r in pipeline_result.results if r.success]

    # 保存比对结果 (json + tsv)
    save_alignment_results(results, prefix, "all")

    # 生成HTML报告
    report_path = os.path.join(out_root, subdir, "report")
    generate_report(
        results, report_path, "html",
        ref_length=len(ref_seq), ref_seq=ref_seq,
        cutsites=cutsites,
        allele_window_start=0, allele_window_end=None,
        allele_top_n=50, version="3.0.0",
    )

    # 也保存JSON报告
    generate_report(
        results, report_path, "json",
        ref_length=len(ref_seq), ref_seq=ref_seq,
        cutsites=cutsites,
        allele_window_start=0, allele_window_end=None,
        allele_top_n=50, version="3.0.0",
    )

    # 统计
    mutated = [r for r in successful if r.stats and r.stats.has_mutation]
    unmutated = [r for r in successful if r.stats and not r.stats.has_mutation]

    summary = {
        "label": label,
        "total": len(results),
        "successful": len(successful),
        "failed": len(results) - len(successful),
        "mutated": len(mutated),
        "unmutated": len(unmutated),
        "mutated_reads": sum(r.query.readCount for r in mutated),
        "total_reads": sum(r.query.readCount for r in successful),
        "time_seconds": round(elapsed, 1),
    }

    with open(os.path.join(out_root, subdir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    log.info("[%s] 完成: %d/%d 成功, 突变率 %.1f%%, 耗时 %.1fs",
             label, len(successful), len(results),
             len(mutated)/len(successful)*100 if successful else 0,
             elapsed)
    return summary


# ══════════════════════════════════════════════════════════════════
# 第一轮：矫正关闭
# ══════════════════════════════════════════════════════════════════
summaries = []
s1 = run_one(
    "矫正关闭 (corrections=off)",
    dict(
        repeat_correction_mode="off",
        enable_target_misalignment_correction=False,
        enable_isolated_match_removal=False,
        enable_dense_mismatch_correction=False,
        enable_point_mutation_filtering=False,
    ),
    "uncorrected",
)
summaries.append(s1)

# ══════════════════════════════════════════════════════════════════
# 第二轮：矫正开启（auto模式）
# ══════════════════════════════════════════════════════════════════
s2 = run_one(
    "矫正开启 (corrections=auto)",
    dict(
        repeat_correction_mode="auto",
        enable_target_misalignment_correction=True,
        enable_isolated_match_removal=True,
        enable_dense_mismatch_correction=True,
        enable_point_mutation_filtering=True,
    ),
    "corrected",
)
summaries.append(s2)

# ══════════════════════════════════════════════════════════════════
# 对比摘要
# ══════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("  对比测试结果摘要")
print("=" * 70)

rows = []
for s in summaries:
    rows.append(f"  {s['label']}")
    rows.append(f"    成功: {s['successful']}/{s['total']}  |  失败: {s['failed']}")
    rows.append(f"    突变: {s['mutated']} ({s['mutated']/s['successful']*100:.1f}%)  |  野生型: {s['unmutated']}")
    if s['total_reads'] > 0:
        rows.append(f"    突变reads: {s['mutated_reads']}/{s['total_reads']} ({s['mutated_reads']/s['total_reads']*100:.1f}%)")
    rows.append(f"    耗时: {s['time_seconds']}s")
    rows.append(f"    输出: results/{'uncorrected' if '关闭' in s['label'] else 'corrected'}/")
    rows.append("")

print("\n".join(rows))
print("=" * 70)

# 写入总体对比
comparison = {
    "modes": [
        {"label": "uncorrected", "path": "results/uncorrected/"},
        {"label": "corrected", "path": "results/corrected/"},
    ],
    "uncorrected": s1,
    "corrected": s2,
}
with open(os.path.join(out_root, "comparison.json"), "w") as f:
    json.dump(comparison, f, indent=2, ensure_ascii=False)

print(f"总体对比: results/comparison.json")
