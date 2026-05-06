#!/usr/bin/env python3
"""批量比对运行器（谱系模式）"""
import sys, json, time
sys.path.insert(0, '.')
import multiprocessing as mp
try:
    mp.set_start_method('fork')
except RuntimeError:
    pass
from concurrent.futures import ProcessPoolExecutor, as_completed
from carlin_tool import _process_align_chunk, generate_report
from affine_gap_alignment import get_amplicon_structure


def main():
    ref = "TATGTGTGGGAGGGCTAAGAGGCCGCCGGACTGCACGACAGTCGACGATGGAGTCGACACGACTCGCGCATACGATGGAGTCGACTACAGTCGCTACGACGATGGAGTCGCGAGCGCTATGAGCGACTATGGAGTCGATACGATACGCGCACGCTATGGAGTCGAGAGCGCGCTCGTCAACGATGGAGTCGCGACTGTACGCACTCGCGATGGAGTCGATAGTATGCGTACACGCGATGGAGTCGACTGCACGACAGTCGACTATGGAGTCGATACGTAGCACGCACATGATGGGAGCTAGCTGTGCCTTCTAGTTGCCAGCCATCTGTTGT"

    # 自动推断cutsite位置
    cutsites = get_amplicon_structure(ref)
    if not cutsites:
        print("错误: 无法推断cutsite位置", file=sys.stderr)
        sys.exit(1)
    print(f"检测到 {len(cutsites)} 个cutsite区域", flush=True)

    queries = []
    with open('example_data/test_queries.tsv') as f:
        header = f.readline().strip().split('\t')
        cm = {n: i for i, n in enumerate(header)}
        for line in f:
            p = line.strip().split('\t')
            queries.append({
                'readName': p[cm['readName']], 'cellBC': p[cm['cellBC']],
                'UMI': p[cm['UMI']], 'readCount': int(p[cm['readCount']]),
                'seq': p[cm['seq']]
            })

    print(f"总序列: {len(queries)}", flush=True)

    CHUNK_SIZE = 300
    chunks = [queries[i:i+CHUNK_SIZE] for i in range(0, len(queries), CHUNK_SIZE)]
    print(f"共 {len(chunks)} 批", flush=True)

    t0 = time.time()
    results = []
    with ProcessPoolExecutor(max_workers=12) as ex:
        futs = {ex.submit(_process_align_chunk, ref, ch,
                          2.0, -3.0, -2.0, -0.1,
                          True, True, cutsites,
                          1.0, 2.0, 2.0, 3, 0.34, 3): ch
                for ch in chunks}
        for f in as_completed(futs):
            res = f.result()
            results.extend(res)
            pct = len(results) * 100 / len(queries)
            elapsed = time.time() - t0
            rate = len(results) / elapsed if elapsed > 0 else 0
            print(f"  进度: {len(results)}/{len(queries)} ({pct:.1f}%), {rate:.0f} seq/s", flush=True)

    # 过滤掉无法锚定的序列
    valid = [r for r in results if "error" not in r]
    discarded = len(results) - len(valid)
    if discarded:
        print(f"丢弃 {discarded} 条无法锚定的序列", flush=True)

    with open('results/correct.json', 'w') as f:
        json.dump(valid, f, indent=2)

    elapsed = time.time() - t0
    print(f"\n完成! {len(valid)} 条有效结果 -> results/correct.json ({elapsed:.1f}s)", flush=True)

    # ── 生成HTML报告 ──
    print("\n生成HTML报告...", flush=True)
    generate_report(
        results=valid,
        output_path='results/report.html',
        format='html',
        ref_length=len(ref),
        ref_seq=ref,
        cutsites=cutsites,
        allele_top_n=50
    )
    print("报告已保存至 results/report.html", flush=True)


if __name__ == '__main__':
    main()
