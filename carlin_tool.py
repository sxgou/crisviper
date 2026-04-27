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

  # 将FASTQ转换为FASTA
  carlin_tool.py convert fastq-to-fasta --fastq reads.fastq.gz --output reads.fasta

  # 并行批量比对TSV文件（自动使用所有CPU核心）
  carlin_tool.py align --reference ref.fasta --queries queries.tsv --output alignments.json

  # 同时输出JSON和TSV，并生成分析报告
  carlin_tool.py align --reference ref.fasta --queries queries.tsv --output results/prefix --format all --report json

  # 生成HTML格式分析报告
  carlin_tool.py align --reference ref.fasta --queries queries.tsv --output results/prefix --format json --report html

  # 指定8个并行进程
  carlin_tool.py align --reference ref.fasta --queries queries.tsv --output alignments.json --threads 8

  # 并行批量比对FASTA文件
  carlin_tool.py align --reference ref.fasta --queries queries.fasta --output alignments.json --threads 4
"""

import argparse
import sys
import os
import gzip
import json
import csv
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple, Optional
from pathlib import Path

# 导入比对算法
try:
    from affine_gap_alignment import affine_gap_alignment, calculate_alignment_stats
except ImportError:
    print("错误: 找不到affine_gap_alignment.py，请确保它在当前目录或PYTHONPATH中", file=sys.stderr)
    sys.exit(1)

# 尝试导入BioPython，如果不可用则给出友好错误
try:
    from Bio import SeqIO
except ImportError:
    print("错误: 需要BioPython库。请安装: pip install biopython", file=sys.stderr)
    sys.exit(1)

__version__ = "2.1.0"


def fastq_to_dataframe(fastq_path: str, sample_name: str = "sample") -> List[Dict]:
    """
    将FASTQ文件转换为字典列表，每个字典代表一条唯一序列

    参数:
        fastq_path: FASTQ文件路径（支持.gz）
        sample_name: 样本名称

    返回:
        字典列表，每个字典包含readName, cellBC, UMI, readCount, seq字段
    """
    counts = {}
    if fastq_path.endswith('.gz'):
        with gzip.open(fastq_path, 'rt') as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = str(record.seq)
                counts[seq] = counts.get(seq, 0) + 1
    else:
        for record in SeqIO.parse(fastq_path, "fastq"):
            seq = str(record.seq)
            counts[seq] = counts.get(seq, 0) + 1

    # 构建字典列表
    rows = []
    for i, (seq, count) in enumerate(counts.items()):
        rows.append({
            "readName": f"{sample_name}_seq{i+1}",
            "cellBC": sample_name,          # 虚拟细胞条形码
            "UMI": f"UMI{i+1}",             # 虚拟UMI
            "readCount": count,             # 该序列的观测次数
            "seq": seq
        })

    print(f"从 {fastq_path} 中读取了 {len(rows)} 条唯一序列（总计 {sum(counts.values())} 条reads）。")
    return rows


def fastq_to_fasta(fastq_path: str, output_fasta: str, sample_name: str = "sample") -> None:
    """
    将FASTQ文件转换为FASTA格式，头部包含元数据

    参数:
        fastq_path: 输入FASTQ文件（支持.gz）
        output_fasta: 输出FASTA文件路径
        sample_name: 样本名称
    """
    counts = {}
    if fastq_path.endswith('.gz'):
        with gzip.open(fastq_path, 'rt') as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = str(record.seq)
                counts[seq] = counts.get(seq, 0) + 1
    else:
        for record in SeqIO.parse(fastq_path, "fastq"):
            seq = str(record.seq)
            counts[seq] = counts.get(seq, 0) + 1

    with open(output_fasta, 'w') as f:
        for i, (seq, count) in enumerate(counts.items()):
            read_name = f"{sample_name}_seq{i+1}"
            f.write(f">{read_name} cellBC={sample_name} UMI=UMI{i+1} readCount={count}\n{seq}\n")

    total_reads = sum(counts.values())
    print(f"已将 {len(counts)} 条唯一序列（{total_reads} 条reads）写入 {output_fasta}")


def save_tsv(rows: List[Dict], output_path: str) -> None:
    """将字典列表保存为TSV文件"""
    fieldnames = ["readName", "cellBC", "UMI", "readCount", "seq"]
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    print(f"数据已保存至 {output_path}")


def read_reference_fasta(fasta_path: str) -> str:
    """读取reference FASTA文件，返回序列字符串"""
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            return str(record.seq).upper()
    except Exception as e:
        print(f"错误: 读取reference文件失败 - {e}", file=sys.stderr)
        sys.exit(1)
    print(f"错误: reference文件为空: {fasta_path}", file=sys.stderr)
    sys.exit(1)


def read_queries_tsv(tsv_path: str) -> List[Dict]:
    """
    读取查询序列TSV文件

    期望格式: readName, cellBC, UMI, readCount, seq
    """
    rows = []
    try:
        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for line_num, row in enumerate(reader, start=2):
                if 'seq' not in row or not row['seq'].strip():
                    print(f"警告: 第{line_num}行缺少序列字段，跳过", file=sys.stderr)
                    continue
                row['seq'] = row['seq'].upper()
                row['readCount'] = int(row.get('readCount', 1))
                rows.append(row)
    except Exception as e:
        print(f"错误: 读取TSV文件失败 - {e}", file=sys.stderr)
        sys.exit(1)

    if not rows:
        print(f"错误: TSV文件为空或不包含有效序列: {tsv_path}", file=sys.stderr)
        sys.exit(1)

    print(f"从 {tsv_path} 中读取了 {len(rows)} 条查询序列。")
    return rows


def read_queries_fasta(fasta_path: str) -> List[Dict]:
    """
    读取FASTA格式查询序列

    头部可包含元数据: >readName cellBC=sample1 UMI=UMI1 readCount=2
    """
    rows = []
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            header = record.description
            seq = str(record.seq).upper()

            read_name = record.id
            cellBC = "unknown"
            UMI = "unknown"
            readCount = 1

            if "cellBC=" in header:
                for part in header.split():
                    if part.startswith("cellBC="):
                        cellBC = part.split("=")[1]
                    elif part.startswith("UMI="):
                        UMI = part.split("=")[1]
                    elif part.startswith("readCount="):
                        try:
                            readCount = int(part.split("=")[1])
                        except ValueError:
                            pass

            rows.append({
                "readName": read_name,
                "cellBC": cellBC,
                "UMI": UMI,
                "readCount": readCount,
                "seq": seq
            })
    except Exception as e:
        print(f"错误: 读取FASTA文件失败 - {e}", file=sys.stderr)
        sys.exit(1)

    if not rows:
        print(f"错误: FASTA文件为空或格式不正确: {fasta_path}", file=sys.stderr)
        sys.exit(1)

    print(f"从 {fasta_path} 中读取了 {len(rows)} 条查询序列。")
    return rows


def _align_single(ref_seq: str, query: Dict,
                  match_score: float, mismatch_penalty: float,
                  gap_open: float, gap_extend: float,
                  semi_global: bool) -> Dict:
    """
    比对单条序列（供并行调用）

    参数:
        ref_seq: reference序列
        query: 单条查询序列字典
        比对参数

    返回:
        比对结果字典
    """
    query_seq = query["seq"]
    try:
        score, aligned_ref, aligned_query, stats = affine_gap_alignment(
            ref_seq, query_seq,
            match_score=match_score,
            mismatch_penalty=mismatch_penalty,
            gap_open=gap_open,
            gap_extend=gap_extend,
            semi_global=semi_global
        )

        return {
            "readName": query["readName"],
            "cellBC": query.get("cellBC", "unknown"),
            "UMI": query.get("UMI", "unknown"),
            "readCount": query.get("readCount", 1),
            "score": score,
            "aligned_ref": aligned_ref,
            "aligned_query": aligned_query,
            "stats": stats
        }

    except Exception as e:
        return {
            "readName": query["readName"],
            "cellBC": query.get("cellBC", "unknown"),
            "UMI": query.get("UMI", "unknown"),
            "readCount": query.get("readCount", 1),
            "error": str(e),
            "score": None,
            "aligned_ref": None,
            "aligned_query": None,
            "stats": None
        }


def parallel_align(reference_seq: str, queries: List[Dict],
                   match_score: float = 2.0, mismatch_penalty: float = -3.0,
                   gap_open: float = -2.0, gap_extend: float = -0.1,
                   semi_global: bool = True,
                   threads: int = None) -> List[Dict]:
    """
    并行批量比对查询序列到reference

    使用多进程并行执行序列比对，充分利用多核CPU加速。

    参数:
        reference_seq: reference序列
        queries: 查询序列字典列表
        比对参数
        threads: 并行进程数（默认：CPU核心数）

    返回:
        比对结果字典列表
    """
    total = len(queries)
    if total == 0:
        return []

    if threads is None or threads < 1:
        threads = multiprocessing.cpu_count()

    # 限制线程数不超过序列数
    threads = min(threads, total)

    print(f"  使用 {threads} 个并行进程处理 {total} 条序列...")

    # 准备参数元组供并行调用
    from functools import partial
    align_func = partial(
        _align_single,
        reference_seq,
        match_score=match_score,
        mismatch_penalty=mismatch_penalty,
        gap_open=gap_open,
        gap_extend=gap_extend,
        semi_global=semi_global
    )

    results = []
    completed = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=threads) as executor:
        future_to_query = {executor.submit(align_func, query): query for query in queries}

        for future in as_completed(future_to_query):
            result = future.result()
            results.append(result)

            completed += 1
            if "error" in result:
                failed += 1

            if completed % 100 == 0 or completed == total:
                pct = completed / total * 100
                print(f"\r  进度: {completed}/{total} ({pct:.1f}%)", end="", file=sys.stderr)

    print(file=sys.stderr)  # 换行
    print(f"批量比对完成: {completed - failed} 成功, {failed} 失败")
    return results


def save_alignment_results(results: List[Dict], output_path: str, format: str = "json") -> None:
    """
    保存比对结果

    参数:
        results: 比对结果列表
        output_path: 输出文件路径
        format: 输出格式 (json, tsv, all)
    """
    if format == "json":
        path = _ensure_extension(output_path, ".json")
        with open(path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"JSON结果已保存至 {path}")

    elif format == "tsv":
        path = _ensure_extension(output_path, ".tsv")
        _save_tsv_results(results, path)

    elif format == "all":
        # 同时输出JSON和TSV
        json_path = _ensure_extension(output_path, ".json")
        tsv_path = _ensure_extension(output_path, ".tsv")

        with open(json_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"JSON结果已保存至 {json_path}")

        _save_tsv_results(results, tsv_path)

    else:
        print(f"错误: 不支持的输出格式: {format}", file=sys.stderr)
        sys.exit(1)


def _ensure_extension(path: str, ext: str) -> str:
    """确保文件路径有指定的扩展名"""
    if not path.endswith(ext):
        return path + ext
    return path


def _save_tsv_results(results: List[Dict], output_path: str) -> None:
    """将比对结果保存为TSV格式，包含对齐后的序列"""
    simplified = []
    for result in results:
        if "error" in result:
            simplified.append({
                "readName": result["readName"],
                "cellBC": result["cellBC"],
                "UMI": result["UMI"],
                "readCount": result["readCount"],
                "score": "NA",
                "matches": "NA",
                "mismatches": "NA",
                "gaps_in_query": "NA",
                "similarity": "NA",
                "aligned_ref": "NA",
                "aligned_query": "NA",
                "error": result["error"]
            })
        else:
            simplified.append({
                "readName": result["readName"],
                "cellBC": result["cellBC"],
                "UMI": result["UMI"],
                "readCount": result["readCount"],
                "score": result["score"],
                "matches": result["stats"]["matches"],
                "mismatches": result["stats"]["mismatches"],
                "gaps_in_query": result["stats"]["gaps_in_query"],
                "similarity": result["stats"]["similarity"],
                "aligned_ref": result.get("aligned_ref", ""),
                "aligned_query": result.get("aligned_query", ""),
                "error": ""
            })

    fieldnames = ["readName", "cellBC", "UMI", "readCount", "score",
                  "matches", "mismatches", "gaps_in_query", "similarity",
                  "aligned_ref", "aligned_query", "error"]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(simplified)

    print(f"TSV结果已保存至 {output_path}")


# ═══════════════════════════════════════════════════════════════
# 分析报告生成
# ═══════════════════════════════════════════════════════════════

def generate_report(results: List[Dict], output_path: str, format: str = "json") -> None:
    """
    生成突变分析报告

    参数:
        results: 比对结果列表
        output_path: 输出文件路径
        format: 报告格式 (json, html)
    """
    total_sequences = len(results)
    successful = [r for r in results if "error" not in r]
    failed = [r for r in results if "error" in r]

    total_successful = len(successful)
    total_failed = len(failed)
    total_reads = sum(r.get("readCount", 1) for r in successful)
    total_reads_all = sum(r.get("readCount", 1) for r in results)

    # 统计突变
    mutated_seqs = []
    unmutated_seqs = []
    for r in successful:
        stats = r["stats"]
        has_mismatch = stats["mismatches"] > 0
        has_deletion = stats["gaps_in_query"] > 0
        has_insertion = stats["gaps_in_ref"] > 0
        if has_mismatch or has_deletion or has_insertion:
            mutated_seqs.append(r)
        else:
            unmutated_seqs.append(r)

    total_mutated = len(mutated_seqs)
    total_unmutated = len(unmutated_seqs)
    mutated_reads = sum(r.get("readCount", 1) for r in mutated_seqs)
    efficiency = total_mutated / total_successful * 100 if total_successful > 0 else 0.0

    # 突变类型统计（序列数 + Reads数）
    only_insertion = {"sequences": 0, "reads": 0}
    only_deletion = {"sequences": 0, "reads": 0}
    only_substitution = {"sequences": 0, "reads": 0}
    insertion_and_deletion = {"sequences": 0, "reads": 0}
    insertion_and_substitution = {"sequences": 0, "reads": 0}
    deletion_and_substitution = {"sequences": 0, "reads": 0}
    all_three = {"sequences": 0, "reads": 0}

    ins_lengths = []     # 插入长度列表
    del_lengths = []     # 删除长度列表
    total_mismatches = 0

    for r in mutated_seqs:
        rc = r.get("readCount", 1)
        stats = r["stats"]
        has_ins = stats["gaps_in_ref"] > 0
        has_del = stats["gaps_in_query"] > 0
        has_sub = stats["mismatches"] > 0

        if has_ins and not has_del and not has_sub:
            only_insertion["sequences"] += 1
            only_insertion["reads"] += rc
        elif has_del and not has_ins and not has_sub:
            only_deletion["sequences"] += 1
            only_deletion["reads"] += rc
        elif has_sub and not has_ins and not has_del:
            only_substitution["sequences"] += 1
            only_substitution["reads"] += rc
        elif has_ins and has_del and not has_sub:
            insertion_and_deletion["sequences"] += 1
            insertion_and_deletion["reads"] += rc
        elif has_ins and has_sub and not has_del:
            insertion_and_substitution["sequences"] += 1
            insertion_and_substitution["reads"] += rc
        elif has_del and has_sub and not has_ins:
            deletion_and_substitution["sequences"] += 1
            deletion_and_substitution["reads"] += rc
        elif has_ins and has_del and has_sub:
            all_three["sequences"] += 1
            all_three["reads"] += rc

        # 收集插入长度
        for block in stats.get("gap_blocks_ref", []):
            if isinstance(block, (list, tuple)) and len(block) >= 2:
                ins_lengths.append(block[1] - block[0])
            elif isinstance(block, int):
                ins_lengths.append(block)

        # 收集删除长度
        for block in stats.get("gap_blocks_query", []):
            if isinstance(block, (list, tuple)) and len(block) >= 2:
                del_lengths.append(block[1] - block[0])
            elif isinstance(block, int):
                del_lengths.append(block)

        total_mismatches += stats["mismatches"]

    avg_insertion_len = sum(ins_lengths) / len(ins_lengths) if ins_lengths else 0.0
    avg_deletion_len = sum(del_lengths) / len(del_lengths) if del_lengths else 0.0
    max_insertion_len = max(ins_lengths) if ins_lengths else 0
    max_deletion_len = max(del_lengths) if del_lengths else 0

    # 构建报告字典
    report = {
        "tool": "CARLIN序列分析工具",
        "version": __version__,
        "summary": {
            "total_sequences": total_sequences,
            "total_reads_all": total_reads_all,
            "successful_alignments": total_successful,
            "failed_alignments": total_failed,
            "total_reads_successful": total_reads,
            "mutated_sequences": total_mutated,
            "unmutated_sequences": total_unmutated,
            "mutated_reads": mutated_reads,
            "editing_efficiency_pct": round(efficiency, 2)
        },
        "mutation_types": {
            "only_insertion": only_insertion,
            "only_deletion": only_deletion,
            "only_substitution": only_substitution,
            "insertion_and_deletion": insertion_and_deletion,
            "insertion_and_substitution": insertion_and_substitution,
            "deletion_and_substitution": deletion_and_substitution,
            "insertion_deletion_substitution": all_three
        },
        "mutation_stats": {
            "total_point_mutations": total_mismatches,
            "total_insertion_events": len(ins_lengths),
            "total_deletion_events": len(del_lengths),
            "avg_insertion_length": round(avg_insertion_len, 2),
            "avg_deletion_length": round(avg_deletion_len, 2),
            "max_insertion_length": max_insertion_len,
            "max_deletion_length": max_deletion_len
        },
        "mutated_sequences_detail": [
            {
                "readName": r["readName"],
                "readCount": r.get("readCount", 1),
                "mismatches": r["stats"]["mismatches"],
                "gaps_in_ref": r["stats"]["gaps_in_ref"],
                "gaps_in_query": r["stats"]["gaps_in_query"],
                "similarity": r["stats"]["similarity"]
            }
            for r in mutated_seqs
        ]
    }

    if format == "json":
        path = _ensure_extension(output_path, ".json")
        with open(path, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        print(f"JSON分析报告已保存至 {path}")

    elif format == "html":
        path = _ensure_extension(output_path, ".html")
        _save_report_html(report, path)
        print(f"HTML分析报告已保存至 {path}")

    else:
        print(f"错误: 不支持的报告格式: {format}", file=sys.stderr)
        sys.exit(1)


def _save_report_html(report: dict, output_path: str) -> None:
    """将报告保存为自包含的HTML文件"""
    s = report["summary"]
    mt = report["mutation_types"]
    ms = report["mutation_stats"]
    detail = report["mutated_sequences_detail"]

    # 生成突变类型表格行（序列数 + Reads数）
    type_rows = ""
    type_labels = [
        ("only_substitution", "仅点突变（替换）"),
        ("only_deletion", "仅删除"),
        ("only_insertion", "仅插入"),
        ("insertion_and_deletion", "插入 + 删除"),
        ("insertion_and_substitution", "插入 + 点突变"),
        ("deletion_and_substitution", "删除 + 点突变"),
        ("insertion_deletion_substitution", "插入 + 删除 + 点突变"),
    ]
    for key, label in type_labels:
        val = mt.get(key, {"sequences": 0, "reads": 0})
        if isinstance(val, dict):
            seq_count = val.get("sequences", 0)
            reads_count = val.get("reads", 0)
        else:
            # 兼容旧格式（纯数字）
            seq_count = val
            reads_count = val
        seq_pct = seq_count / s["mutated_sequences"] * 100 if s["mutated_sequences"] > 0 else 0
        reads_pct = reads_count / s["mutated_reads"] * 100 if s["mutated_reads"] > 0 else 0
        type_rows += (
            f"<tr>"
            f"<td>{label}</td>"
            f"<td>{seq_count}</td><td>{seq_pct:.1f}%</td>"
            f"<td>{reads_count}</td><td>{reads_pct:.1f}%</td>"
            f"</tr>\n"
        )

    # 生成详细序列表格行
    detail_rows = ""
    for i, d in enumerate(detail[:100]):  # 最多显示100条
        detail_rows += (
            f"<tr>"
            f"<td>{i+1}</td>"
            f"<td>{d['readName']}</td>"
            f"<td>{d['readCount']}</td>"
            f"<td>{d['mismatches']}</td>"
            f"<td>{d['gaps_in_ref']}</td>"
            f"<td>{d['gaps_in_query']}</td>"
            f"<td>{d['similarity']}</td>"
            f"</tr>\n"
        )
    if len(detail) > 100:
        detail_rows += f"<tr><td colspan='7' style='text-align:center;color:#888;'>... 还有 {len(detail)-100} 条突变序列未显示</td></tr>\n"

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>CARLIN 序列分析报告</title>
<style>
  body {{ font-family: -apple-system, 'Segoe UI', Helvetica, Arial, sans-serif; margin: 0; padding: 20px; background: #f5f7fa; color: #333; }}
  .container {{ max-width: 960px; margin: 0 auto; }}
  h1 {{ color: #1a1a2e; border-bottom: 3px solid #4a6cf7; padding-bottom: 10px; }}
  h2 {{ color: #1a1a2e; margin-top: 30px; }}
  .meta {{ color: #666; font-size: 14px; margin-bottom: 20px; }}
  .summary-cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px; margin: 20px 0; }}
  .card {{ background: #fff; border-radius: 10px; padding: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); text-align: center; }}
  .card .value {{ font-size: 28px; font-weight: 700; color: #4a6cf7; }}
  .card .label {{ font-size: 13px; color: #888; margin-top: 5px; }}
  .card.highlight {{ background: #4a6cf7; color: #fff; }}
  .card.highlight .value {{ color: #fff; }}
  .card.highlight .label {{ color: rgba(255,255,255,0.85); }}
  table {{ width: 100%; border-collapse: collapse; background: #fff; border-radius: 10px; overflow: hidden; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin: 15px 0; }}
  th {{ background: #4a6cf7; color: #fff; padding: 12px 15px; text-align: left; font-weight: 600; }}
  td {{ padding: 10px 15px; border-bottom: 1px solid #eee; }}
  tr:last-child td {{ border-bottom: none; }}
  tr:hover td {{ background: #f0f4ff; }}
  .section {{ background: #fff; border-radius: 10px; padding: 20px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin: 15px 0; }}
  .stat-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 12px; }}
  .stat-item {{ padding: 8px 12px; background: #f8f9fc; border-radius: 6px; }}
  .stat-item .s-label {{ font-size: 12px; color: #888; }}
  .stat-item .s-value {{ font-size: 18px; font-weight: 600; color: #333; }}
  footer {{ text-align: center; color: #aaa; font-size: 12px; margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; }}
</style>
</head>
<body>
<div class="container">
<h1>🔬 CARLIN 序列分析报告</h1>
<div class="meta">生成时间: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')} | 工具版本: {report["version"]}</div>

<div class="summary-cards">
  <div class="card">
    <div class="value">{s["total_sequences"]}</div>
    <div class="label">总序列数</div>
  </div>
  <div class="card">
    <div class="value">{s["total_reads_all"]}</div>
    <div class="label">总Reads数</div>
  </div>
  <div class="card">
    <div class="value">{s["successful_alignments"]}</div>
    <div class="label">成功比对</div>
  </div>
  <div class="card highlight">
    <div class="value">{s["editing_efficiency_pct"]}%</div>
    <div class="label">编辑效率</div>
  </div>
</div>

<div class="summary-cards">
  <div class="card">
    <div class="value">{s["mutated_sequences"]}</div>
    <div class="label">突变序列数</div>
  </div>
  <div class="card">
    <div class="value">{s["mutated_reads"]}</div>
    <div class="label">突变Reads数</div>
  </div>
  <div class="card">
    <div class="value">{s["unmutated_sequences"]}</div>
    <div class="label">未突变序列数</div>
  </div>
  <div class="card">
    <div class="value">{s["failed_alignments"]}</div>
    <div class="label">比对失败</div>
  </div>
</div>

<h2>📊 突变类型分布</h2>
<table>
<tr><th>突变类型</th><th>序列数</th><th>序列占比</th><th>Reads数</th><th>Reads占比</th></tr>
{type_rows}
</table>

<h2>📏 突变统计指标</h2>
<div class="section">
<div class="stat-grid">
  <div class="stat-item">
    <div class="s-label">点突变总数</div>
    <div class="s-value">{ms["total_point_mutations"]}</div>
  </div>
  <div class="stat-item">
    <div class="s-label">插入事件数</div>
    <div class="s-value">{ms["total_insertion_events"]}</div>
  </div>
  <div class="stat-item">
    <div class="s-label">删除事件数</div>
    <div class="s-value">{ms["total_deletion_events"]}</div>
  </div>
  <div class="stat-item">
    <div class="s-label">平均插入长度</div>
    <div class="s-value">{ms["avg_insertion_length"]} bp</div>
  </div>
  <div class="stat-item">
    <div class="s-label">平均删除长度</div>
    <div class="s-value">{ms["avg_deletion_length"]} bp</div>
  </div>
  <div class="stat-item">
    <div class="s-label">最大插入长度</div>
    <div class="s-value">{ms["max_insertion_length"]} bp</div>
  </div>
  <div class="stat-item">
    <div class="s-label">最大删除长度</div>
    <div class="s-value">{ms["max_deletion_length"]} bp</div>
  </div>
</div>
</div>

<h2>📋 突变序列明细（前100条）</h2>
<table>
<tr><th>#</th><th>序列名</th><th>Reads</th><th>点突变</th><th>插入事件</th><th>删除事件</th><th>相似度</th></tr>
{detail_rows}
</table>

<footer>Generated by CARLIN 序列分析工具 v{report["version"]}</footer>
</div>
</body>
</html>"""
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(
        description="CARLIN序列分析命令行工具",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    subparsers = parser.add_subparsers(dest="command", required=True, help="命令")

    # ── 子命令: convert ──────────────────────────────────────────
    convert_parser = subparsers.add_parser("convert", help="转换FASTQ文件格式")
    convert_subparsers = convert_parser.add_subparsers(dest="convert_command", required=True, help="转换子命令")

    # fastq-to-tsv
    tsv_parser = convert_subparsers.add_parser("fastq-to-tsv", help="将FASTQ转换为TSV格式")
    tsv_parser.add_argument("--fastq", required=True, help="输入FASTQ文件（支持.gz压缩）")
    tsv_parser.add_argument("--output", required=True, help="输出TSV文件路径")
    tsv_parser.add_argument("--sample-name", default="sample", help="样本名称")

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
                              help="输出格式: json（默认）、tsv、all（同时输出json和tsv）")

    # 比对参数
    align_parser.add_argument("--match-score", type=float, default=2.0, help="匹配得分（默认：2.0）")
    align_parser.add_argument("--mismatch-penalty", type=float, default=-3.0, help="错配惩罚（默认：-3.0）")
    align_parser.add_argument("--gap-open", type=float, default=-2.0, help="gap开启惩罚（默认：-2.0）")
    align_parser.add_argument("--gap-extend", type=float, default=-0.1, help="gap延伸惩罚（默认：-0.1）")
    align_parser.add_argument("--global", action="store_true", dest="use_global",
                              help="使用全局比对（默认使用半全局比对）")

    # 并行参数
    align_parser.add_argument("--threads", "-t", type=int, default=None,
                              help="并行进程数（默认：自动使用所有CPU核心）")

    # 报告参数
    align_parser.add_argument("--report", choices=["json", "html"], default=None,
                              help="生成分析报告（json或html格式）")
    align_parser.add_argument("--report-output", default=None,
                              help="报告输出路径（默认：基于--output文件名）")

    args = parser.parse_args()

    # ═══════════════════════════════════════════════════════════════
    # 执行命令
    # ═══════════════════════════════════════════════════════════════
    if args.command == "convert":
        if args.convert_command == "fastq-to-tsv":
            rows = fastq_to_dataframe(args.fastq, args.sample_name)
            save_tsv(rows, args.output)

        elif args.convert_command == "fastq-to-fasta":
            fastq_to_fasta(args.fastq, args.output, args.sample_name)

    elif args.command == "align":
        # 读取reference
        print(f"读取reference序列: {args.reference}")
        reference_seq = read_reference_fasta(args.reference)
        print(f"reference序列长度: {len(reference_seq)} bp")

        # 读取查询序列
        queries_path = args.queries
        if queries_path.endswith('.tsv'):
            queries = read_queries_tsv(queries_path)
        elif queries_path.endswith(('.fasta', '.fa', '.fas')):
            queries = read_queries_fasta(queries_path)
        else:
            print(f"错误: 不支持的查询文件格式: {queries_path}", file=sys.stderr)
            print("请使用.tsv或.fasta/.fa/.fas扩展名", file=sys.stderr)
            sys.exit(1)

        # 运行并行批量比对
        if args.use_global:
            align_type = "全局"
        else:
            align_type = "半全局"

        if args.threads:
            print(f"开始并行批量比对（{align_type}，{args.threads} 线程）...")
        else:
            cpu_count = multiprocessing.cpu_count()
            print(f"开始并行批量比对（{align_type}，自动检测到 {cpu_count} 个CPU核心）...")

        results = parallel_align(
            reference_seq, queries,
            match_score=args.match_score,
            mismatch_penalty=args.mismatch_penalty,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend,
            semi_global=not args.use_global,
            threads=args.threads
        )

        # 保存结果
        save_alignment_results(results, args.output, args.format)

        # 生成分析报告
        if args.report:
            report_output = args.report_output
            if report_output is None:
                # 默认：基于--output生成报告文件名
                base = args.output
                # 去掉可能存在的扩展名
                for ext in [".json", ".tsv"]:
                    if base.endswith(ext):
                        base = base[:-len(ext)]
                        break
                report_output = f"{base}_report"

            print(f"生成{args.report.upper()}格式分析报告...")
            generate_report(results, report_output, args.report)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
