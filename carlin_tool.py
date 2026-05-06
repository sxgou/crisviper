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
import multiprocessing as mp
try:
    mp.set_start_method('fork')
except RuntimeError:
    pass
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple, Optional
from pathlib import Path

# 导入比对算法
try:
    from affine_gap_alignment import (
        affine_gap_alignment, calculate_alignment_stats,
        lineage_tracer_align, get_amplicon_structure,
        CutsiteRegion, build_gap_penalty_profile,
        correct_repetitive_misalignment, remove_isolated_matches,
        correct_target_misalignments
    )
except ImportError:
    print("错误: 找不到affine_gap_alignment.py，请确保它在当前目录或PYTHONPATH中", file=sys.stderr)
    sys.exit(1)

# 尝试导入BioPython，如果不可用则给出友好错误
try:
    from Bio import SeqIO
except ImportError:
    print("错误: 需要BioPython库。请安装: pip install biopython", file=sys.stderr)
    sys.exit(1)

# 尝试导入matplotlib用于可视化（可选）
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    _HAS_MPL = True
except ImportError:
    _HAS_MPL = False
    print("提示: 未安装matplotlib，图表功能将不可用。请安装: pip install matplotlib", file=sys.stderr)

import io
import base64
from collections import Counter

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


def _process_align_chunk(ref_seq, queries_chunk,
                          match_score, mismatch_penalty,
                          gap_open, gap_extend,
                          semi_global, lineage_mode, cutsites,
                          cutsite_gap_scale, flank_gap_scale, far_gap_scale,
                          flank_width, mismatch_density_threshold, mutation_window):
    """处理一批序列（供并行调用）"""
    chunk_results = []
    for query in queries_chunk:
        result = _align_single(
            ref_seq, query,
            match_score=match_score, mismatch_penalty=mismatch_penalty,
            gap_open=gap_open, gap_extend=gap_extend, semi_global=semi_global,
            lineage_mode=lineage_mode, cutsites=cutsites,
            cutsite_gap_scale=cutsite_gap_scale, flank_gap_scale=flank_gap_scale,
            far_gap_scale=far_gap_scale, flank_width=flank_width,
            mismatch_density_threshold=mismatch_density_threshold,
            mutation_window=mutation_window
        )
        chunk_results.append(result)
    return chunk_results


def _align_single(ref_seq: str, query: Dict,
                  match_score: float, mismatch_penalty: float,
                  gap_open: float, gap_extend: float,
                  semi_global: bool,
                  lineage_mode: bool = False,
                  cutsites: list = None,
                  cutsite_gap_scale: float = 1.0,
                  flank_gap_scale: float = 2.0,
                  far_gap_scale: float = 6.0,
                  flank_width: int = 3,
                  mismatch_density_threshold: float = 0.34,
                  mutation_window: int = 3) -> Dict:
    """
    比对单条序列（供并行调用）

    参数:
        ref_seq: reference序列
        query: 单条查询序列字典
        比对参数
        lineage_mode: 使用谱系示踪比对模式（结构感知）
        cutsites: CutsiteRegion列表（lineage_mode需要）

    返回:
        比对结果字典
    """
    query_seq = query["seq"]

    try:
        primer5_len, primer3_len = 23, 33
        q_seq = query_seq
        r_seq = ref_seq

        # ── 双端锚定: 检查Primer5和Primer3是否匹配 ──
        p5_matches = sum(1 for a, b in zip(q_seq[:primer5_len], r_seq[:primer5_len]) if a == b)
        q3_check = q_seq[-primer3_len:] if len(q_seq) >= primer3_len else ""
        p3_matches = sum(1 for a, b in zip(q3_check, r_seq[-primer3_len:]) if a == b)

        if not (len(q_seq) >= primer5_len + primer3_len and
                p5_matches >= primer5_len - 4 and
                p3_matches >= primer3_len - 4):
            return {
                "readName": query["readName"],
                "cellBC": query.get("cellBC", "unknown"),
                "UMI": query.get("UMI", "unknown"),
                "readCount": query.get("readCount", 1),
                "error": "无法双端锚定: Primer5或Primer3匹配度不足",
                "score": None,
                "aligned_ref": None,
                "aligned_query": None,
                "stats": None
            }

        # ── 双端锚定通过，进行比对 ──
        if lineage_mode and cutsites:
            # ═══════════════════════════════════════════════
            #  谱系模式：使用结构感知比对（专为CARLIN设计）
            #  - 位置依赖gap惩罚（cutsite鼓励gap，保守区抑制）
            #  - 高密度mismatch→insertion转换
            #  - 区域感知点突变过滤
            # ═══════════════════════════════════════════════
            from affine_gap_alignment import CutsiteRegion
            int_q = q_seq[primer5_len:-primer3_len]
            int_r = r_seq[primer5_len:-primer3_len]

            # 调整cutsite坐标到内部区域
            int_cutsites = []
            for cs in cutsites:
                ns = cs.start - primer5_len
                ne = cs.end - primer5_len
                if ns < 0 and ne < 0:
                    continue
                if ne >= len(int_r):
                    ne = len(int_r) - 1
                if ns < 0:
                    ns = 0
                if ns <= ne:
                    int_cutsites.append(CutsiteRegion(
                        name=cs.name, start=ns, end=ne
                    ))

            if int_q:
                # ── CGCCG前缀剥离处理 ──
                # 保证5bp prefix(CGCCG)始终锚定在内部区域起始位置
                PREFIX = 'CGCCG'
                int_q_stripped = int_q
                int_r_stripped = int_r
                prefix_matched = 0
                if int_q.startswith(PREFIX):
                    int_q_stripped = int_q[len(PREFIX):]
                    int_r_stripped = int_r[len(PREFIX):]
                    prefix_matched = len(PREFIX)
                    prefix_score = 5 * match_score

                # 调整cutsite坐标到剥离后的区域
                int_cutsites_stripped = []
                for cs in int_cutsites:
                    ns = max(0, cs.start - prefix_matched)
                    ne = max(0, cs.end - prefix_matched)
                    if ns <= ne:
                        int_cutsites_stripped.append(CutsiteRegion(
                            name=cs.name, start=ns, end=ne
                        ))

                if int_q_stripped:
                    score_int, ar_int, aq_int, stats = lineage_tracer_align(
                        int_r_stripped, int_q_stripped, int_cutsites_stripped,
                        match_score=match_score,
                        mismatch_penalty=mismatch_penalty,
                        base_gap_open=gap_open,
                        base_gap_extend=gap_extend,
                        semi_global='anchor5',
                        cutsite_gap_scale=cutsite_gap_scale,
                        flank_gap_scale=flank_gap_scale,
                        far_gap_scale=far_gap_scale,
                        flank_width=flank_width,
                        mismatch_density_threshold=mismatch_density_threshold,
                        mutation_window=mutation_window
                    )
                else:
                    score_int, ar_int, aq_int = 0.0, "", ""
                    stats = {"matches": 0, "mismatches": 0, "gaps_in_ref": 0,
                             "gaps_in_query": 0, "gap_blocks_ref": [],
                             "gap_blocks_query": [], "alignment_length": 0,
                             "similarity": 0.0, "identity": 0.0}

                if prefix_matched > 0:
                    ar_int = PREFIX + ar_int
                    aq_int = PREFIX + aq_int
                    score_int += prefix_score
                    stats['matches'] = stats.get('matches', 0) + prefix_matched
                    stats['alignment_length'] = stats.get('alignment_length', 0) + prefix_matched
                    if stats.get('alignment_length', 0) > 0:
                        stats['similarity'] = stats['matches'] / stats['alignment_length']
                        stats['identity'] = stats['matches'] / stats['alignment_length']
            else:
                score_int, ar_int, aq_int = 0.0, "", ""
                stats = {"matches": 0, "mismatches": 0, "gaps_in_ref": 0,
                         "gaps_in_query": 0, "gap_blocks_ref": [],
                         "gap_blocks_query": [], "alignment_length": 0,
                         "similarity": 0.0, "identity": 0.0}

            # 确定ar_int在int_r中的对应位置（半全局可能跳过首尾ref碱基）
            ar_bases = ar_int.replace('-', '')
            if len(ar_bases) > 0:
                n_anchor = min(20, len(ar_bases))
                anchor = ar_bases[:n_anchor]
                best_pos, best_cnt = 0, -1
                for i in range(len(int_r) - n_anchor + 1):
                    cnt = sum(1 for k in range(n_anchor) if int_r[i+k] == anchor[k])
                    if cnt > best_cnt:
                        best_cnt = cnt
                        best_pos = i
                start_in_int = best_pos
                end_in_int = start_in_int + len(ar_bases)
            else:
                start_in_int = 0
                end_in_int = 0

            # 用position tracking组装全长比对（首尾锚定）
            ar_parts = [r_seq[:primer5_len]]
            aq_parts = [q_seq[:primer5_len]]
            if start_in_int > 0:
                ar_parts.append(int_r[:start_in_int])
                aq_parts.append('-' * start_in_int)
            ri = primer5_len + start_in_int
            qi = primer5_len
            for ac, qc in zip(ar_int, aq_int):
                if ac != '-':
                    ar_parts.append(r_seq[ri])
                    ri += 1
                else:
                    ar_parts.append('-')
                if qc != '-':
                    aq_parts.append(q_seq[qi])
                    qi += 1
                else:
                    aq_parts.append('-')
            if end_in_int < len(int_r):
                remaining = len(int_r) - end_in_int
                ar_parts.append(int_r[end_in_int:])
                aq_parts.append('-' * remaining)
            ar_parts.append(r_seq[primer5_len + len(int_r):])
            aligned_ref = ''.join(ar_parts)
            aq_parts.append(q_seq[qi:])
            aligned_query = ''.join(aq_parts)

            # 矫正重复序列导致的错误比对
            aligned_ref, aligned_query, _ = correct_repetitive_misalignment(
                aligned_ref, aligned_query, r_seq
            )

            # 矫正跨Target的小片段错误比对（TAGTAT, 单碱基A）
            aligned_ref, aligned_query, _ = correct_target_misalignments(
                aligned_ref, aligned_query, r_seq
            )

            # 清除孤立匹配（合并被切断的deletion片段）
            aligned_ref, aligned_query, _ = remove_isolated_matches(
                aligned_ref, aligned_query
            )

            # 重新计分（全长）
            score = score_int + (p5_matches + p3_matches) * match_score
            stats = calculate_alignment_stats(aligned_ref, aligned_query)
            stats['score'] = score

        else:
            # ═══════════════════════════════════════════════
            #  标准模式：Gotoh仿射gap比对（内部区域）
            # ═══════════════════════════════════════════════
            int_q = q_seq[primer5_len:-primer3_len]
            int_r = r_seq[primer5_len:-primer3_len]
            if int_q:
                score2, ar_int, aq_int, _ = affine_gap_alignment(
                    int_r, int_q, match_score=match_score,
                    mismatch_penalty=mismatch_penalty, gap_open=gap_open,
                    gap_extend=gap_extend, semi_global=True
                )
            else:
                score2, ar_int, aq_int = 0, "", ""

            # 检测内部比对是否因重复元件(GAGTCGAT)导致片段化或位置偏移
            if len(aq_int) > 0:
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

                first_nongap = blocks[0][0] if blocks else 0
                last_block_end = blocks[-1][1] if blocks else 0

                shifted = first_nongap > 20
                scattered = (len(blocks) >= 2 and
                             last_block_end < len(int_r) - 5)
                needs_anchoring = shifted or scattered or len(blocks) > 2

                if needs_anchoring and len(int_q) <= len(int_r) - 10:
                    if len(blocks) >= 2:
                        front_q = ''.join(b[2] for b in blocks[:-1])
                        back_q = blocks[-1][2]
                    else:
                        back_len = min(12, len(int_q) // 3)
                        front_q = int_q[:-back_len]
                        back_q = int_q[-back_len:]

                    if len(front_q) + len(back_q) <= len(int_r):
                        gap_len = len(int_r) - len(front_q) - len(back_q)
                        new_aq = front_q + '-' * gap_len + back_q
                        f_matches = 0
                        for i in range(len(front_q)):
                            if i < len(int_r) and front_q[i] == int_r[i]:
                                f_matches += 1
                        b_start_q = len(int_q) - len(back_q)
                        b_start_r = len(int_r) - len(back_q)
                        b_matches = 0
                        for i in range(len(back_q)):
                            if (b_start_q + i < len(int_q) and
                                b_start_r + i < len(int_r) and
                                int_q[b_start_q + i] == int_r[b_start_r + i]):
                                b_matches += 1
                        total_matches = f_matches + b_matches
                        total_mismatches = len(int_q) - total_matches
                        score2 = (total_matches * match_score +
                                  total_mismatches * mismatch_penalty)
                        aq_int = new_aq
                        ar_int = int_r
                elif len(blocks) > 2:
                    fb_start, fb_end, fb_seq = blocks[0]
                    lb_start, lb_end, lb_seq = blocks[-1]
                    merged_query = ''.join(b[2] for b in blocks[:-1])
                    merged_end = fb_start + len(merged_query)
                    new_aq = (aq_int[:fb_start] + merged_query +
                              '-' * (lb_start - merged_end) + lb_seq)
                    if len(new_aq) < len(ar_int):
                        new_aq += '-' * (len(ar_int) - len(new_aq))
                    aq_int = new_aq[:len(ar_int)]

            # 补齐ar_int中缺失的ref碱基
            ar_bases = ar_int.replace('-', '')
            n_anchor = min(20, len(ar_bases))
            if n_anchor > 0:
                anchor = ar_bases[:n_anchor]
                best_pos, best_cnt = 0, -1
                for i in range(len(int_r) - n_anchor + 1):
                    cnt = sum(1 for k in range(n_anchor) if int_r[i+k] == anchor[k])
                    if cnt > best_cnt:
                        best_cnt = cnt
                        best_pos = i
                start_in_int = best_pos
                end_in_int = start_in_int + len(ar_bases)

            # 用position tracking组装全长比对
            ar_parts = [r_seq[:primer5_len]]
            aq_parts = [q_seq[:primer5_len]]
            if start_in_int > 0:
                ar_parts.append(int_r[:start_in_int])
                aq_parts.append('-' * start_in_int)
            ri = primer5_len + start_in_int
            qi = primer5_len
            for ac, qc in zip(ar_int, aq_int):
                if ac != '-':
                    ar_parts.append(r_seq[ri])
                    ri += 1
                else:
                    ar_parts.append('-')
                if qc != '-':
                    aq_parts.append(q_seq[qi])
                    qi += 1
                else:
                    aq_parts.append('-')
            if end_in_int < len(int_r):
                ar_parts.append(int_r[end_in_int:])
                aq_parts.append('-' * (len(int_r) - end_in_int))
            ar_parts.append(r_seq[primer5_len + len(int_r):])
            aligned_ref = ''.join(ar_parts)
            aq_parts.append(q_seq[qi:])
            aligned_query = ''.join(aq_parts)

            # 矫正重复序列导致的错误比对
            aligned_ref, aligned_query, _ = correct_repetitive_misalignment(
                aligned_ref, aligned_query, r_seq
            )

            # 矫正跨Target的小片段错误比对（TAGTAT, 单碱基A）
            aligned_ref, aligned_query, _ = correct_target_misalignments(
                aligned_ref, aligned_query, r_seq
            )

            # 清除孤立匹配（合并被切断的deletion片段）
            aligned_ref, aligned_query, _ = remove_isolated_matches(
                aligned_ref, aligned_query
            )

            score = score2 + (p5_matches + p3_matches) * match_score
            stats = calculate_alignment_stats(aligned_ref, aligned_query)
            stats['score'] = score

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
                   threads: int = None,
                   lineage_mode: bool = False,
                   cutsites: list = None,
                   cutsite_gap_scale: float = 1.0,
                   flank_gap_scale: float = 2.0,
                   far_gap_scale: float = 6.0,
                   flank_width: int = 3,
                   mismatch_density_threshold: float = 0.34,
                   mutation_window: int = 3) -> List[Dict]:
    """
    并行批量比对查询序列到reference

    使用多进程并行执行序列比对，充分利用多核CPU加速。

    参数:
        reference_seq: reference序列
        queries: 查询序列字典列表
        比对参数
        threads: 并行进程数（默认：CPU核心数）
        lineage_mode: 使用谱系示踪比对模式
        cutsites: 结构配置（CutsiteRegion列表）

    返回:
        比对结果字典列表
    """
    total = len(queries)
    if total == 0:
        return []

    if threads is None or threads < 1:
        threads = mp.cpu_count()

    # 限制线程数不超过序列数
    threads = min(threads, total)

    mode_label = "谱系示踪比对" if lineage_mode else "标准比对"
    print(f"  使用 {threads} 个并行进程处理 {total} 条序列 ({mode_label})...")

    # 分批处理以减少进程间通信开销
    CHUNK_SIZE = 500
    chunks = [queries[i:i + CHUNK_SIZE] for i in range(0, total, CHUNK_SIZE)]
    print(f"  共 {len(chunks)} 个批次 (每批最多 {CHUNK_SIZE} 条)...")

    from functools import partial
    chunk_func = partial(
        _process_align_chunk, reference_seq,
        match_score=match_score, mismatch_penalty=mismatch_penalty,
        gap_open=gap_open, gap_extend=gap_extend, semi_global=semi_global,
        lineage_mode=lineage_mode, cutsites=cutsites,
        cutsite_gap_scale=cutsite_gap_scale, flank_gap_scale=flank_gap_scale,
        far_gap_scale=far_gap_scale, flank_width=flank_width,
        mismatch_density_threshold=mismatch_density_threshold,
        mutation_window=mutation_window
    )

    results = []
    completed = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(chunk_func, ch): ch for ch in chunks}

        for future in as_completed(futures):
            chunk_results = future.result()
            results.extend(chunk_results)
            for r in chunk_results:
                if "error" in r:
                    failed += 1
                completed += 1
            pct = completed / total * 100
            print(f"\r  进度: {completed}/{total} ({pct:.1f}%)", end="", file=sys.stderr)

    print(file=sys.stderr)
    # 过滤掉无法锚定/比对的序列，不出现在输出中
    filtered = [r for r in results if "error" not in r]
    n_discarded = len(results) - len(filtered)
    if n_discarded > 0:
        print(f"已丢弃 {n_discarded} 条无法锚定的序列")
    print(f"批量比对完成: {len(filtered)} 条有效结果")
    return filtered


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
# 可视化图表生成
# ═══════════════════════════════════════════════════════════════

def _img_to_b64(fig) -> str:
    """将matplotlib图形转换为base64 HTML可嵌入字符串"""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight', facecolor='white')
    buf.seek(0)
    img_b64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    return img_b64


def _gen_reads_distribution(results: List[Dict]) -> str:
    """Reads数分布直方图"""
    read_counts = [r.get("readCount", 1) for r in results]
    if not read_counts:
        return ""
    fig, ax = plt.subplots(figsize=(6, 3))
    max_rc = max(read_counts)
    bins = min(50, max_rc) if max_rc > 1 else 1
    ax.hist(read_counts, bins=bins, color='#4a6cf7', edgecolor='white', alpha=0.8)
    ax.set_xlabel('Reads per Unique Sequence')
    ax.set_ylabel('Sequence Count')
    ax.set_title('Reads Count Distribution')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_alignment_rate(results: List[Dict]) -> str:
    """比对成功率饼图"""
    successful = sum(1 for r in results if "error" not in r)
    failed = len(results) - successful
    if failed == 0 and successful == 0:
        return ""
    values = [successful, failed] if failed > 0 else [successful]
    labels = ['Successful', 'Failed'] if failed > 0 else ['Successful']
    colors = ['#4a6cf7', '#ff6b6b'] if failed > 0 else ['#4a6cf7']
    fig, ax = plt.subplots(figsize=(4, 3))
    wedges, texts, autotexts = ax.pie(
        values, labels=labels, autopct='%1.1f%%',
        colors=colors, startangle=90, textprops={'fontsize': 10}
    )
    for t in autotexts:
        t.set_fontweight('bold')
    ax.set_title('Alignment Rate')
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_length_distribution(results: List[Dict], ref_length: int = 0) -> str:
    """查询序列长度分布直方图"""
    lengths = []
    for r in results:
        if "error" in r:
            continue
        seq = r.get("seq", "")
        if seq:
            lengths.append(len(seq))
        else:
            stats = r.get("stats")
            if not stats:
                continue
            alen = stats.get("alignment_length", 0)
            if alen:
                lengths.append(int(alen))
    if not lengths:
        return ""
    fig, ax = plt.subplots(figsize=(6, 3))
    bins = min(30, max(lengths) - min(lengths) + 1) if max(lengths) > min(lengths) else 1
    ax.hist(lengths, bins=bins, color='#20c997', edgecolor='white', alpha=0.8)
    if ref_length > 0:
        ax.axvline(ref_length, color='#dc3545', linestyle='--', linewidth=1.5, label=f'Ref ({ref_length}bp)')
        ax.legend(fontsize=9)
    ax.set_xlabel('Sequence Length (bp)')
    ax.set_ylabel('Count')
    ax.set_title('Sequence Length Distribution')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_mutation_overview(report_data: dict) -> Tuple[str, str]:
    """突变总览：序列角度和Reads角度的饼图"""
    s = report_data["summary"]
    mutated_seq = s.get("mutated_sequences", 0)
    unmutated_seq = s.get("unmutated_sequences", 0)
    mutated_reads = s.get("mutated_reads", 0)
    total_reads = s.get("total_reads_successful", 0)
    unmutated_reads = total_reads - mutated_reads

    def make_pie(values, labels, title):
        # 过滤掉零值
        nz = [(v, l) for v, l in zip(values, labels) if v > 0]
        if len(nz) <= 1:
            return ""
        nz_values = [v for v, l in nz]
        nz_labels = [l for v, l in nz]
        fig, ax = plt.subplots(figsize=(4, 3))
        colors_list = ['#ff6b6b', '#adb5bd', '#20c997', '#4a6cf7']
        wedges, texts, autotexts = ax.pie(
            nz_values, labels=nz_labels, autopct='%1.1f%%',
            colors=colors_list[:len(nz_values)], startangle=90,
            textprops={'fontsize': 10}
        )
        for t in autotexts:
            t.set_fontweight('bold')
        ax.set_title(title)
        fig.tight_layout()
        return _img_to_b64(fig)

    seq_pie = make_pie(
        [mutated_seq, unmutated_seq],
        [f'Mutated ({mutated_seq})', f'Unmutated ({unmutated_seq})'],
        'Mutation Rate (by Sequence)'
    )
    reads_pie = make_pie(
        [mutated_reads, unmutated_reads],
        [f'Mutated ({mutated_reads})', f'Unmutated ({unmutated_reads})'],
        'Mutation Rate (by Reads)'
    )
    return seq_pie, reads_pie


def _gen_mutation_type_chart(report_data: dict) -> str:
    """突变类型分布柱状图（序列数和Reads数并排）"""
    mt = report_data["mutation_types"]
    labels = ['Only\nSubstitution', 'Only\nDeletion', 'Only\nInsertion',
              'Ins+Del', 'Ins+Sub', 'Del+Sub', 'All Three']
    keys = ['only_substitution', 'only_deletion', 'only_insertion',
            'insertion_and_deletion', 'insertion_and_substitution',
            'deletion_and_substitution', 'insertion_deletion_substitution']

    seq_counts = []; reads_counts = []
    for k in keys:
        v = mt.get(k, {"sequences": 0, "reads": 0})
        if isinstance(v, dict):
            seq_counts.append(v.get("sequences", 0))
            reads_counts.append(v.get("reads", 0))
        else:
            seq_counts.append(v); reads_counts.append(v)

    if all(s == 0 for s in seq_counts):
        return ""

    fig, ax = plt.subplots(figsize=(7, 3.5))
    x = range(len(labels))
    w = 0.35
    ax.bar([i - w/2 for i in x], seq_counts, w, label='Sequences', color='#4a6cf7', alpha=0.85)
    ax.bar([i + w/2 for i in x], reads_counts, w, label='Reads', color='#ff6b6b', alpha=0.85)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('Count')
    ax.set_title('Mutation Type Distribution')
    ax.legend(fontsize=9)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_indel_length_chart(report_data: dict) -> str:
    """Insertion和Deletion长度分布直方图并排"""
    ms = report_data.get("mutation_stats", {})
    ins_len = ms.get("avg_insertion_length", 0)
    del_len = ms.get("avg_deletion_length", 0)
    max_ins = ms.get("max_insertion_length", 0)
    max_del = ms.get("max_deletion_length", 0)

    detail = report_data.get("mutated_sequences_detail", [])
    ins_lengths = []; del_lengths = []
    for d in detail:
        for _ in range(d.get("gaps_in_ref", 0)):
            ins_lengths.append(1)
        for _ in range(d.get("gaps_in_query", 0)):
            del_lengths.append(1)
    if not ins_lengths and not del_lengths:
        return ""

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    if ins_lengths:
        axes[0].hist(ins_lengths, bins=min(20, max(ins_lengths)),
                     color='#20c997', edgecolor='white', alpha=0.8)
        axes[0].set_xlabel('Length (bp)'); axes[0].set_ylabel('Count')
        axes[0].set_title(f'Insertion (avg={ins_len:.1f}bp, max={max_ins}bp)')
        axes[0].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    else:
        axes[0].text(0.5, 0.5, 'No insertions', ha='center', va='center', transform=axes[0].transAxes)
        axes[0].set_title('Insertion')
    if del_lengths:
        axes[1].hist(del_lengths, bins=min(20, max(del_lengths)),
                     color='#ff6b6b', edgecolor='white', alpha=0.8)
        axes[1].set_xlabel('Length (bp)'); axes[1].set_ylabel('Count')
        axes[1].set_title(f'Deletion (avg={del_len:.1f}bp, max={max_del}bp)')
        axes[1].yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    else:
        axes[1].text(0.5, 0.5, 'No deletions', ha='center', va='center', transform=axes[1].transAxes)
        axes[1].set_title('Deletion')
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_indel_position_chart(results: List[Dict], ref_length: int) -> str:
    """基于参考序列位置的Indel分布柱状图"""
    if not results or ref_length == 0:
        return ""

    ins_positions = []; del_positions = []
    for r in results:
        if "error" in r:
            continue
        ar = r.get("aligned_ref", ""); aq = r.get("aligned_query", "")
        if not ar or not aq:
            continue
        ref_pos = 0
        for a, b in zip(ar, aq):
            if a == '-' and b != '-':
                ins_positions.append(ref_pos)
            elif a != '-' and b == '-':
                del_positions.append(ref_pos)
            if a != '-':
                ref_pos += 1

    if not ins_positions and not del_positions:
        return ""

    fig, ax = plt.subplots(figsize=(10, 3.5))
    bins = min(100, ref_length)
    if del_positions:
        ax.hist(del_positions, bins=bins, alpha=0.6, color='#ff6b6b',
                label=f'Deletion ({len(del_positions)})', edgecolor=None)
    if ins_positions:
        ax.hist(ins_positions, bins=bins, alpha=0.6, color='#20c997',
                label=f'Insertion ({len(ins_positions)})', edgecolor=None)
    ax.set_xlabel('Reference Position (bp)'); ax.set_ylabel('Indel Count')
    ax.set_title('Position-based Indel Distribution')
    ax.legend(fontsize=9)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_mutation_composition_chart(report_data: dict) -> str:
    """点突变/插入/删除总体比例图"""
    ms = report_data.get("mutation_stats", {})
    values = [ms.get("total_point_mutations", 0),
              ms.get("total_insertion_events", 0),
              ms.get("total_deletion_events", 0)]
    labels = ['Substitution', 'Insertion', 'Deletion']
    colors = ['#ffd43b', '#20c997', '#ff6b6b']
    # 过滤掉零值
    non_zero = [(v, l, c) for v, l, c in zip(values, labels, colors) if v > 0]
    if len(non_zero) <= 1:
        return ""
    fig, ax = plt.subplots(figsize=(5, 3))
    nz_values = [v for v, l, c in non_zero]
    nz_labels = [l for v, l, c in non_zero]
    nz_colors = [c for v, l, c in non_zero]
    wedges, texts, autotexts = ax.pie(
        nz_values, labels=nz_labels, autopct='%1.1f%%',
        colors=nz_colors, startangle=90, textprops={'fontsize': 10}
    )
    for t in autotexts:
        t.set_fontweight('bold')
    centre = plt.Circle((0, 0), 0.5, fc='white')
    ax.add_artist(centre)
    ax.set_title('Mutation Composition')
    fig.tight_layout()
    return _img_to_b64(fig)


def _gen_allele_heatmap(results: List[Dict], ref_seq: str,
                         cutsites: list = None, top_n: int = 50,
                         window_start: int = 0,
                         window_end: int = None) -> str:
    """
    Allele热图：每行一个allele，每列一个碱基位置。

    展示top_n个最频繁allele的碱基级比对，支持自定义显示区域。
    window_start/window_end = None 时显示全长参考序列。

    参数:
        results: 比对结果列表
        ref_seq: 参考序列
        cutsites: 用于标注cutsite区域
        top_n: 展示的top allele数
        window_start: 显示起始位置（0-indexed，默认0）
        window_end: 显示结束位置（含，默认ref_seq全长）

    返回: base64 PNG
    """
    if not results or not ref_seq:
        return ""

    ref_len = len(ref_seq)
    if window_end is None:
        window_end = ref_len - 1
    window_end = min(window_end, ref_len - 1)
    window_len = window_end - window_start + 1
    if window_len <= 0:
        return ""

    # 1. 聚合allele
    allele_counts = {}
    allele_data = {}
    for r in results:
        if "error" in r:
            continue
        aq = r.get("aligned_query", "")
        if not aq:
            continue
        rc = r.get("readCount", 1)
        allele_counts[aq] = allele_counts.get(aq, 0) + rc
        if aq not in allele_data:
            allele_data[aq] = {
                "aligned_ref": r.get("aligned_ref", ""),
                "readCount": rc
            }
    if not allele_counts:
        return ""

    # 2. 排序取top_n
    sorted_alleles = sorted(allele_counts.items(), key=lambda x: -x[1])
    top_alleles = sorted_alleles[:top_n]
    total_reads = sum(allele_counts.values())

    # 3. 从对齐序列中提取窗口碱基
    # 先找一个非错误的 aligned_ref 来确定列范围
    sample_aligned_ref = ""
    for r in results:
        if "error" not in r:
            sample_aligned_ref = r.get("aligned_ref", "")
            if sample_aligned_ref:
                break
    if not sample_aligned_ref:
        return ""

    # 找到窗口对应的列范围
    ref_pos = 0
    col_start, col_end = None, None
    for col, c in enumerate(sample_aligned_ref):
        if c != '-':
            if ref_pos == window_start:
                col_start = col
            if ref_pos == window_end:
                col_end = col
                break
            ref_pos += 1
    if col_start is None:
        col_start = 0
    if col_end is None:
        col_end = min(col_start + window_len * 2, len(sample_aligned_ref) - 1)

    n_cols = col_end - col_start + 1

    def slice_align(align_str):
        if col_end <= len(align_str):
            return align_str[col_start:col_end + 1]
        s = align_str[col_start:]
        while len(s) < n_cols:
            s += '-'
        return s[:n_cols]

    # 窗口参考碱基
    window_ref_bases = ref_seq[window_start:window_end + 1]

    # 4. 碱基配色
    base_colors = {
        'A': '#cbe9cb', 'C': '#fee5ce', 'G': '#ffffd6',
        'T': '#e5deed', 'U': '#e5deed', 'N': '#e9e9e9',
        'a': '#cbe9cb', 'c': '#fee5ce', 'g': '#ffffd6',
        't': '#e5deed', 'u': '#e5deed', 'n': '#e9e9e9',
        '-': '#f0f0f0',
    }

    n_rows = len(top_alleles)

    # 5. 自适应单元格尺寸
    if window_len > 200:
        cell_size = 10
        font_size = 5.5
    elif window_len > 100:
        cell_size = 14
        font_size = 6.5
    elif window_len > 60:
        cell_size = 18
        font_size = 7
    else:
        cell_size = 22
        font_size = 9

    label_w = 180  # 右侧标注宽度(px)
    fig_w_px = n_cols * cell_size + label_w
    fig_h_px = (n_rows + 5.0) * cell_size

    dpi = 100
    fig_w = max(8, fig_w_px / dpi)
    fig_h = max(3, fig_h_px / dpi)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, n_cols)
    ax.set_ylim(-2.5, n_rows + 2.5)
    ax.set_aspect('equal')
    ax.axis('off')

    # 6. 参考序列行
    ref_y = n_rows + 1.5
    for ci in range(n_cols):
        base = window_ref_bases[ci] if ci < len(window_ref_bases) else ' '
        x = ci + 0.5
        color = base_colors.get(base, '#e9e9e9')
        ax.add_patch(plt.Rectangle((ci, ref_y - 0.5), 1, 1,
                                    facecolor=color, edgecolor='white', linewidth=0.3))
        ax.text(x, ref_y, base, ha='center', va='center',
                fontsize=font_size, fontweight='normal')

    # 标注cutsite区域（灰色半透明方框）
    if cutsites:
        for cs in cutsites:
            if cs.start > window_end or cs.end < window_start:
                continue
            cs_c0 = max(cs.start - window_start, 0)
            cs_c1 = min(cs.end - window_start, n_cols - 1)
            if cs_c0 <= cs_c1:
                ax.add_patch(plt.Rectangle(
                    (cs_c0, ref_y - 0.5), cs_c1 - cs_c0 + 1, 1,
                    facecolor='#888888', alpha=0.3, edgecolor=None, zorder=0
                ))

    # 8. sgRNA区域示意（半高灰色方块，紧贴reference下方）
    if cutsites:
        # 建立 ref_pos → display column 映射
        ref_pos_to_col = {}
        rp = window_start
        for ci in range(n_cols):
            if sample_aligned_ref[col_start + ci] != '-':
                ref_pos_to_col[rp] = ci
                rp += 1

        grna_y = ref_y - 0.75  # 半高行中心
        for cs in cutsites:
            if not cs.name.startswith("Target"):
                continue
            tgt_start = cs.start - 13
            tgt_end = cs.end
            if tgt_start > window_end or tgt_end < window_start:
                continue
            col0 = ref_pos_to_col.get(max(tgt_start, window_start))
            col1 = ref_pos_to_col.get(min(tgt_end, window_end))
            if col0 is None or col1 is None:
                continue
            # 绘制灰色方块
            ax.add_patch(plt.Rectangle(
                (col0, grna_y - 0.3), col1 - col0 + 1, 0.6,
                facecolor='#bbbbbb', edgecolor='#999999', linewidth=0.3, zorder=0
            ))
            # 在前方标注Target名称
            ax.text(col0 - 0.1, grna_y, cs.name.replace('Target', 'T'),
                    ha='right', va='center', fontsize=max(6, font_size - 0.5),
                    color='#666', fontweight='bold')

    # 7. Allele行
    for ai, (aligned_query, total_rc) in enumerate(top_alleles):
        y = n_rows - 1 - ai
        ad = allele_data.get(aligned_query, {})
        aligned_ref = ad.get("aligned_ref", "")
        win_bases = slice_align(aligned_query)
        win_ref_b = slice_align(aligned_ref) if aligned_ref else window_ref_bases

        while len(win_bases) < n_cols:
            win_bases += '-'
        while len(win_ref_b) < n_cols:
            win_ref_b += '-'
        win_bases = win_bases[:n_cols]
        win_ref_b = win_ref_b[:n_cols]

        for ci in range(n_cols):
            base = win_bases[ci]
            ref_base = win_ref_b[ci]
            is_del = (base == '-')
            is_ins = (base != '-' and ref_base == '-')
            is_sub = (not is_del and not is_ins and base != ref_base)

            if is_ins:
                color = '#ffcccc'
            elif is_del:
                color = '#f0f0f0'
            else:
                color = base_colors.get(base, '#e9e9e9')

            ax.add_patch(plt.Rectangle((ci, y - 0.5), 1, 1,
                                        facecolor=color, edgecolor='white', linewidth=0.2))

            if is_ins:
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, color='red', fontweight='bold')
                ax.add_patch(plt.Rectangle(
                    (ci, y - 0.5), 1, 1,
                    facecolor='none', edgecolor='red', linewidth=1.2, zorder=5
                ))
            elif is_del:
                ax.text(ci + 0.5, y, '-', ha='center', va='center',
                        fontsize=font_size, color='#666', fontweight='bold')
            else:
                w = 'bold' if is_sub else 'normal'
                ax.text(ci + 0.5, y, base, ha='center', va='center',
                        fontsize=font_size, fontweight=w)

        # 右侧label
        pct = total_rc / total_reads * 100
        ax.text(n_cols + 0.3, y, f"{pct:.1f}% ({total_rc})",
                ha='left', va='center', fontsize=max(8, font_size + 1))

    # 图例
    lgy = -1.2
    items = [('A', '#cbe9cb'), ('C', '#fee5ce'), ('G', '#ffffd6'),
             ('T', '#e5deed'), ('-', '#f0f0f0'), ('ins', '#ffcccc')]
    lx = 0
    leg_font = max(9, font_size + 2)
    for label, color in items:
        ax.add_patch(plt.Rectangle((lx, lgy - 0.35), 0.7, 0.7,
                                    facecolor=color, edgecolor='#666', linewidth=0.6))
        ax.text(lx + 0.9, lgy, label, ha='left', va='center', fontsize=leg_font)
        lx += 2.2

    ax.set_title(f'Top {top_n} Alleles — ref {window_start + 1}–{window_end + 1}bp ({window_len}bp)',
                 fontsize=max(10, font_size + 3), pad=4)
    fig.tight_layout()
    return _img_to_b64(fig)


def generate_charts(results: List[Dict], report_data: dict = None,
                     ref_length: int = 0, ref_seq: str = "",
                     cutsites: list = None,
                     allele_window_start: int = 0,
                     allele_window_end: int = None,
                     allele_top_n: int = 50) -> Dict[str, str]:
    """
    生成所有图表，返回 {图表名: base64图片} 字典。

    参数:
        results: 比对结果列表
        report_data: 报告数据字典
        ref_length: 参考序列长度
        ref_seq: 参考序列字符串（用于allele热图）
        cutsites: CutsiteRegion列表（用于allele热图标注）
        allele_window_start: Allele热图显示起始位置（默认0）
        allele_window_end: Allele热图显示结束位置（含，默认全长）
        allele_top_n: Allele热图展示的top allele数（默认50）
    """
    charts = {}
    if not _HAS_MPL:
        return charts
    try:
        charts['reads_dist'] = _gen_reads_distribution(results)
        charts['align_rate'] = _gen_alignment_rate(results)
        charts['length_dist'] = _gen_length_distribution(results, ref_length)
        if report_data:
            seq_pie, reads_pie = _gen_mutation_overview(report_data)
            charts['mutation_seq_pie'] = seq_pie
            charts['mutation_reads_pie'] = reads_pie
            charts['mutation_type'] = _gen_mutation_type_chart(report_data)
            charts['indel_length'] = _gen_indel_length_chart(report_data)
            charts['indel_position'] = _gen_indel_position_chart(results, ref_length)
            charts['mutation_composition'] = _gen_mutation_composition_chart(report_data)
        # Allele热图（默认显示全长，窗口可调）
        if ref_seq and ref_length > 0:
            ws = allele_window_start
            we = allele_window_end if allele_window_end is not None else ref_length - 1
            charts['allele_heatmap'] = _gen_allele_heatmap(
                results, ref_seq, cutsites=cutsites, top_n=allele_top_n,
                window_start=ws, window_end=we
            )
    except Exception as e:
        print(f"警告: 图表生成失败 - {e}", file=sys.stderr)
    return charts


# ═══════════════════════════════════════════════════════════════
# 分析报告生成
# ═══════════════════════════════════════════════════════════════

def generate_report(results: List[Dict], output_path: str, format: str = "json",
                     ref_length: int = 0, ref_seq: str = "",
                     cutsites: list = None,
                     allele_window_start: int = 0,
                     allele_window_end: int = None,
                     allele_top_n: int = 50) -> None:
    """
    生成突变分析报告

    参数:
        results: 比对结果列表
        output_path: 输出文件路径
        format: 报告格式 (json, html)
        ref_length: 参考序列长度（用于图表）
        ref_seq: 参考序列（用于allele热图）
        cutsites: cutsite区域列表（用于allele热图标注）
        allele_window_start: Allele热图显示起始位置
        allele_window_end: Allele热图显示结束位置（含）
        allele_top_n: Allele热图展示的top allele数（默认50）
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
            ins_lengths.append(block)

        # 收集删除长度
        for block in stats.get("gap_blocks_query", []):
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
        charts = generate_charts(results, report, ref_length,
                                  ref_seq=ref_seq, cutsites=cutsites,
                                  allele_window_start=allele_window_start,
                                  allele_window_end=allele_window_end,
                                  allele_top_n=allele_top_n)
        _save_report_html(report, path, charts)
        print(f"HTML分析报告已保存至 {path}")

    else:
        print(f"错误: 不支持的报告格式: {format}", file=sys.stderr)
        sys.exit(1)


def _save_report_html(report: dict, output_path: str, charts: dict = None) -> None:
    """将报告保存为自包含的HTML文件（可选嵌入图表）"""
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

    # 图表HTML生成辅助（可点击放大）
    def _img_html(key):
        b64 = (charts or {}).get(key, '')
        if not b64:
            return '<p style="color:#999;font-style:italic;">图表不可用（需安装matplotlib）</p>'
        klass = 'chart-img' if key != 'allele_heatmap' else 'chart-img allele-img'
        return f'<div class="chart-scroll"><img class="{klass}" src="data:image/png;base64,{b64}" onclick="openModal(this)" /></div>'

    chart_reads_dist = _img_html('reads_dist')
    chart_align_rate = _img_html('align_rate')
    chart_length_dist = _img_html('length_dist')
    chart_mut_seq = _img_html('mutation_seq_pie')
    chart_mut_reads = _img_html('mutation_reads_pie')
    chart_mut_type = _img_html('mutation_type')
    chart_mut_comp = _img_html('mutation_composition')
    chart_indel_len = _img_html('indel_length')
    chart_indel_pos = _img_html('indel_position')
    chart_allele = _img_html('allele_heatmap')

    has_charts = bool(charts)

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
  .chart-box {{ background: #fff; border-radius: 10px; padding: 15px; box-shadow: 0 2px 8px rgba(0,0,0,0.08); margin: 15px 0; text-align: center; }}
  .chart-box .chart-img {{ max-width: 100%; height: auto; border-radius: 6px; cursor: zoom-in; }}
  .chart-row {{ display: flex; flex-wrap: wrap; gap: 15px; margin: 15px 0; }}
  .chart-row .chart-box {{ flex: 1; min-width: 300px; }}
  footer {{ text-align: center; color: #aaa; font-size: 12px; margin-top: 40px; padding-top: 20px; border-top: 1px solid #ddd; }}
  /* scrollable container */
  .chart-scroll {{ overflow-x: auto; overflow-y: hidden; max-width: 100%; }}
  .chart-scroll .allele-img {{ max-width: none; cursor: zoom-in; }}
  /* modal overlay */
  .modal-overlay {{ display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.85); z-index: 9999; cursor: zoom-out; }}
  .modal-overlay.active {{ display: flex; align-items: center; justify-content: center; }}
  .modal-content {{ position: relative; max-width: 95vw; max-height: 95vh; overflow: auto; cursor: default; background: transparent; }}
  .modal-content img {{ display: block; max-width: none; width: auto; height: auto; border-radius: 4px; box-shadow: 0 4px 30px rgba(0,0,0,0.3); }}
  .modal-close {{ position: absolute; top: -36px; right: 0; color: #fff; font-size: 28px; cursor: pointer; font-weight: bold; line-height: 1; background: none; border: none; padding: 4px 12px; border-radius: 4px; }}
  .modal-close:hover {{ background: rgba(255,255,255,0.15); }}
  .modal-hint {{ position: absolute; bottom: -28px; left: 0; color: rgba(255,255,255,0.5); font-size: 13px; }}
</style>
<script>
  function openModal(imgEl) {{
    var overlay = document.getElementById('modalOverlay');
    var modalImg = document.getElementById('modalImg');
    modalImg.src = imgEl.src;
    modalImg.style.transform = 'scale(1)';
    modalImg.style.marginLeft = '0px';
    modalImg.style.marginTop = '0px';
    overlay.classList.add('active');
    document.body.style.overflow = 'hidden';
  }}
  function closeModal() {{
    document.getElementById('modalOverlay').classList.remove('active');
    document.body.style.overflow = '';
  }}
  document.addEventListener('DOMContentLoaded', function() {{
    var overlay = document.getElementById('modalOverlay');
    var modalImg = document.getElementById('modalImg');
    var scale = 1;
    overlay.addEventListener('click', function(e) {{
      if (e.target === overlay) closeModal();
    }});
    overlay.addEventListener('wheel', function(e) {{
      e.preventDefault();
      var delta = e.deltaY > 0 ? -0.1 : 0.1;
      scale = Math.max(0.2, Math.min(10, scale + delta));
      modalImg.style.transform = 'scale(' + scale + ')';
    }}, {{ passive: false }});
    overlay.addEventListener('dblclick', function() {{
      scale = 1;
      modalImg.style.transform = 'scale(1)';
      modalImg.style.marginLeft = '0px';
      modalImg.style.marginTop = '0px';
    }});
    var isDragging = false, startX, startY, origX = 0, origY = 0;
    modalImg.addEventListener('mousedown', function(e) {{
      if (e.button !== 0) return;
      isDragging = true;
      startX = e.clientX - origX;
      startY = e.clientY - origY;
      modalImg.style.cursor = 'grabbing';
      e.preventDefault();
    }});
    document.addEventListener('mousemove', function(e) {{
      if (!isDragging) return;
      origX = e.clientX - startX;
      origY = e.clientY - startY;
      modalImg.style.marginLeft = origX + 'px';
      modalImg.style.marginTop = origY + 'px';
    }});
    document.addEventListener('mouseup', function() {{
      isDragging = false;
      modalImg.style.cursor = 'default';
    }});
  }});
</script>
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

<!-- 图表：数据总览 -->
{has_charts and '<h2>Reads & Alignment Overview</h2>' or ''}
<div class="chart-row">
  <div class="chart-box">{chart_reads_dist}</div>
  <div class="chart-box">{chart_align_rate}</div>
</div>
<div class="chart-box">{chart_length_dist}</div>

<h2>📊 突变类型分布</h2>
<table>
<tr><th>突变类型</th><th>序列数</th><th>序列占比</th><th>Reads数</th><th>Reads占比</th></tr>
{type_rows}
</table>

<!-- 图表：突变总览 -->
{has_charts and '<h2>Mutation Overview</h2>' or ''}
<div class="chart-row">
  <div class="chart-box">{chart_mut_seq}</div>
  <div class="chart-box">{chart_mut_reads}</div>
</div>
{has_charts and '<div class="chart-row">' or ''}
  <div class="chart-box">{chart_mut_type}</div>
  <div class="chart-box">{chart_mut_comp}</div>
{has_charts and '</div>' or ''}

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

{has_charts and '<h2>Indel Length & Position Distribution</h2>' or ''}
<div class="chart-box">{chart_indel_len}</div>
<div class="chart-box">{chart_indel_pos}</div>

{has_charts and '<h2>Allele Heatmap (Top 20)</h2>' or ''}
<div class="chart-box">{chart_allele}</div>

<h2>📋 突变序列明细（前100条）</h2>
<table>
<tr><th>#</th><th>序列名</th><th>Reads</th><th>点突变</th><th>插入事件</th><th>删除事件</th><th>相似度</th></tr>
{detail_rows}
</table>

<footer>Generated by CARLIN 序列分析工具 v{report["version"]}</footer>
</div>

<!-- Modal -->
<div id="modalOverlay" class="modal-overlay">
  <div class="modal-content" onclick="event.stopPropagation()">
    <button class="modal-close" onclick="closeModal()">&times;</button>
    <img id="modalImg" src="" />
    <div class="modal-hint">滚轮缩放 · 拖拽平移 · 双击复原</div>
  </div>
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

    # 谱系示踪模式 ──────────────────────────────────────────
    align_parser.add_argument("--lineage", action="store_true",
                              help="启用谱系示踪比对模式（结构感知，优先在cutsite开gap，自动过滤假阳性突变）")
    align_parser.add_argument("--cutsite-scale", type=float, default=1.0,
                              help="cutsite区域gap惩罚倍率（默认：1.0，越小越容易开gap）")
    align_parser.add_argument("--flank-scale", type=float, default=2.0,
                              help="cutsite侧翼区gap惩罚倍率（默认：2.0）")
    align_parser.add_argument("--far-scale", type=float, default=2.0,
                              help="远离cutsite区域gap惩罚倍率（默认：2.0，越大越难开gap）")
    align_parser.add_argument("--flank-width", type=int, default=3,
                              help="cutsite侧翼范围bp（默认：3）")
    align_parser.add_argument("--mutation-window", type=int, default=3,
                              help="保留点突变的cutsite窗口半径bp（默认：3，即cutsite±3bp内保留点突变）")
    align_parser.add_argument("--density-threshold", type=float, default=0.34,
                              help="mismatch密度阈值，超此阈值视为insertion（默认：0.34）")
    align_parser.add_argument("--cutsites", default=None,
                              help="cutsite位置配置文件（JSON格式），不指定则自动推断标准CARLIN结构")

    # 并行参数
    align_parser.add_argument("--threads", "-t", type=int, default=None,
                              help="并行进程数（默认：自动使用所有CPU核心）")

    # 报告参数
    align_parser.add_argument("--report", choices=["json", "html"], default=None,
                              help="生成分析报告（json或html格式）")
    align_parser.add_argument("--report-output", default=None,
                              help="报告输出路径（默认：基于--output文件名）")
    align_parser.add_argument("--allele-window-start", type=int, default=0,
                              help="Allele热图显示起始位置（0-indexed，默认：0）")
    align_parser.add_argument("--allele-window-end", type=int, default=None,
                              help="Allele热图显示结束位置（含，0-indexed，默认：全长）")
    align_parser.add_argument("--allele-top-n", type=int, default=50,
                              help="Allele热图展示的top allele数（默认：50）")

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

        # 尝试自动检测cutsite区域（用于allele热图的gRNA标注，不依赖比对模式）
        display_cutsites = get_amplicon_structure(reference_seq)
        if display_cutsites:
            print(f"自动检测到 {len(display_cutsites)} 个cutsite区域（用于热图标注）")
        else:
            display_cutsites = None

        # 结构感知比对模式
        cutsites = None
        if args.lineage:
            if args.cutsites:
                # 从JSON文件读取cutsite位置
                try:
                    import json as _json
                    with open(args.cutsites) as f:
                        cs_data = _json.load(f)
                    cutsites = [
                        CutsiteRegion(name=c.get("name", f"Target{i+1}"),
                                      start=c["start"], end=c["end"])
                        for i, c in enumerate(cs_data.get("cutsites", cs_data))
                    ]
                    print(f"从配置文件加载 {len(cutsites)} 个cutsite区域")
                except Exception as e:
                    print(f"错误: 读取cutsite配置文件失败 - {e}", file=sys.stderr)
                    sys.exit(1)
            else:
                cutsites = display_cutsites
                if not cutsites:
                    print("错误: 无法自动推断cutsite位置，请使用 --cutsites 手动指定", file=sys.stderr)
                    sys.exit(1)

            print(f"  谱系示踪参数: cutsite倍率={args.cutsite_scale}, "
                  f"侧翼倍率={args.flank_scale}, 远端倍率={args.far_scale}")
            print(f"  突变窗口: cutsite±{args.mutation_window}bp, "
                  f"mismatch密度阈值: {args.density_threshold}")

        # 运行并行批量比对
        if args.use_global:
            align_type = "全局"
        else:
            align_type = "半全局"

        if args.threads:
            print(f"开始并行批量比对（{align_type}，{args.threads} 线程）...")
        else:
            cpu_count = mp.cpu_count()
            print(f"开始并行批量比对（{align_type}，自动检测到 {cpu_count} 个CPU核心）...")

        results = parallel_align(
            reference_seq, queries,
            match_score=args.match_score,
            mismatch_penalty=args.mismatch_penalty,
            gap_open=args.gap_open,
            gap_extend=args.gap_extend,
            semi_global=not args.use_global,
            threads=args.threads,
            lineage_mode=args.lineage,
            cutsites=cutsites,
            cutsite_gap_scale=args.cutsite_scale,
            flank_gap_scale=args.flank_scale,
            far_gap_scale=args.far_scale,
            flank_width=args.flank_width,
            mismatch_density_threshold=args.density_threshold,
            mutation_window=args.mutation_window
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
            generate_report(results, report_output, args.report,
                             ref_length=len(reference_seq),
                             ref_seq=reference_seq,
                             cutsites=display_cutsites,
                             allele_window_start=args.allele_window_start,
                             allele_window_end=args.allele_window_end,
                             allele_top_n=args.allele_top_n)

    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
