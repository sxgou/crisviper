"""I/O functions for FASTQ/TSV/FASTA file handling.

Includes single-cell 10x and InDrops/BGI FASTQ parsing
(MATLAB equivalents from @SCFastQData).
"""

import sys
import os
import gzip
import csv
import shutil
import subprocess
import tempfile
from typing import List, Dict, Optional, Tuple
from collections import Counter
from Bio import SeqIO
from crisviper.logging_config import get_logger

log = get_logger(__name__)


def _read_fastq_counts(fastq_path: str) -> Dict[str, int]:
    """读取FASTQ文件，返回 {序列: 出现次数} 字典"""
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
    return counts


def fastq_to_dataframe(fastq_path: str, sample_name: str = "sample") -> List[Dict]:
    """
    将FASTQ文件转换为字典列表，每个字典代表一条唯一序列

    参数:
        fastq_path: FASTQ文件路径（支持.gz）
        sample_name: 样本名称

    返回:
        字典列表，每个字典包含readName, cellBC, UMI, readCount, seq字段
    """
    counts = _read_fastq_counts(fastq_path)

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

    log.info("从 %s 中读取了 %d 条唯一序列（总计 %d 条reads）。", fastq_path, len(rows), sum(counts.values()))
    return rows


def fastq_to_fasta(fastq_path: str, output_fasta: str, sample_name: str = "sample") -> None:
    """
    将FASTQ文件转换为FASTA格式，头部包含元数据

    参数:
        fastq_path: 输入FASTQ文件（支持.gz）
        output_fasta: 输出FASTA文件路径
        sample_name: 样本名称
    """
    counts = _read_fastq_counts(fastq_path)

    with open(output_fasta, 'w') as f:
        for i, (seq, count) in enumerate(counts.items()):
            read_name = f"{sample_name}_seq{i+1}"
            f.write(f">{read_name} cellBC={sample_name} UMI=UMI{i+1} readCount={count}\n{seq}\n")

    total_reads = sum(counts.values())
    log.info("已将 %d 条唯一序列（%d 条reads）写入 %s", len(counts), total_reads, output_fasta)


def fastq_to_fasta_from_rows(rows: List[Dict], output_path: str) -> None:
    """将字典列表（merge_paired_end/fastq_to_dataframe 输出格式）写入FASTA。

    参数:
        rows: 字典列表，包含 readName, cellBC, UMI, readCount, seq
        output_path: 输出FASTA文件路径
    """
    with open(output_path, 'w') as f:
        for row in rows:
            f.write(f">{row['readName']} cellBC={row['cellBC']} "
                    f"UMI={row['UMI']} readCount={row['readCount']}\n"
                    f"{row['seq']}\n")
    log.info("数据已保存至 %s", output_path)


def _check_fastp() -> None:
    """检查fastp是否可用，否则退出。"""
    if shutil.which("fastp") is None:
        log.error("fastp 未安装。双端合并需要 fastp，请执行: "
                  "brew install fastp  # macOS 或 "
                  "conda install -c bioconda fastp  # conda")
        sys.exit(1)


def merge_paired_end(
    fastq1: str,
    fastq2: str,
    min_overlap: int = 10,
    max_mismatch_rate: int = 20,
    max_mismatch_diff: int = 5,
    require_qual: int = 15,
    sample_name: str = "sample",
) -> List[Dict]:
    """调用 fastp 合并双端测序 reads，返回去重后的字典列表。

    参数:
        fastq1: R1 FASTQ 文件路径
        fastq2: R2 FASTQ 文件路径
        min_overlap: 最小 overlap 长度（fastp: overlap_len_require）
        max_mismatch_rate: Overlap 错配比例上限 %（fastp: overlap_diff_percent_limit）
        max_mismatch_diff: Overlap 最大错配绝对数（fastp: overlap_diff_limit）
        require_qual: 碱基质量阈值（fastp: -q）
        sample_name: 样本名称

    返回:
        字典列表（与 fastq_to_dataframe 格式一致）
    """
    _check_fastp()

    with tempfile.TemporaryDirectory() as tmpdir:
        merged_out = os.path.join(tmpdir, "merged.fastq")
        unmerged1 = os.path.join(tmpdir, "unmerged_R1.fastq")
        unmerged2 = os.path.join(tmpdir, "unmerged_R2.fastq")

        cmd = [
            "fastp",
            "-i", fastq1,
            "-I", fastq2,
            "-m",
            "--merged_out", merged_out,
            "--out1", unmerged1,
            "--out2", unmerged2,
            "--overlap_len_require", str(min_overlap),
            "--overlap_diff_percent_limit", str(max_mismatch_rate),
            "--overlap_diff_limit", str(max_mismatch_diff),
            "-q", str(require_qual),
            "-j", "/dev/null",
            "-h", "/dev/null",
        ]

        log.info("运行 fastp 合并双端 reads...")
        log.debug("  cmd: %s", " ".join(cmd))

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if result.stdout:
                for line in result.stdout.splitlines():
                    log.debug("  fastp: %s", line)
        except subprocess.CalledProcessError as e:
            log.error("fastp 运行失败 (exit %d)", e.returncode)
            if e.stderr:
                for line in e.stderr.splitlines():
                    log.error("  %s", line)
            sys.exit(1)

        # 读取合并结果
        counts = {}
        if not os.path.exists(merged_out) or os.path.getsize(merged_out) == 0:
            log.warning("fastp 未产生任何合并 reads！")
            return []

        with open(merged_out) as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = str(record.seq)
                counts[seq] = counts.get(seq, 0) + 1

        if not counts:
            log.warning("fastp 未产生任何合并 reads！")
            return []

        rows = []
        for i, (seq, count) in enumerate(counts.items()):
            rows.append({
                "readName": f"{sample_name}_seq{i+1}",
                "cellBC": sample_name,
                "UMI": f"UMI{i+1}",
                "readCount": count,
                "seq": seq,
            })

        log.info("fastp 合并完成: %d 条唯一序列（总计 %d 条reads）",
                 len(rows), sum(counts.values()))
        return rows


def save_tsv(rows: List[Dict], output_path: str) -> None:
    """将字典列表保存为TSV文件"""
    fieldnames = ["readName", "cellBC", "UMI", "readCount", "seq"]
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    log.info("数据已保存至 %s", output_path)


def read_reference_fasta(fasta_path: str) -> str:
    """读取reference FASTA文件，返回序列字符串"""
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            return str(record.seq).upper()
    except Exception as e:
        log.error("读取reference文件失败 - %s", e)
        sys.exit(1)
    log.error("reference文件为空: %s", fasta_path)
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
                    log.warning("第%d行缺少序列字段，跳过", line_num)
                    continue
                row['seq'] = row['seq'].upper()
                row['readCount'] = int(row.get('readCount', 1))
                rows.append(row)
    except Exception as e:
        log.error("读取TSV文件失败 - %s", e)
        sys.exit(1)

    if not rows:
        log.error("TSV文件为空或不包含有效序列: %s", tsv_path)
        sys.exit(1)

    log.info("从 %s 中读取了 %d 条查询序列。", tsv_path, len(rows))
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
        log.error("读取FASTA文件失败 - %s", e)
        sys.exit(1)

    if not rows:
        log.error("FASTA文件为空或格式不正确: %s", fasta_path)
        sys.exit(1)

    log.info("从 %s 中读取了 %d 条查询序列。", fasta_path, len(rows))
    return rows


# ═══════════════════════════════════════════════════════════════
# 单细胞 FASTQ 解析 (10x / InDrops / BGI)
# ═══════════════════════════════════════════════════════════════

def parse_10x_provenance(
    cb_reads: List[str],
    qc_reads: List[str],
    cb_length: int = 16,
    umi_length: int = 12,
) -> Tuple[List[str], List[str], List[str]]:
    """Parse 10x cell barcode and UMI from provenance read (R1).

    MATLAB equivalent: @SCFastQData/parse_10x_provenance

    Extracts CB from the first `cb_length` bases and UMI from the
    remaining bases (up to cb_length + umi_length) of each R1 read.

    Args:
        cb_reads: List of R1 read sequences (CB+UMI).
        qc_reads: List of R1 quality strings.
        cb_length: Length of cell barcode (default 16 for 10x).
        umi_length: Length of UMI (default 12 for 10x).

    Returns:
        Tuple of (cb_list, umi_list, qc_list).
    """
    prov_len = cb_length + umi_length
    cbs = []
    umis = []
    qcs = []
    for cb, qc in zip(cb_reads, qc_reads):
        cbs.append(cb[:min(len(cb), cb_length)])
        umis.append(cb[min(len(cb), cb_length):min(len(cb), prov_len)])
        qcs.append(qc[:min(len(qc), prov_len)])
    return cbs, umis, qcs


def parse_10x_fastq(
    r1_fastq: str,
    r2_fastq: str,
    cb_length: int = 16,
    umi_length: int = 12,
) -> Dict:
    """Parse 10x paired FASTQ files (R1=CB/UMI, R2=sequence).

    MATLAB equivalent: @SCFastQData/parse_10x_fastq

    Args:
        r1_fastq: R1 FASTQ path (cell barcode + UMI).
        r2_fastq: R2 FASTQ path (genomic sequence).
        cb_length: Cell barcode length in bp.
        umi_length: UMI length in bp.

    Returns:
        Dict with keys: CB, read_CB, UMI, read_UMI, SEQ, read_SEQ, QC, Nreads.
    """
    r1_records = list(SeqIO.parse(_maybe_gzopen(r1_fastq), "fastq"))
    r2_records = list(SeqIO.parse(_maybe_gzopen(r2_fastq), "fastq"))

    assert len(r1_records) == len(r2_records), \
        f"R1 ({len(r1_records)}) and R2 ({len(r2_records)}) have different read counts"

    N = len(r1_records)
    raw_cb = [str(r.seq) for r in r1_records]
    raw_qc = [str(r.letter_annotations["phred_quality"]) for r in r1_records]
    raw_seq = [str(r.seq) for r in r2_records]

    # Parse CB and UMI from R1
    cbs, umis, qcs = parse_10x_provenance(raw_cb, raw_qc, cb_length, umi_length)

    # Frequency-based dedup
    cb_counter = Counter(cbs)
    cb_unique = list(cb_counter.keys())
    cb_read_idx = [cb_unique.index(c) for c in cbs]

    umi_counter = Counter(umis)
    umi_unique = list(umi_counter.keys())
    umi_read_idx = [umi_unique.index(u) for u in umis]

    seq_counter = Counter(raw_seq)
    seq_unique = list(seq_counter.keys())
    seq_read_idx = [seq_unique.index(s) for s in raw_seq]

    return {
        "CB": cb_unique,
        "read_CB": cb_read_idx,
        "UMI": umi_unique,
        "read_UMI": umi_read_idx,
        "SEQ": seq_unique,
        "read_SEQ": seq_read_idx,
        "QC": qcs,
        "Nreads": N,
    }


def filter_sc_cbs_and_umis(
    cbs: List[str],
    read_cb: List[int],
    umis: List[str],
    read_umi: List[int],
    qcs: List[str],
    cb_length: int = 16,
    umi_length: int = 12,
    platform: str = "10x",
    min_qscore: int = 20,
) -> Dict:
    """Filter single-cell CBs and UMIs by quality metrics.

    MATLAB equivalent: @SCFastQData/filter_sc_CBs_and_UMIs

    Filters based on:
    - QC length matches (CB + UMI length)
    - CB and UMI correct length
    - No 'N' bases in CB or UMI
    - Quality score >= min_qscore

    Args:
        cbs: List of unique cell barcodes.
        read_cb: CB index per read.
        umis: List of unique UMIs.
        read_umi: UMI index per read.
        qcs: List of quality strings per read.
        cb_length: Expected CB length.
        umi_length: Expected UMI length.
        platform: '10x' or 'inDrops'.
        min_qscore: Minimum Phred quality score.

    Returns:
        Dict of boolean masks (as index arrays).
    """
    N = len(read_cb)
    masks_raw = {}

    # QC length match: len(CB) + len(UMI) == len(QC)
    cb_by_read = [cbs[i] for i in read_cb]
    umi_by_read = [umis[i] for i in read_umi]
    qc_lens = [len(q) for q in qcs]
    cb_umi_lens = [len(cb_by_read[i]) + len(umi_by_read[i]) for i in range(N)]
    masks_raw["QC_length_match"] = [cb_umi_lens[i] == qc_lens[i] for i in range(N)]

    # CB correct length
    masks_raw["CB_correct_length"] = [len(cb) == cb_length for cb in cb_by_read]

    # UMI correct length
    masks_raw["UMI_correct_length"] = [len(u) == umi_length for u in umi_by_read]

    # No N bases
    masks_raw["CB_no_N"] = ['N' not in cb for cb in cb_by_read]
    masks_raw["UMI_no_N"] = ['N' not in u for u in umi_by_read]

    # Quality score >= min_qscore
    masks_raw["good_CB_UMI_QC"] = [
        all(ord(c) - 33 >= min_qscore for c in qcs[i])
        for i in range(N)
    ]

    masks = {}
    for key, arr in masks_raw.items():
        masks[key] = [i for i, v in enumerate(arr) if v]

    masks["valid_provenance_structure"] = [
        i for i in range(N)
        if all(masks_raw[k][i] for k in masks_raw)
    ]

    return masks


def parse_indrops_provenance(
    headers: List[str],
    platform: str = "inDrops",
) -> Tuple[List[str], List[str], List[str]]:
    """Parse InDrops cell barcode and UMI from read headers.

    MATLAB equivalent: @SCFastQData/parse_indrops_provenance

    InDrops encodes CB/UMI in read headers.

    Args:
        headers: List of FASTQ read headers.
        platform: 'inDrops' or 'bgi'.

    Returns:
        Tuple of (cb_list, umi_list, qc_list).
    """
    cbs = []
    umis = []
    qcs = []
    for h in headers:
        parts = h.split()
        cb = ""
        umi = ""
        for p in parts:
            if p.startswith("CB:"):
                cb = p.split(":")[-1]
            elif p.startswith("UMI:"):
                umi = p.split(":")[-1]
            elif p.startswith("CR:"):
                cb = p.split(":")[-1]
            elif p.startswith("CY:"):
                qcs.append(p.split(":")[-1])
        cbs.append(cb)
        umis.append(umi)
    return cbs, umis, qcs


def _maybe_gzopen(path: str):
    """Open a file, transparently decompressing .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)
