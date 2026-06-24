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
    """Read a FASTQ file and return a {sequence: occurrence_count} dictionary."""
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


def _read_fastq_groups(fastq_path: str) -> Dict[str, List[str]]:
    """Read a FASTQ file and return {sequence: [read_name, ...]} preserving original read names.

    Note: For large FASTQ files (millions of reads), this keeps all read names in memory.
    Consider memory constraints when enabling keep_read_names with high-depth datasets.
    """
    groups = {}
    if fastq_path.endswith('.gz'):
        with gzip.open(fastq_path, 'rt') as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = str(record.seq)
                groups.setdefault(seq, []).append(record.id)
    else:
        for record in SeqIO.parse(fastq_path, "fastq"):
            seq = str(record.seq)
            groups.setdefault(seq, []).append(record.id)
    return groups


def fastq_to_dataframe(fastq_path: str, sample_name: str = "sample",
                        keep_read_names: bool = False) -> List[Dict]:
    """
    Convert a FASTQ file to a list of deduplicated sequence dicts.

    Each unique sequence is represented as a single entry with its
    readCount set to the number of times it appeared in the FASTQ.

    When keep_read_names=True, each row also contains an 'original_read_names'
    field listing the original FASTQ read identifiers for that sequence.

    Args:
        fastq_path: Path to FASTQ file (supports .gz compression).
        sample_name: Sample name for read labeling.
        keep_read_names: Whether to preserve original read names.
            Note: enabling this keeps all read names in memory, which may
            consume significant memory on large datasets.

    Returns:
        List of dicts with keys: readName, cellBC, UMI, readCount, seq.
        Also includes 'original_read_names' when keep_read_names=True.
    """
    if keep_read_names:
        groups = _read_fastq_groups(fastq_path)
        rows = []
        for i, (seq, read_names) in enumerate(groups.items()):
            rows.append({
                "readName": f"{sample_name}_seq{i+1}",
                "cellBC": sample_name,
                "UMI": f"UMI{i+1}",
                "readCount": len(read_names),
                "seq": seq,
                "original_read_names": read_names,
            })
        total_reads = sum(len(n) for n in groups.values())
        log.info("Read %d unique sequences (%d total reads) from %s (names preserved)",
                 len(rows), total_reads, fastq_path)
        return rows

    counts = _read_fastq_counts(fastq_path)

    # Build list of dicts
    rows = []
    for i, (seq, count) in enumerate(counts.items()):
        rows.append({
            "readName": f"{sample_name}_seq{i+1}",
            "cellBC": sample_name,          # Virtual cell barcode (sample-level)
            "UMI": f"UMI{i+1}",             # Virtual UMI (sequence-level)
            "readCount": count,             # Observation count for this sequence
            "seq": seq
        })

    log.info("Read %d unique sequences (%d total reads) from %s", len(rows), sum(counts.values()), fastq_path)
    return rows


def fastq_to_fasta(fastq_path: str, output_fasta: str, sample_name: str = "sample") -> None:
    """
    Convert a FASTQ file to FASTA format with metadata in headers.

    Each FASTA header includes: readName, cellBC, UMI, readCount.

    Args:
        fastq_path: Input FASTQ file path (supports .gz).
        output_fasta: Output FASTA file path.
        sample_name: Sample name for read labeling.
    """
    counts = _read_fastq_counts(fastq_path)

    with open(output_fasta, 'w') as f:
        for i, (seq, count) in enumerate(counts.items()):
            read_name = f"{sample_name}_seq{i+1}"
            f.write(f">{read_name} cellBC={sample_name} UMI=UMI{i+1} readCount={count}\n{seq}\n")

    total_reads = sum(counts.values())
    log.info("Wrote %d unique sequences (%d reads) to %s", len(counts), total_reads, output_fasta)


def fastq_to_fasta_from_rows(rows: List[Dict], output_path: str) -> None:
    """Write a list of row dicts (from merge_paired_end/fastq_to_dataframe) to FASTA format.

    Args:
        rows: List of dicts with keys: readName, cellBC, UMI, readCount, seq.
        output_path: Output FASTA file path.
    """
    with open(output_path, 'w') as f:
        for row in rows:
            f.write(f">{row['readName']} cellBC={row['cellBC']} "
                    f"UMI={row['UMI']} readCount={row['readCount']}\n"
                    f"{row['seq']}\n")
    log.info("Data saved to %s", output_path)


def _check_fastp() -> None:
    """Check that fastp is installed, exit with error if not."""
    if shutil.which("fastp") is None:
        log.error("fastp is not installed. Paired-end merging requires fastp. "
                  "Install with: brew install fastp  # macOS, or "
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
    keep_read_names: bool = False,
) -> List[Dict]:
    """Merge paired-end reads using fastp and return deduplicated row dicts.

    When keep_read_names=True, each row also contains an 'original_read_names'
    field listing the original FASTQ read identifiers for that sequence.

    Args:
        fastq1: R1 FASTQ file path.
        fastq2: R2 FASTQ file path.
        min_overlap: Minimum overlap length (fastp: overlap_len_require).
        max_mismatch_rate: Max mismatch rate %% in overlap (fastp: overlap_diff_percent_limit).
        max_mismatch_diff: Max absolute mismatches in overlap (fastp: overlap_diff_limit).
        require_qual: Minimum base quality threshold (fastp: -q).
        sample_name: Sample name for read labeling.
        keep_read_names: Whether to preserve original read names.

    Returns:
        List of dicts (same format as fastq_to_dataframe).
        Includes 'original_read_names' when keep_read_names=True.
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
            "-j", os.devnull,
            "-h", os.devnull,
        ]

        log.info("Running fastp to merge paired-end reads...")
        log.debug("  cmd: %s", " ".join(cmd))

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            if result.stdout:
                for line in result.stdout.splitlines():
                    log.debug("  fastp: %s", line)
        except subprocess.CalledProcessError as e:
            log.error("fastp failed (exit %d)", e.returncode)
            if e.stderr:
                for line in e.stderr.splitlines():
                    log.error("  %s", line)
            sys.exit(1)

        # Read merged results
        if not os.path.exists(merged_out) or os.path.getsize(merged_out) == 0:
            log.warning("fastp produced no merged reads!")
            return []

        if keep_read_names:
            groups = {}
            with open(merged_out) as f:
                for record in SeqIO.parse(f, "fastq"):
                    seq = str(record.seq)
                    groups.setdefault(seq, []).append(record.id)
            if not groups:
                log.warning("fastp produced no merged reads!")
                return []
            rows = []
            for i, (seq, read_names) in enumerate(groups.items()):
                rows.append({
                    "readName": f"{sample_name}_seq{i+1}",
                    "cellBC": sample_name,
                    "UMI": f"UMI{i+1}",
                    "readCount": len(read_names),
                    "seq": seq,
                    "original_read_names": read_names,
                })
            total = sum(len(n) for n in groups.values())
            log.info("fastp merge complete: %d unique sequences (%d total reads, names preserved)",
                     len(rows), total)
            return rows

        counts = {}
        with open(merged_out) as f:
            for record in SeqIO.parse(f, "fastq"):
                seq = str(record.seq)
                counts[seq] = counts.get(seq, 0) + 1

        if not counts:
            log.warning("fastp produced no merged reads!")
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

        log.info("fastp merge complete: %d unique sequences (%d total reads)",
                 len(rows), sum(counts.values()))
        return rows


def save_tsv(rows: List[Dict], output_path: str) -> None:
    """Save a list of row dicts to a TSV file."""
    if not rows:
        log.warning("No rows to save, writing empty TSV to %s", output_path)
        fieldnames = ["readName", "cellBC", "UMI", "readCount", "seq"]
    else:
        fieldnames = list(rows[0].keys())
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    log.info("Data saved to %s", output_path)


def read_reference_fasta(fasta_path: str) -> str:
    """Read a reference sequence from a FASTA file and return the upper-case sequence string."""
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            return str(record.seq).upper()
    except Exception as e:
        log.error("Failed to read reference file - %s", e)
        sys.exit(1)
    log.error("Reference file is empty: %s", fasta_path)
    sys.exit(1)


def read_queries_tsv(tsv_path: str) -> List[Dict]:
    """
    Read query sequences from a TSV file.

    Expected format: readName, cellBC, UMI, readCount, seq
    """
    rows = []
    try:
        with open(tsv_path, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for line_num, row in enumerate(reader, start=2):
                if 'seq' not in row or not row['seq'].strip():
                    log.warning("Line %d missing sequence field, skipping", line_num)
                    continue
                row['seq'] = row['seq'].upper()
                for key in ('readName', 'cellBC', 'UMI'):
                    if key in row:
                        row[key] = row[key].strip()
                row['readCount'] = int(row.get('readCount', 1))
                rows.append(row)
    except Exception as e:
        log.error("Failed to read TSV file - %s", e)
        sys.exit(1)

    if not rows:
        log.error("TSV file is empty or contains no valid sequences: %s", tsv_path)
        sys.exit(1)

    log.info("Read %d query sequences from %s", len(rows), tsv_path)
    return rows


def read_queries_fasta(fasta_path: str) -> List[Dict]:
    """
    Read query sequences from a FASTA file.

    FASTA headers can contain metadata:
    >readName cellBC=sample1 UMI=UMI1 readCount=2
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
        log.error("Failed to read FASTA file - %s", e)
        sys.exit(1)

    if not rows:
        log.error("FASTA file is empty or invalid: %s", fasta_path)
        sys.exit(1)

    log.info("Read %d query sequences from %s", len(rows), fasta_path)
    return rows


# ═══════════════════════════════════════════════════════════════
# Single-cell FASTQ parsing (10x / InDrops / BGI)
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

    if len(r1_records) != len(r2_records):
        raise ValueError(
            f"R1 ({len(r1_records)}) and R2 ({len(r2_records)}) have different read counts"
        )

    N = len(r1_records)
    raw_cb = [str(r.seq) for r in r1_records]
    # Convert Phred quality scores (list of ints) to proper ASCII quality strings
    # r.letter_annotations["phred_quality"] returns e.g. [40, 30, 20, ...]
    # We encode as ASCII: chr(score + 33) per the Sanger/Illumina 1.8+ FASTQ standard.
    # Note: Phred+64 (Illumina 1.3/1.5) is NOT supported; encountering such files
    # will produce garbage quality strings.
    raw_qc = ["".join(chr(q + 33) for q in r.letter_annotations["phred_quality"])
              for r in r1_records]
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
        len(qcs[i]) > 0 and all(ord(c) - 33 >= min_qscore for c in qcs[i])
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
        qc = ""
        for p in parts:
            if p.startswith("CB:"):
                cb = p.split(":")[-1]
            elif p.startswith("UMI:"):
                umi = p.split(":")[-1]
            elif p.startswith("CR:"):
                cb = p.split(":")[-1]
            elif p.startswith("CY:"):
                qc = p.split(":")[-1]
        cbs.append(cb)
        umis.append(umi)
        qcs.append(qc)
    return cbs, umis, qcs


def _maybe_gzopen(path: str):
    """Open a file, transparently decompressing .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)
