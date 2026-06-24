"""Tests for single-cell FASTQ parsing and I/O (crisviper/io.py)."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
from crisviper import (
    parse_10x_provenance, parse_10x_fastq,
    filter_sc_cbs_and_umis, parse_indrops_provenance,
    read_queries_tsv, read_reference_fasta, save_tsv,
)


class TestParse10xProvenance:
    def test_basic(self):
        r1_reads = ["ACGTACGTACGTACGT" + "ACGTACGTACGT"]  # 16bp CB + 12bp UMI
        qc_reads = ["~" * 28]  # high quality
        cbs, umis, qcs = parse_10x_provenance(r1_reads, qc_reads, 16, 12)
        assert cbs == ["ACGTACGTACGTACGT"]
        assert umis == ["ACGTACGTACGT"]
        assert len(qcs) == 1

    def test_empty_lists(self):
        cbs, umis, qcs = parse_10x_provenance([], [], 16, 12)
        assert cbs == []
        assert umis == []
        assert len(qcs) == 0

    def test_qc_length_mismatch_uses_shortest(self):
        """When qc_reads length differs from r1_reads, should handle gracefully."""
        r1 = ["ACGTACGTACGTACGT" + "ACGTACGTACGT"]
        qc = ["~" * 28, "~" * 28]
        cbs, umis, qcs = parse_10x_provenance(r1, qc, 16, 12)
        assert len(cbs) == len(qcs)

    def test_short_read(self):
        r1_reads = ["ACGT"]  # shorter than cb_length
        qc_reads = ["~~~"]
        cbs, umis, qcs = parse_10x_provenance(r1_reads, qc_reads, 16, 12)
        assert cbs == ["ACGT"]
        assert umis == [""]

    def test_multiple(self):
        r1_reads = ["A" * 28, "C" * 28]
        qc_reads = ["~" * 28, "~" * 28]
        cbs, umis, qcs = parse_10x_provenance(r1_reads, qc_reads, 16, 12)
        assert cbs == ["A" * 16, "C" * 16]
        assert umis == ["A" * 12, "C" * 12]


class TestFilterScCbsAndUmis:
    def test_all_good(self):
        cbs = ["ACGTACGTACGTACGT"]
        umis = ["ACGTACGTACGT"]
        read_cb = [0, 0, 0]
        read_umi = [0, 0, 0]
        qcs = ["~" * 28, "~" * 28, "~" * 28]
        masks = filter_sc_cbs_and_umis(cbs, read_cb, umis, read_umi, qcs)
        assert len(masks["valid_provenance_structure"]) == 3

    def test_bad_qc_filtered(self):
        cbs = ["ACGTACGTACGTACGT"]
        umis = ["ACGTACGTACGT"]
        read_cb = [0, 0]
        read_umi = [0, 0]
        qcs = ["~" * 28, "!" * 28]  # second read has very low quality
        masks = filter_sc_cbs_and_umis(cbs, read_cb, umis, read_umi, qcs, min_qscore=20)
        assert len(masks["good_CB_UMI_QC"]) == 1

    def test_n_in_cb_filtered(self):
        cbs = ["ACGTACGTACGTACGN"]
        umis = ["ACGTACGTACGT"]
        read_cb = [0, 0]
        read_umi = [0, 0]
        qcs = ["~" * 28, "~" * 28]
        masks = filter_sc_cbs_and_umis(cbs, read_cb, umis, read_umi, qcs)
        assert len(masks["CB_no_N"]) == 0

    def test_wrong_length_filtered(self):
        cbs = ["ACGT"]  # too short
        umis = ["ACGTACGTACGT"]
        read_cb = [0]
        read_umi = [0]
        qcs = ["~" * 16]
        masks = filter_sc_cbs_and_umis(cbs, read_cb, umis, read_umi, qcs, cb_length=16)
        assert len(masks["CB_correct_length"]) == 0


class TestParseIndropsProvenance:
    def test_basic(self):
        headers = ["@read1 CB:ACGTACGTACGTACGT UMI:ACGTACGTACGT"]
        cbs, umis, qcs = parse_indrops_provenance(headers)
        assert cbs == ["ACGTACGTACGTACGT"]
        assert umis == ["ACGTACGTACGT"]

    def test_empty_header(self):
        headers = ["@read1"]
        cbs, umis, qcs = parse_indrops_provenance(headers)
        assert cbs == [""]
        assert umis == [""]


class TestParse10xFastq:
    def test_missing_files(self):
        with pytest.raises((FileNotFoundError, OSError)):
            parse_10x_fastq("/nonexistent/R1.fastq", "/nonexistent/R2.fastq")


# ═══════════════════════════════════════════════════════════════
# read_queries_tsv tests
# ═══════════════════════════════════════════════════════════════

class TestReadQueriesTsv:
    def test_basic(self, tmp_path):
        tsv = tmp_path / "queries.tsv"
        tsv.write_text("readName\tcellBC\tUMI\treadCount\tseq\n"
                       "r1\tCB1\tUMI1\t10\tACGT\n"
                       "r2\tCB2\tUMI2\t5\tTGCA\n")
        rows = read_queries_tsv(str(tsv))
        assert len(rows) == 2
        assert rows[0]["readName"] == "r1"
        assert rows[0]["readCount"] == 10
        assert rows[0]["seq"] == "ACGT"

    def test_readcount_converted_to_int(self, tmp_path):
        tsv = tmp_path / "q.tsv"
        tsv.write_text("readName\tcellBC\tUMI\treadCount\tseq\nr1\tCB\tU\t10\tA\n")
        rows = read_queries_tsv(str(tsv))
        assert isinstance(rows[0]["readCount"], int)

    def test_empty_file_exits(self, tmp_path):
        tsv = tmp_path / "empty.tsv"
        tsv.write_text("readName\tcellBC\tUMI\treadCount\tseq\n")
        with pytest.raises(SystemExit):
            read_queries_tsv(str(tsv))

    def test_whitespace_in_fields_stripped(self, tmp_path):
        """Whitespace around readName, cellBC, UMI should be stripped."""
        tsv = tmp_path / "ws.tsv"
        tsv.write_text("readName\tcellBC\tUMI\treadCount\tseq\n"
                       "  r1  \t CB1 \t UMI1 \t10\tACGT\n")
        rows = read_queries_tsv(str(tsv))
        assert rows[0]["readName"] == "r1"
        assert rows[0]["cellBC"] == "CB1"
        assert rows[0]["UMI"] == "UMI1"

    def test_file_not_found(self):
        with pytest.raises(SystemExit):
            read_queries_tsv("/nonexistent/path.tsv")


# ═══════════════════════════════════════════════════════════════
# read_reference_fasta tests
# ═══════════════════════════════════════════════════════════════

class TestReadReferenceFasta:
    def test_basic(self, tmp_path):
        fa = tmp_path / "ref.fa"
        fa.write_text(">ref_name\nACGTACGT\n")
        seq = read_reference_fasta(str(fa))
        assert seq == "ACGTACGT"

    def test_multiline(self, tmp_path):
        fa = tmp_path / "ref.fa"
        fa.write_text(">ref\nACGT\nACGT\n")
        seq = read_reference_fasta(str(fa))
        assert seq == "ACGTACGT"

    def test_file_not_found(self):
        with pytest.raises(SystemExit):
            read_reference_fasta("/nonexistent.fa")

    def test_empty_fasta_file(self, tmp_path):
        fa = tmp_path / "empty.fa"
        fa.write_text("")
        with pytest.raises(SystemExit):
            read_reference_fasta(str(fa))

    def test_multi_sequence_fasta_ignores_second(self, tmp_path):
        """Only first sequence should be returned."""
        fa = tmp_path / "multi.fa"
        fa.write_text(">seq1\nACGT\n>seq2\nTGCA\n")
        seq = read_reference_fasta(str(fa))
        assert seq == "ACGT"


# ═══════════════════════════════════════════════════════════════
# save_tsv tests
# ═══════════════════════════════════════════════════════════════

class TestSaveTsv:
    def test_basic(self, tmp_path):
        out = tmp_path / "out.tsv"
        rows = [{"readName": "r1", "cellBC": "CB", "UMI": "U",
                 "readCount": 5, "seq": "ACGT"}]
        save_tsv(rows, str(out))
        content = out.read_text()
        assert "readName" in content
        assert "r1" in content
        assert "ACGT" in content

    def test_empty_rows(self, tmp_path):
        out = tmp_path / "empty.tsv"
        save_tsv([], str(out))
        content = out.read_text()
        assert "readName" in content  # header only
        assert len(content.strip().split("\n")) == 1
