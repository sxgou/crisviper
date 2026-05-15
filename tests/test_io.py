"""Tests for single-cell FASTQ parsing (crisviper/io.py)."""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
from crisviper import (
    parse_10x_provenance, parse_10x_fastq,
    filter_sc_cbs_and_umis, parse_indrops_provenance,
)


class TestParse10xProvenance:
    def test_basic(self):
        r1_reads = ["ACGTACGTACGTACGT" + "ACGTACGTACGT"]  # 16bp CB + 12bp UMI
        qc_reads = ["~" * 28]  # high quality
        cbs, umis, qcs = parse_10x_provenance(r1_reads, qc_reads, 16, 12)
        assert cbs == ["ACGTACGTACGTACGT"]
        assert umis == ["ACGTACGTACGT"]
        assert len(qcs) == 1

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
