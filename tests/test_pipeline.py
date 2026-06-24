"""End-to-end integration test for the full crisviper pipeline.

Tests the complete flow:
  1. Creates synthetic CARLIN-like reference and query sequences
  2. Converts to TSV format
  3. Runs the alignment pipeline
  4. Runs with standard and lineage modes
  5. Generates reports (JSON and HTML)
  6. Validates output statistics
"""

import sys
import os
import json
import tempfile
import shutil
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pytest

from crisviper import (
    Pipeline, PipelineConfig, QueryRecord,
    align_single,
    generate_report,
    get_amplicon_structure,
    read_queries_tsv,
    MutationType,
)


from shared import CARLIN_REF


def _make_single_mutated_query(ref_seq: str, pos: int, new_base: str) -> str:
    """Create a point mutation at the specified position"""
    lst = list(ref_seq)
    lst[pos] = new_base
    return ''.join(lst)


def _make_deletion_query(ref_seq: str, start: int, length: int) -> str:
    """Create a deletion at the specified position"""
    return ref_seq[:start] + ref_seq[start + length:]


def _make_insertion_query(ref_seq: str, pos: int, insert: str) -> str:
    """Create an insertion at the specified position"""
    return ref_seq[:pos] + insert + ref_seq[pos:]


def _build_test_queries(ref_seq: str) -> list:
    """Build a list of test sequences containing various mutations"""
    queries = [
        # Wildtype (unmutated)
        QueryRecord(readName="wt_high",  cellBC="test", UMI="UMI1", readCount=50, seq=ref_seq),
        QueryRecord(readName="wt_low",   cellBC="test", UMI="UMI2", readCount=1,  seq=ref_seq),

        # Point mutation inside cutsite window
        QueryRecord(readName="sub_in_window", cellBC="test", UMI="UMI3", readCount=15,
                    seq=_make_single_mutated_query(ref_seq, 43, 'T')),  # Inside Target1 cutsite

        # Point mutation outside window (should be corrected/filtered)
        QueryRecord(readName="sub_out_window", cellBC="test", UMI="UMI4", readCount=12,
                    seq=_make_single_mutated_query(ref_seq, 30, 'T')),  # Conserved region

        # Small deletion in cutsite
        QueryRecord(readName="del_cutsite_3bp", cellBC="test", UMI="UMI5", readCount=8,
                    seq=_make_deletion_query(ref_seq, 43, 3)),

        # Large deletion spanning a Target
        QueryRecord(readName="del_large_20bp", cellBC="test", UMI="UMI6", readCount=5,
                    seq=_make_deletion_query(ref_seq, 41, 20)),

        # Insertion in conserved region
        QueryRecord(readName="ins_small", cellBC="test", UMI="UMI7", readCount=4,
                    seq=_make_insertion_query(ref_seq, 50, "AAA")),

        # Primer3 anchoring failure (modified Primer3 region, last 33 bases)
        QueryRecord(readName="primer3_bad", cellBC="test", UMI="UMI8", readCount=3,
                    seq=ref_seq[:299] + "N" * 33),  # Completely replace Primer3 region

        # Full-length deletion
        QueryRecord(readName="del_full_target", cellBC="test", UMI="UMI9", readCount=6,
                    seq=_make_deletion_query(ref_seq, 41, 7)),  # Delete Target1 cutsite

        # Point mutation + low readCount (use 3bp mutation to avoid alignment algorithm gap preference)
        QueryRecord(readName="sub_low_rc", cellBC="test", UMI="UMI10", readCount=2,
                    seq=_make_single_mutated_query(
                        _make_single_mutated_query(
                            _make_single_mutated_query(ref_seq, 100, 'G'), 101, 'A'), 102, 'T')),

        # Complex mutation (deletion + insertion)
        QueryRecord(readName="complex", cellBC="test", UMI="UMI11", readCount=7,
                    seq=_make_deletion_query(
                        _make_insertion_query(ref_seq, 50, "ACG"), 44, 5
                    )),
    ]
    return queries


# ═══════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════

class TestEndToEndPipeline:
    """End-to-end pipeline test"""

    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup before each test"""
        self.ref_seq = CARLIN_REF
        self.queries = _build_test_queries(self.ref_seq)
        self.temp_dir = tempfile.mkdtemp()
        yield
        shutil.rmtree(self.temp_dir)

    def _run_pipeline(self, lineage_mode=False):
        """Run the pipeline and return the result"""
        config = PipelineConfig(
            lineage_mode=lineage_mode,
            primer5_len=23,
            primer3_len=33,
            primer5_threshold=19,
            primer3_threshold=29,
            min_reads_sub=10,
            min_reads_indel=3,
            threads=1,
        )
        pipeline = Pipeline(config=config, ref_seq=self.ref_seq)
        if lineage_mode:
            pipeline.load_cutsites()
        return pipeline.run(self.queries)

    def test_pipeline_runs_successfully(self):
        """Pipeline runs correctly and returns results"""
        result = self._run_pipeline()
        assert result is not None
        assert len(result.results) == len(self.queries)
        assert result.stats.total_queries == len(self.queries)

    def test_primer_anchoring_detects_failures(self):
        """Primer3 anchoring failure should be correctly identified"""
        result = self._run_pipeline()
        failed = result.get_failed()
        failed_names = [r.query.readName for r in failed]
        assert "primer3_bad" in failed_names

    def test_wildtype_high_rc_passes(self):
        """Wildtype (high readCount) should pass"""
        result = self._run_pipeline()
        successful = result.get_successful()
        names = [r.query.readName for r in successful]
        assert "wt_high" in names

    def test_wildtype_low_rc_passes(self):
        """Wildtype (low readCount but no mutation) should pass"""
        result = self._run_pipeline()
        successful = result.get_successful()
        names = [r.query.readName for r in successful]
        assert "wt_low" in names

    def test_snv_in_window_passes(self):
        """Point mutation inside cutsite window should pass"""
        result = self._run_pipeline()
        successful = result.get_successful()
        names = [r.query.readName for r in successful]
        assert "sub_in_window" in names

    def test_snv_low_rc_filtered(self):
        """Point mutation + low readCount should be filtered"""
        result = self._run_pipeline()
        successful = result.get_successful()
        successful_names = [r.query.readName for r in successful]
        assert "sub_low_rc" not in successful_names
        # Cross-validate: should appear in failed results
        failed = result.get_failed()
        failed_names = [r.query.readName for r in failed]
        assert "sub_low_rc" in failed_names, \
            f"sub_low_rc not in failed list either: {failed_names}"

    def test_mutations_detected_in_results(self):
        """Successfully aligned results should contain mutation events"""
        result = self._run_pipeline()
        for r in result.get_successful():
            # All successful results should have a mutations field
            assert hasattr(r, 'mutations')
            if r.query.readName == "sub_in_window":
                assert len(r.mutations) > 0
                assert r.mutations[0].type == MutationType.SUBSTITUTION

    def test_stats_collected(self):
        """Pipeline statistics should be reasonable"""
        result = self._run_pipeline()
        stats = result.stats
        assert stats.total_queries > 0
        assert stats.successful <= stats.total_queries
        assert stats.failed >= 0
        assert stats.successful + stats.failed == stats.total_queries

    def test_lineage_mode_produces_results(self):
        """Lineage tracing alignment mode also works correctly"""
        result = self._run_pipeline(lineage_mode=True)
        assert result is not None
        assert len(result.results) > 0

    def test_lineage_mode_detects_cutsites(self):
        """Lineage tracing mode automatically detects cutsites"""
        # Directly test get_amplicon_structure
        cutsites = get_amplicon_structure(self.ref_seq)
        assert cutsites is not None
        assert len(cutsites) == 10  # CARLIN has 10 Targets

    def test_pipeline_stats_consistency(self):
        """Statistics should be consistent (mutated + unmutated = successful)"""
        result = self._run_pipeline()
        stats = result.stats
        assert stats.mutated_sequences + stats.unmutated_sequences == stats.successful
        # Cross-validate: stats.failed must match actual failed count
        actual_failed = sum(1 for r in result.results if not r.success)
        assert stats.failed == actual_failed, \
            f"stats.failed={stats.failed} but {actual_failed} results have success=False"

    def test_pipeline_result_conversion_to_dict(self):
        """AlignmentResult.to_dict() is compatible with the old format"""
        result = self._run_pipeline()
        for r in result.get_successful()[:3]:
            d = r.to_dict()
            assert "readName" in d
            assert "score" in d
            assert "aligned_ref" in d
            assert "aligned_query" in d
            assert "stats" in d

    def test_error_result_conversion(self):
        """to_dict() of error results should include an error field"""
        result = self._run_pipeline()
        for r in result.get_failed():
            d = r.to_dict()
            assert d["error"] != ""
            assert d["score"] is None

    def test_save_and_read_tsv_roundtrip(self):
        """TSV write and read-back should be consistent"""
        import tempfile
        save_path = os.path.join(self.temp_dir, "test_queries.tsv")

        # Write TSV first
        from crisviper import save_tsv
        dict_rows = [{"readName": q.readName, "cellBC": q.cellBC,
                      "UMI": q.UMI, "readCount": q.readCount, "seq": q.seq}
                     for q in self.queries]
        save_tsv(dict_rows, save_path)

        # Read it back
        read_back = read_queries_tsv(save_path)
        assert len(read_back) == len(self.queries)
        assert read_back[0]["readName"] == self.queries[0].readName

    def test_report_json_generation(self):
        """JSON report generation should succeed"""
        result = self._run_pipeline()
        output_results = [r.to_dict() for r in result.results]
        report_path = os.path.join(self.temp_dir, "test_report")

        generate_report(output_results, report_path, fmt="json",
                         ref_length=len(self.ref_seq))

        # Verify JSON file exists and contains required fields
        json_path = report_path + ".json"
        assert os.path.exists(json_path)
        with open(json_path) as f:
            report = json.load(f)
        assert "summary" in report
        assert "editing_efficiency_pct" in report["summary"]

    def test_html_report_generation(self):
        """HTML report generation should succeed (skip matplotlib check)"""
        result = self._run_pipeline()
        output_results = [r.to_dict() for r in result.results]
        report_path = os.path.join(self.temp_dir, "test_html_report")

        generate_report(output_results, report_path, fmt="html",
                         ref_length=len(self.ref_seq))

        html_path = report_path + ".html"
        assert os.path.exists(html_path)
        # Verify HTML contains key content
        with open(html_path) as f:
            content = f.read()
        assert "CARLIN" in content
        assert "Editing Efficiency" in content
        # Plots may be empty (no matplotlib), but HTML structure should be complete
        assert "<html" in content
        assert "modal" in content

    def test_pipeline_with_no_cutsites(self):
        """Standard mode (without dual cutsite) also works"""
        config = PipelineConfig(
            lineage_mode=False,
            primer5_len=23,
            primer3_len=33,
            min_reads_sub=10,
            min_reads_indel=3,
            threads=1,
        )
        pipeline = Pipeline(config=config, ref_seq=self.ref_seq)
        result = pipeline.run(self.queries)
        assert len(result.results) > 0

    def test_align_single_function(self):
        """align_single function can be called directly"""
        from crisviper import align_single, PipelineConfig
        config = PipelineConfig()
        result = align_single(self.queries[0], self.ref_seq, config)
        assert result is not None
        assert result.success, f"Wildtype alignment should succeed, got: {result.error}"
        assert len(result.mutations) == 0, f"Wildtype should have 0 mutations, got {len(result.mutations)}"

    def test_mutations_in_del_query(self):
        """Deletion mutation alignment results should contain MutationEvent"""
        result = self._run_pipeline()
        del_results = [r for r in result.get_successful()
                       if "del_" in r.query.readName]
        for r in del_results:
            if r.mutations:
                has_deletion = any(m.type == MutationType.DELETION
                                   for m in r.mutations)
                if not has_deletion:
                    # May have been merged into a complex event
                    has_indel = any(m.type == MutationType.INDEL
                                      for m in r.mutations)
                    assert has_indel or not r.stats.has_indel

    def test_parallel_vs_sequential_consistency(self):
        """threads=1 serial vs threads=2 parallel produce identical stats and per-query results."""
        # Run with threads=1
        result_serial = self._run_pipeline()
        # Run the same queries with threads=2
        config2 = PipelineConfig(
            lineage_mode=False,
            primer5_len=23, primer3_len=33,
            primer5_threshold=19, primer3_threshold=29,
            min_reads_sub=10, min_reads_indel=3,
            threads=2,
        )
        pipeline2 = Pipeline(config=config2, ref_seq=self.ref_seq)
        result_parallel = pipeline2.run(self.queries)

        # Compare pipeline-level statistics
        s = result_serial.stats
        p = result_parallel.stats
        assert p.total_queries == s.total_queries
        assert p.successful == s.successful
        assert p.failed == s.failed
        assert p.total_reads == s.total_reads
        assert p.mutated_sequences == s.mutated_sequences
        assert p.unmutated_sequences == s.unmutated_sequences
        assert p.mutated_reads == s.mutated_reads
        assert p.n_anchor_failed == s.n_anchor_failed
        assert p.n_noise_filtered == s.n_noise_filtered

        # Compare per-query results (sorted by readName for determinism)
        serial_by_name = {r.query.readName: r for r in result_serial.results}
        parallel_by_name = {r.query.readName: r for r in result_parallel.results}
        assert set(serial_by_name.keys()) == set(parallel_by_name.keys())
        for name in serial_by_name:
            sr = serial_by_name[name]
            pr = parallel_by_name[name]
            assert pr.success == sr.success, f"{name}: success mismatch"
            assert len(pr.mutations) == len(sr.mutations), f"{name}: mutation count mismatch"
            if not sr.success:
                assert pr.error == sr.error, f"{name}: error mismatch"


class TestBoundaryInputs:
    """Boundary and edge-case input handling."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.ref_seq = CARLIN_REF
        self.cfg = PipelineConfig(
            primer5_len=23, primer3_len=33,
            primer5_threshold=19, primer3_threshold=29,
            min_reads_sub=10, min_reads_indel=3,
        )

    def test_empty_query_sequence(self):
        """Empty query should fail with alignment error."""
        q = QueryRecord(readName="empty", cellBC="test", UMI="U1", readCount=1, seq="")
        result = align_single(q, self.ref_seq, self.cfg)
        assert not result.success
        assert "alignment" in result.error.lower() or "failed" in result.error.lower()

    def test_all_n_query(self):
        """All-N query should still align (valid characters, just ambiguous)."""
        q = QueryRecord(readName="alln", cellBC="test", UMI="U1", readCount=50, seq="N" * len(self.ref_seq))
        result = align_single(q, self.ref_seq, self.cfg)
        # Should either succeed (many mismatches) or fail primer anchoring
        if result.success:
            assert len(result.mutations) > 0 or not result.stats.has_indel
        else:
            assert result.error != ""

    def test_query_shorter_than_primers(self):
        """Query shorter than primer5 should fail primer anchoring."""
        q = QueryRecord(readName="short", cellBC="test", UMI="U1", readCount=1, seq="ACGT")
        result = align_single(q, self.ref_seq, self.cfg)
        assert not result.success

    def test_reference_shorter_than_primers(self, tmp_path):
        """Pipeline should handle a very short reference gracefully."""
        cfg = PipelineConfig(primer5_len=50, primer3_len=50, min_reads_sub=10, min_reads_indel=3)
        q = QueryRecord(readName="q", cellBC="test", UMI="U1", readCount=1, seq="A" * 10)
        result = align_single(q, "ACGT", cfg)
        assert not result.success
