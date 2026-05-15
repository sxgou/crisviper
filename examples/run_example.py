"""CrisViper example workflow.

Demonstrates the complete analysis pipeline:
  1. Convert FASTQ to TSV
  2. Standard (single-target) alignment
  3. Lineage (multi-target) alignment with DP-native features
  4. Generate HTML report

Usage:
  python examples/run_example.py                          # full data
  python examples/run_example.py --subset 1000            # first 1000 reads only
"""

import argparse
import os
import sys
import time

# Ensure the package is importable (works when run from repo root)
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from crisviper import (
    PipelineConfig, Pipeline, QueryRecord,
    read_reference_fasta, fastq_to_dataframe,
    get_amplicon_structure,
    setup_logging, get_logger,
)

log = get_logger("example")

EXAMPLE_DIR = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(EXAMPLE_DIR, "reference.fa")
FASTQ = os.path.join(EXAMPLE_DIR, "test.fastq.gz")
OUTPUT_DIR = os.path.join(EXAMPLE_DIR, "output")


def main():
    setup_logging(verbose=True)
    parser = argparse.ArgumentParser(description="CrisViper example workflow")
    parser.add_argument("--subset", type=int, default=0,
                        help="Only process first N reads (default: all)")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Step 1: Load reference ──────────────────────────────────
    log.info("Loading reference: %s", REFERENCE)
    ref_seq = read_reference_fasta(REFERENCE)
    log.info("Reference length: %d bp", len(ref_seq))

    # ── Step 2: Convert FASTQ → TSV ──────────────────────────
    log.info("Converting FASTQ: %s", FASTQ)
    records = fastq_to_dataframe(FASTQ, sample_name="example")

    if args.subset:
        records = records[:args.subset]
    records = [QueryRecord(**r) for r in records]
    log.info("Reads loaded: %d", len(records))

    # ── Step 3: Standard alignment ───────────────────────────
    log.info("=" * 50)
    log.info("Standard (single-target) alignment")
    log.info("=" * 50)

    config_standard = PipelineConfig(threads=4)
    pipeline_standard = Pipeline(config_standard, ref_seq)

    t0 = time.time()
    result_standard = pipeline_standard.run(records)
    t1 = time.time()

    success = result_standard.get_successful()
    mutated = result_standard.get_mutated()
    log.info("Processed %d / %d successful  |  %.1f reads/sec",
             len(success), len(records), len(success) / (t1 - t0) if (t1 - t0) > 0 else 0)
    log.info("Mutated: %d / %d (%.1f%%)",
             len(mutated), len(success),
             len(mutated) / len(success) * 100 if success else 0)

    # ── Step 4: Lineage mode alignment ───────────────────────
    log.info("=" * 50)
    log.info("Lineage (multi-target) alignment")
    log.info("=" * 50)

    cutsites = get_amplicon_structure(ref_seq)
    if cutsites:
        log.info("Detected %d cut sites", len(cutsites))
        for cs in cutsites:
            log.info("  %s: %d-%d", cs.name, cs.start, cs.end)

    config_lineage = PipelineConfig(
        lineage_mode=True,
        gap_exit_bonus=-1.0,
        short_match_window=3,
        short_match_discount=0.5,
        dense_mismatch_penalty=-2.0,
        homology_penalty=-1.0,
        isolated_base_penalty=-2.0,
        threads=4,
    )
    pipeline_lineage = Pipeline(config_lineage, ref_seq)

    t0 = time.time()
    result_lineage = pipeline_lineage.run(records)
    t1 = time.time()

    success_l = result_lineage.get_successful()
    mutated_l = result_lineage.get_mutated()
    log.info("Processed %d / %d successful  |  %.1f reads/sec",
             len(success_l), len(records), len(success_l) / (t1 - t0) if (t1 - t0) > 0 else 0)
    log.info("Mutated: %d / %d (%.1f%%)",
             len(mutated_l), len(success_l),
             len(mutated_l) / len(success_l) * 100 if success_l else 0)

    n_corrected = sum(
        1 for r in result_lineage.results if r.success and r.stats
        and getattr(r.stats, "n_mutations_corrected", 0) > 0
    )
    log.info("Corrected alignments: %d", n_corrected)

    # ── Step 5: Save results and report ──────────────────────
    output_json = os.path.join(OUTPUT_DIR, "alignments.json")
    output_html = os.path.join(OUTPUT_DIR, "report.html")

    from crisviper import save_alignment_results, generate_report
    save_alignment_results(
        [r.to_dict() for r in result_lineage.results],
        output_json, "json",
    )
    generate_report(
        [r.to_dict() for r in result_lineage.results],
        output_html, "html",
        ref_length=len(ref_seq),
        ref_seq=ref_seq,
        cutsites=cutsites,
    )

    log.info("Results saved to: %s", output_json)
    log.info("Report saved to:  %s", output_html)
    log.info("Done.")


if __name__ == "__main__":
    main()
