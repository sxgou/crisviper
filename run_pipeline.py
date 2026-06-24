"""Run CrisViper pipeline on example data, output to results/"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Use available CPUs, cap at 12 for reasonable parallelism
CPU_COUNT = min(12, os.cpu_count() or 1)

from crisviper import (
    PipelineConfig, Pipeline, QueryRecord,
    read_reference_fasta, fastq_to_dataframe, get_amplicon_structure,
    save_alignment_results, save_summary_tables, generate_report,
    setup_logging, get_logger,
)

log = get_logger("run")


def main():
    setup_logging(verbose=True)

    EXAMPLE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "examples")
    RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Load reference
    ref_seq = read_reference_fasta(os.path.join(EXAMPLE_DIR, "reference.fa"))
    log.info("Reference length: %d bp", len(ref_seq))

    # Load FASTQ
    records = fastq_to_dataframe(os.path.join(EXAMPLE_DIR, "test.fastq.gz"), sample_name="example")
    records = [QueryRecord(**r) for r in records]
    log.info("Reads loaded: %d", len(records))

    # Detect cutsites
    cutsites = get_amplicon_structure(ref_seq)
    if cutsites:
        log.info("Detected %d cut sites", len(cutsites))
        for cs in cutsites:
            log.info("  %s: %d-%d", cs.name, cs.start, cs.end)

    # ── Lineage mode ──
    log.info("=" * 60)
    log.info("Lineage (multi-target) alignment")
    log.info("=" * 60)

    config_l = PipelineConfig(
        lineage_mode=True,
        gap_exit_strength=-1.0,
        short_match_window=3, short_match_discount=0.5,
        dense_mismatch_penalty=-2.0,
        homology_penalty=-1.0,
        isolated_base_penalty=-2.0,
        threads=CPU_COUNT,
    )
    pipeline_l = Pipeline(config_l, ref_seq)
    t0 = time.time()
    result_l = pipeline_l.run(records)
    t1 = time.time()

    success_l = result_l.get_successful()
    mutated_l = result_l.get_mutated()
    log.info("Lineage: %d/%d success, %d mutated (%.1f%%) | %.1f reads/sec",
             len(success_l), len(records), len(mutated_l),
             len(mutated_l)/len(success_l)*100 if success_l else 0,
             len(success_l)/(t1-t0) if (t1-t0) > 0 else 0)

    # Save results
    output_json = os.path.join(RESULTS_DIR, "alignments.json")
    output_html = os.path.join(RESULTS_DIR, "report.html")

    result_dicts = [r.to_dict() for r in result_l.results]
    save_alignment_results(result_dicts, output_json, "json")
    generate_report(
        result_dicts, output_html, "html",
        ref_length=len(ref_seq), ref_seq=ref_seq, cutsites=cutsites,
    )
    save_summary_tables(result_dicts, RESULTS_DIR, ref_seq=ref_seq, cutsites=cutsites)

    log.info("Results saved to %s/", RESULTS_DIR)
    log.info("Done.")


if __name__ == "__main__":
    main()
