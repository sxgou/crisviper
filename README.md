# CrisViper

**Affine-gap alignment for CRISPR lineage tracing.** Align amplicon sequencing reads to a reference, detect indels and SNVs at cut sites, and trace edited cell lineages — across bulk, multi-target (e.g. CARLIN), and single-cell RNA-seq lineage data.

```bash
pip install crisviper
```

## Quick start

```bash
# Standard amplicon: align reads, detect mutations
crisviper align --reference ref.fa --queries reads.tsv --output results.json

# Multi-target lineage mode: structure-aware alignment with correction
crisviper align --reference ref.fa --queries reads.tsv --output results.json --lineage --report html

# scRNA-seq lineage: from 10x FASTQ to mutation report
crisviper convert fastq-to-tsv --fastq reads.fastq.gz --output reads.tsv
crisviper align --reference ref.fa --queries reads.tsv --output results.json --lineage --report html
```

## Features

**Core algorithm.** Gotoh affine-gap (Needleman-Wunsch extension) with semi-global alignment by default. NumPy-vectorized DP delivers ~4.3x per-sequence speedup.

**Lineage mode (`--lineage`).** Position-aware gap penalties: low at cut sites, high in conserved regions. Post-alignment correction pipeline filters false-positive point mutations, converts dense mismatch blocks to indels, realigns repetitive elements between duplicated targets, and merges isolated matches that split continuous deletions.

**DP-native features.** Gap exit bonus, short match discount, dense mismatch penalty, homology penalty, and isolated base penalty — all integrated into the DP recurrence so the aligner makes better decisions during alignment, not after.

**Parallel processing.** ProcessPoolExecutor with automatic core detection (capped at 12 threads). Falls back to single-thread on worker failure. OMP_NUM_THREADS=1 set automatically to avoid NumPy thread conflicts under fork.

**Output.** JSON, TSV, or both. Optional HTML report with mutation classification, length distributions, editing efficiency, per-target editing charts, mutation segment plots, cross-target chord diagrams, and allele heatmaps. Summary tables (TSV) for allele frequency, per-target editing, filter reasons, indel length distributions, and event-level details.

**Allele confidence filter.** Point-mutation-only alleles require `readCount > 5` (`--min-reads-sub`). Indel-containing alleles are unfiltered by default (`--min-reads-indel 0`). Adjust thresholds for higher precision or recall.

**Summary tables.** The pipeline generates 6 TSV tables alongside the main output: `allele_frequency.tsv`, `per_target_editing.tsv`, `filter_reason.tsv`, `deletion_length.tsv`, `insertion_length.tsv`, and `event_level_details.tsv` — covering aggregate statistics, per-target editing rates, filter causes, length distributions, and per-event details with target overlap info.

## Docs

- [Algorithm](docs/ALGORITHM.md) — DP formulation, vectorization, native features, correction pipeline
- [User guide](docs/USER_GUIDE.md) — installation, commands, parameters, examples

## License

MIT
