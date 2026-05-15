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

**Output.** JSON, TSV, or both. Optional HTML report with mutation classification, length distributions, and editing efficiency.

## Docs

- [Algorithm](docs/ALGORITHM.md) — DP formulation, vectorization, native features, correction pipeline
- [User guide](docs/USER_GUIDE.md) — installation, commands, parameters, examples

## License

MIT
