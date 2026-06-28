# CrisViper

**Affine-gap alignment for CRISPR lineage tracing.**

Align amplicon sequencing reads to a reference, detect indels and SNVs at cut sites, and trace edited cell lineages across bulk, multi-target (e.g. CARLIN), and single-cell RNA-seq lineage data.

Python 3.8+ | MIT License

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

Many additional options are available: `--config` (YAML config), `--format {json,tsv,all}`,
`--threads N`, `--cutsites`, `--fastq1` / `--fastq2` (paired-end input), `--report {json,html}`,
`--read-to-allele`, and 20+ alignment parameters. Run `crisviper align --help` for the full list.

## Installation

From source:

```bash
git clone https://github.com/sxgou/crisviper
cd crisviper
python -m venv venv
source venv/bin/activate
pip install -e .
```

Docker:

```bash
docker build -t crisviper .
docker run --rm -v $(pwd)/data:/data crisviper align --reference /data/ref.fa --queries /data/q.tsv --output /data/out.json
```

Dependencies: numpy>=1.19.0, biopython>=1.79, matplotlib>=3.3.0, jinja2>=3.0, numba>=0.55.0, pyyaml>=5.1

## Features

### Alignment engine
- **Global affine-gap DP** — Gotoh-style full-length alignment with match, Ix (insertion in query), Iy (insertion in reference) states. Numba JIT accelerated (28x speedup), with automatic pure-NumPy fallback when Numba is unavailable.
- **Position-aware gap penalties** — Gradient-based profile centered on cut sites for lineage-tracer data (lineage mode). DP-native features: gap exit bonus, short match discount, dense mismatch penalty, homology penalty, isolated base penalty.
- **Lineage-tracer mode** — Structure-aware alignment with gradient-scaled gap penalties, homology penalty profiles, and amplicon-structure auto-detection for multi-target designs (e.g. CARLIN 10-target amplicon).

### Input &amp; output
- **Input formats**: FASTA reference, TSV query table, single-end FASTQ, paired-end FASTQ (R1+R2 with automatic overlap-based merging). UMI/barcode extraction for scRNA-seq (10x Genomics, InDrops).
- **Output formats**: JSON (default), TSV, or both (`--format all`). HTML analysis report with allele heatmaps and mutation summary.
- **Summary tables** — 6 automatically-generated TSV tables: allele frequency, per-target editing rates, filter reason counts, indel length distributions (deletion + insertion), and event-level details.
- **Read-to-allele mapping** — Optional per-read allele assignment output (TSV) for single-cell UMI-level analysis.

### Post-processing pipeline
- **Primer quality check** — Validates 5'/3' primer anchoring for each alignment.
- **Background substitution correction** — Removes systematic substitution noise.
- **DEL→INS→DEL merge** — Collapses interrupted indel events into unified mutation records.
- **Mutation extraction &amp; classification** — SNVs, insertions, deletions classified per cut site.
- **Allele confidence filtering** — Threshold-based per-allele QC (min read counts for substitutions vs. indels, primer match thresholds).

### Downstream analysis modules
- **Denoiser** — UMI/cell barcode denoising via directional-adjacency top-down clustering (Hamming distance 1).
- **Allele caller** — Coarse-grain and exact consensus allele calling from aligned reads.
- **Metrics** — Diversity and heterogeneity metrics: Shannon effective alleles, diversity index, alleles per cell, singleton rates.
- **Threshold** — Statistical read-count threshold computation for UMI/cell barcode filtering (MATLAB-compatible heuristic).

### Parallelism
- `ProcessPoolExecutor`-based parallel batch alignment with configurable threads and chunk size.
- Automatic `OMP_NUM_THREADS=1` / `MKL_NUM_THREADS=1` to prevent thread contention with BLAS.

### YAML configuration

All pipeline parameters, amplicon structure, and cut-site coordinates can be specified via a YAML config file (`--config`). CLI arguments take precedence over YAML values. See [`crisviper_config.yaml`](crisviper_config.yaml) for all available options.

## Docs
- [Algorithm](docs/ALGORITHM.md)
- [User guide](docs/USER_GUIDE.md)
- [Changelog](CHANGELOG.md)
- [Contributing](CONTRIBUTING.md)

## License
MIT
