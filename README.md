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

Gotoh affine-gap with semi-global alignment. NumPy-vectorized DP. Position-aware gap penalties in lineage mode. DP-native features (gap exit bonus, short match discount, dense mismatch penalty, homology penalty, isolated base penalty). Parallel ProcessPoolExecutor. JSON/TSV/HTML output with allele heatmaps.

## Docs
- [Algorithm](docs/ALGORITHM.md)
- [User guide](docs/USER_GUIDE.md)

## License
MIT
