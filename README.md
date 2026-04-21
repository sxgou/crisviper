# LineageTracer-Amplicon

CRISPR lineage tracing amplicon analysis software for structured CRISPR arrays.

## Overview

This software analyzes amplicon sequencing data from CRISPR-Cas9 lineage tracing experiments with synthetic target arrays. It accurately identifies editing events including large deletions, insertions, and point mutations while filtering out sequencing artifacts.

## Key Features

- **Paired-end read merging**: Support for paired-end sequencing data using fastp with configurable minimum overlap
- **Primer validation**: Only reads containing both Primer5 and Primer3 (with configurable mismatch tolerance) are analyzed
- **Anchor-guided alignment**: Uses conserved regions as anchors for accurate alignment even with large deletions
- **Fragment classification**: Correctly classifies "fragment" sequences within deletions as MMEJ or NHEJ repair products
- **Position-filtered SNV detection**: Only point mutations within cutsite ±3bp windows are considered valid editing events
- **Target state determination**: Generates mutation state matrix for each target (WT, DEL, INS, SNV, COMPLEX)
- **Enrichment analysis**: Validates cutting specificity by comparing breakpoint distribution to expected cutsites

## Installation

### From source

```bash
git clone <repository-url>
cd lineage_tracer_amplicon
pip install -e .
```

### Dependencies

- Python >= 3.9
- biopython >= 1.81
- pandas >= 2.0
- numpy >= 1.24
- scipy >= 1.10
- matplotlib >= 3.7
- pyyaml >= 6.0
- edlib >= 1.3.9
- tqdm >= 4.65
- pysam >= 0.21

**For paired-end merging**:
- [fastp](https://github.com/OpenGene/fastp) >= 0.23.0 (install via `conda install -c bioconda fastp`)

## Quick Start

### 1. Prepare configuration file

Define your amplicon structure in JSON format:

```json
{
  "reference": "ATCG...full sequence...",
  "features": [
    {"name": "Primer5", "start": 0, "end": 23, "type": "primer"},
    {"name": "prefix", "start": 23, "end": 28, "type": "prefix"},
    {"name": "Target1", "start": 28, "end": 48, "type": "target",
     "conserved": [28, 41], "cutsite": [41, 48]},
    ...
  ]
}
```

### 2. Run analysis

**Single-end analysis**:
```bash
lineage-tracer analyze \
  --fastq_r1 sample_R1.fastq.gz \
  --config amplicon_structure.json \
  --output_dir ./results \
  --sample_name Embryo_Day5 \
  --threads 8
```

**Paired-end analysis with merging**:
```bash
lineage-tracer analyze \
  --fastq_r1 sample_R1.fastq.gz \
  --fastq_r2 sample_R2.fastq.gz \
  --config amplicon_structure.json \
  --output_dir ./results \
  --sample_name Embryo_Day5 \
  --min_overlap 30 \
  --max_mismatch_rate 0.2 \
  --threads 8
```

**Paired-end analysis without merging** (use R1 only):
```bash
lineage-tracer analyze \
  --fastq_r1 sample_R1.fastq.gz \
  --fastq_r2 sample_R2.fastq.gz \
  --config amplicon_structure.json \
  --output_dir ./results \
  --skip_merge
```

### 3. View results

- `per_read_annotation.tsv`: Mutation states for each read
- `barcode_frequencies.tsv`: Frequency of each barcode pattern
- `statistics.json`: Summary statistics and editing efficiency
- `intermediate/` directory (for paired-end): Contains fastp merge reports and merged reads

## Configuration Wizard

Generate configuration interactively:

```bash
lineage-tracer config-wizard --output my_config.json
```

## Paired-End Sequencing Support

The software supports paired-end sequencing data through integration with [fastp](https://github.com/OpenGene/fastp). When `--fastq_r2` is provided, reads are automatically merged before analysis.

### Key Parameters for Merging

- `--min_overlap`: Minimum overlap length required for merging (default: 30)
- `--max_mismatch_rate`: Maximum mismatch rate allowed in overlap region (default: 0.2)
- `--skip_merge`: Skip merging even for paired-end data (analyze R1 only)
- `--min_quality`: Minimum base quality for filtering (default: 20)

### Merging Process

1. **Quality filtering**: Remove low-quality bases and reads
2. **Adapter trimming**: Automatically detect and remove adapters
3. **Overlap detection**: Find overlapping regions between R1 and R2
4. **Merging**: Combine overlapping reads into single consensus sequences
5. **Error correction**: Correct sequencing errors in overlap regions

### Output Files for Paired-End Analysis

- `intermediate/{sample}_merged.fastq.gz`: Merged reads
- `intermediate/{sample}_fastp.json`: JSON report with merge statistics
- `intermediate/{sample}_fastp.html`: HTML report for visualization

## Detailed Documentation

### Amplicon Structure

The software supports structured CRISPR arrays with:
- Primer5 (5' primer, 23bp)
- Prefix (5bp)
- 10 targets, each with:
  - Conserved region (13bp)
  - Cutsite region (7bp, contains PAM upstream 3bp)
- PAM_Linker between targets (3bp PAM + 4bp linker)
- Postfix (8bp)
- Primer3 (3' primer, 33bp)

### Algorithm Overview

1. **Read preprocessing**: For paired-end data, merge reads using fastp
2. **Primer validation**: Detect Primer5 at 5' end and Primer3 at 3' end with configurable mismatch tolerance
3. **Anchor matching**: Use k-mer indexing to find conserved anchor sequences in reads
4. **Deletion inference**: Identify gaps between matched anchors as candidate deletions
5. **Fragment classification**: Analyze sequences within deletions as pure deletions, MMEJ, or NHEJ with insertion
6. **SNV filtering**: Only retain point mutations within cutsite ±3bp windows
7. **Target state determination**: Map events to each target and assign mutation state
8. **Barcode generation**: Create barcode string from target states for lineage tracing

### Output Files

- **per_read_annotation.tsv**: Tab-separated file with read ID, validity, target states, barcode, and events JSON
- **barcode_frequencies.tsv**: Frequency table of unique barcodes
- **statistics.json**: JSON summary with editing efficiency, per-target statistics, and enrichment analysis
- **breakpoint_density.bedgraph**: Breakpoint density for visualization in IGV
- **indel_length_distribution.tsv**: Distribution of indel lengths
- **snv_distribution.tsv**: Distribution of valid point mutations

## Testing

Run the test suite:

```bash
cd lineage_tracer_amplicon
python test/test_suite.py
```

## Performance

- Designed for high-throughput analysis (millions of reads)
- Multi-threading support for both fastp merging and analysis
- Memory-efficient processing with streaming
- C++ core for performance-critical operations (planned)

## Citation

If you use this software in your research, please cite:

[Citation information will be added]

## License

[License information will be added]

## Support

For issues and feature requests, please open an issue on GitHub.

## Acknowledgements

This software was developed based on the detailed specification document "CRISPR谱系追踪扩增子分析软件 - 完整开发文档.md".