# LineageTracer-Amplicon: Paired-End Sequencing Usage Examples

## Overview

This document provides practical examples for using LineageTracer-Amplicon with paired-end sequencing data.

## Example Data Structure

Assume you have the following files:
```
data/
├── sample_R1.fastq.gz    # Forward reads
├── sample_R2.fastq.gz    # Reverse reads
└── amplicon_structure.json  # Configuration file
```

## Basic Usage Examples

### 1. Single-End Analysis (R1 only)

```bash
lineage-tracer analyze \
  --fastq_r1 data/sample_R1.fastq.gz \
  --config data/amplicon_structure.json \
  --output_dir results_single_end \
  --sample_name MySample \
  --threads 4
```

### 2. Paired-End Analysis with Default Merging

```bash
lineage-tracer analyze \
  --fastq_r1 data/sample_R1.fastq.gz \
  --fastq_r2 data/sample_R2.fastq.gz \
  --config data/amplicon_structure.json \
  --output_dir results_paired_end \
  --sample_name MySample \
  --threads 4
```

**Default merging parameters:**
- `--min_overlap 30`: Require at least 30bp overlap
- `--max_mismatch_rate 0.2`: Allow up to 20% mismatches in overlap region
- `--min_quality 20`: Minimum base quality score

### 3. Paired-End Analysis with Custom Overlap Requirements

```bash
lineage-tracer analyze \
  --fastq_r1 data/sample_R1.fastq.gz \
  --fastq_r2 data/sample_R2.fastq.gz \
  --config data/amplicon_structure.json \
  --output_dir results_custom \
  --sample_name MySample \
  --min_overlap 20 \
  --max_mismatch_rate 0.15 \
  --min_quality 30 \
  --threads 8
```

### 4. Paired-End Analysis Without Merging (Analyze R1 only)

```bash
lineage-tracer analyze \
  --fastq_r1 data/sample_R1.fastq.gz \
  --fastq_r2 data/sample_R2.fastq.gz \
  --config data/amplicon_structure.json \
  --output_dir results_no_merge \
  --sample_name MySample \
  --skip_merge \
  --threads 4
```

## Output Files Structure

### For Paired-End Analysis with Merging:

```
results_paired_end/
├── intermediate/
│   ├── MySample_merged.fastq.gz      # Merged reads
│   ├── MySample_fastp.json           # JSON report with statistics
│   └── MySample_fastp.html           # HTML report (visualization)
├── per_read_annotation.tsv           # Mutation states per read
├── barcode_frequencies.tsv           # Barcode frequency table
└── statistics.json                   # Summary statistics
```

### For Single-End or No-Merge Analysis:

```
results_single_end/
├── per_read_annotation.tsv
├── barcode_frequencies.tsv
└── statistics.json
```

## Fastp Merge Statistics

The `intermediate/MySample_fastp.json` file contains detailed merging statistics:

```json
{
  "summary": {
    "before_filtering": {
      "total_reads": 100000,
      "total_bases": 15000000
    },
    "after_filtering": {
      "total_reads": 95000,
      "total_bases": 14250000
    }
  },
  "merging": {
    "merged_reads": 85000,
    "merge_rate": 0.8947
  }
}
```

## Troubleshooting

### 1. Fastp Not Installed

**Error:**
```
Error: fastp is required for paired-end merging but not found.
```

**Solution:**
```bash
# Install fastp using conda
conda install -c bioconda fastp

# Or install from source
# See: https://github.com/OpenGene/fastp
```

### 2. Insufficient Overlap

**Symptoms:** Low merge rate in fastp report

**Solutions:**
- Decrease `--min_overlap` value (e.g., `--min_overlap 15`)
- Increase `--max_mismatch_rate` (e.g., `--max_mismatch_rate 0.3`)
- Check if your library preparation generates sufficient overlap

### 3. Low Quality Reads

**Symptoms:** Many reads filtered out

**Solutions:**
- Adjust `--min_quality` parameter
- Check raw data quality with FastQC
- Consider trimming adapters before analysis

## Best Practices

1. **Always check fastp reports**: Review `fastp.html` and `fastp.json` to ensure merging worked correctly
2. **Test different overlap parameters**: If merge rate is low, try adjusting `--min_overlap` and `--max_mismatch_rate`
3. **Compare with single-end**: Run both merged and unmerged analyses to compare results
4. **Monitor resource usage**: Use `--threads` appropriately for your system

## Performance Tips

- Use `--threads` to parallelize both fastp merging and analysis
- For large datasets, consider using `--max_reads` for initial testing
- Intermediate files can be deleted after analysis to save space
- Use `--verbose` flag for detailed logging during troubleshooting