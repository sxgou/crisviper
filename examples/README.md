# Example Data / 示例数据

## Files

| File | Description |
|------|-------------|
| `reference.fa` | CARLIN reference sequence (332 bp, 10 targets) |
| `EPSC2_L1_1.fq.gz` | Paired-end R1 FASTQ (mouse EPSC2 L1 sample) |
| `EPSC2_L1_2.fq.gz` | Paired-end R2 FASTQ (mouse EPSC2 L1 sample) |
| `test.tsv` | Small TSV query file for quick testing |

## Usage

### Full pipeline test (FASTQ input)

```bash
bash run_lineage_test.sh
```

### Quick test (TSV input)

```bash
crisviper align --reference examples/reference.fa \
  --queries examples/test.tsv \
  --output results/test
```

### Verify installation

```bash
crisviper --version
crisviper --help
```
