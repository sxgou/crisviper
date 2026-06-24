#!/usr/bin/env bash
# Use available CPUs, cap at 12
CPU_COUNT=$(python3 -c "import os; print(min(12, os.cpu_count() or 1))")

python -m crisviper.cli align \
--reference examples/reference.fa \
--fastq1 examples/EPSC2_L1_1.fq.gz \
--fastq2 examples/EPSC2_L1_2.fq.gz \
--min-overlap 10 \
--output results/test \
--read-to-allele \
--lineage \
--format all \
--report html \
--threads "$CPU_COUNT" \
--correct-bg-sub \
--allele-top-n 50 \
--min-reads-sub 5 \
--min-reads-indel 5 \
--max-scale 6 \
--gap-exit-strength -1.0 \
--short-match-window 3 \
--short-match-discount 0.5 \
--dense-mismatch-penalty -2.0 \
--homology-window 4 \
--homology-penalty -1.0 \
--isolated-base-penalty -2.0 \
--primer5-threshold 17 \
--primer3-threshold 25
