lineage-tracer analyze \
--fastq_r1 test_data/EPSC2_L1_1.fq.gz \
--fastq_r2 test_data/EPSC2_L1_2.fq.gz \
--config test_data/EPSC2_config.json \
--output_dir ./results \
--sample_name EPSC2 \
--kmer_size 6 \
--min_overlap 10 \
--threads 10 \
--no_trim_primers
