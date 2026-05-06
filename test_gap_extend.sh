python3 carlin_tool.py align \
     --threads 16 \
     --gap-extend -1 \
     --lineage \
     --reference example_data/reference.fa \
     --queries example_data/test_queries.tsv \
     --output results/alignment_results_v2 \
     --format all \
     --report html \
     --report-output results/analysis_report_v2 \
     --allele-top-n 50 2>results/log.txt


