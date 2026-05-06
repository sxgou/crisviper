python3 carlin_tool.py align \
     --threads 16 \
     --global \
     --lineage \
     --reference example_data/reference.fa \
     --queries example_data/test_queries.tsv \
     --output results/alignment_results \
     --format all \
     --report html \
     --report-output results/analysis_report \
     --allele-top-n 50 2>results/log.txt


