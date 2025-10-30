seq 1 50 | xargs -t -I{} bash -c "head -n25 candidates/introgression_candidates_{}.tsv | shuf | head -5 > selected/introgression_selected_{}.tsv"
seq 1 50 | xargs -t -I{} bash -c "tail -n+25 candidates/introgression_candidates_{}.tsv | head -n25 | shuf | head -5 >> selected/introgression_selected_{}.tsv"
