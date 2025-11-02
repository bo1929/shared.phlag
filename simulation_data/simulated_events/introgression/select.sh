seq 1 50 | xargs -t -I{} bash -c "sort -gk4 candidates/introgression_candidates_{}.tsv | tail -10 > selected/introgression_selected_{}.tsv"
seq 1 50 | xargs -t -I{} bash -c "sort -gk3 candidates/introgression_candidates_{}.tsv | head -10 >> selected/introgression_selected_{}.tsv"
seq 1 50 | xargs -t -I{} bash -c "sort selected/introgression_selected_{}.tsv | uniq > selected/introgression_selected_{}.tsv.tmp && mv selected/introgression_selected_{}.tsv.tmp selected/introgression_selected_{}.tsv"

seq 1 50 | xargs -t -I{} bash -c "head -n25 candidates/introgression_candidates_{}.tsv | shuf | head -5 > selected-difficult/introgression_selected_{}.tsv"
seq 1 50 | xargs -t -I{} bash -c "tail -n+25 candidates/introgression_candidates_{}.tsv | head -n25 | shuf | head -5 >> selected-difficult/introgression_selected_{}.tsv"
