#!/bin/bash
OUTPUT_DIR=../simulated_events/single_event/population_increase_10X/
mkdir -p ${OUTPUT_DIR}
NUM_THREADS=8

for p in 0.05 0.10 0.15 0.20 0.25 0.30 0.35; do
  for d in 0.75 1.00; do
    p_lbl="${p//./}"
    d_lbl="${d//./}"
    TMPDIR=$(mktemp -d -u -p ${OUTPUT_DIR} p${p_lbl}-d${d_lbl}-XXXXX)
    echo ${TMPDIR}
    seq 1 50 \
    | xargs -I{} -S 1024 -P ${NUM_THREADS} -t bash -c "
    python simulate_mixture_condition.py \
    -x ../gene_trees/SR201_default_condition/{}_1X_S201_0_haploid.gtrees \
    -y ../gene_trees/SR201_10X_population/{}_1X_S201_0_1e6_haploid.gtrees \
    -p ${p} -d ${d}  \
    -o ${TMPDIR}-{}_1X_S201"
  done
done
