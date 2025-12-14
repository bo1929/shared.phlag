#!/bin/bash
NUM_THREADS=16

# OUTPUT_DIR=../simulated_events/single_event/recombination_suppression_fixed/
# mkdir -p ${OUTPUT_DIR}
# 
# for p in 0.05 0.10 0.15 0.20 0.25 0.30 0.35; do
#   for r in 0.80 0.90 0.99; do
#     p_lbl="${p//./}"
#     r_lbl="${r//./}"
#     TMPDIR=$(mktemp -d -u -p ${OUTPUT_DIR} p${p_lbl}-r${r_lbl}-XXXXX)
#     echo ${TMPDIR}
#     seq 1 50 \
#     | xargs -I{} -P ${NUM_THREADS} -t bash -c "
#     python simulate_recombination_supression.py \
#     -g ../gene_trees/SR201_default_condition/{}_1X_S201_0_haploid.gtrees \
#     -p ${p} -r ${r} \
#     -o ${TMPDIR}-{}_1X_S201"
#   done
# done
# 
# OUTPUT_DIR=../simulated_events/single_event/recombination_suppression_random/
# mkdir -p ${OUTPUT_DIR}
# 
# for p in 0.05 0.10 0.15 0.20 0.25 0.30 0.35; do
#   for r in 0.80 0.90 0.99; do
#     p_lbl="${p//./}"
#     r_lbl="${r//./}"
#     TMPDIR=$(mktemp -d -u -p ${OUTPUT_DIR} p${p_lbl}-r${r_lbl}-XXXXX)
#     echo ${TMPDIR}
#     seq 1 50 \
#     | xargs -I{} -P ${NUM_THREADS} -t bash -c "
#     python simulate_recombination_supression.py \
#     -g ../gene_trees/SR201_default_condition/{}_1X_S201_0_haploid.gtrees \
#     -p ${p} -r ${r} \
#     -o ${TMPDIR}-{}_1X_S201 \
#     --option random"
#   done
# done
# 
# OUTPUT_DIR=../simulated_events/single_event/recombination_suppression_support-population_10X
# mkdir -p ${OUTPUT_DIR}
# 
# for p in 0.05 0.10 0.15 0.20 0.25 0.30 0.35; do
#   for r in 0.33 0.66 0.99; do
#     p_lbl="${p//./}"
#     r_lbl="${r//./}"
#     TMPDIR=$(mktemp -d -u -p ${OUTPUT_DIR} p${p_lbl}-r${r_lbl}-XXXXX)
#     echo ${TMPDIR}
#     seq 1 50 \
#     | xargs -I{} -P ${NUM_THREADS} -t bash -c "
#     python simulate_recombination_supression.py \
#     -g ../gene_trees/SR201_10X_population/{}_1X_S201_0_1e6_haploid.gtrees \
#     -p ${p} -r ${r} \
#     -o ${TMPDIR}-{}_1X_S201 \
#     --option support"
#   done
# done

OUTPUT_DIR=../simulated_events/single_event/recombination_suppression_support-default/
mkdir -p ${OUTPUT_DIR}

for p in 0.05 0.10 0.15 0.20 0.25 0.30; do
  for r in 0.33 0.66 0.99; do
    p_lbl="${p//./}"
    r_lbl="${r//./}"
    TMPDIR=$(mktemp -d -u -p ${OUTPUT_DIR} p${p_lbl}-r${r_lbl}-XXXXX)
    echo ${TMPDIR}
    seq 1 50 \
    | xargs -I{} -P ${NUM_THREADS} -t bash -c "
    python simulate_recombination_supression.py \
    -g ../gene_trees/SR201_default_condition/{}_1X_S201_0_haploid.gtrees \
    -p ${p} -r ${r} \
    -o ${TMPDIR}-{}_1X_S201 \
    --option support"
  done
done
