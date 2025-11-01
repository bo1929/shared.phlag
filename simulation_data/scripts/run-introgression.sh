#!/bin/bash
NUM_THREADS=32

run_introgression_cases () {
  OUTPUT_DIR=../simulated_events/single_event/introgression-separate_blocks/
  mkdir -p ${OUTPUT_DIR}
  for p in 0.05 0.10 0.15 0.20 0.25 0.30; do
    for b in 1; do
      p_lbl="${p//./}"
      b_lbl="${b//./}"
      TMPDIR=$(mktemp -d -u -p ${OUTPUT_DIR} p${p_lbl}-b${b_lbl}-XXXXX)
      echo ${TMPDIR}
      exec < ../simulated_events/introgression/selected/introgression_selected_$1.tsv || exit 1a
      while IFS='	' read -r recipient donor v1 v2; do
        python simulate_introgression.py \
          -s ../species_trees/estimated/SR201_default_condition/$1_1X_S201_0_haploid.caster-pair \
          -g ../gene_trees/SR201_default_condition/$1_1X_S201_0_haploid.gtrees \
          -p ${p} -b ${b} -d ${donor} -r ${recipient} \
          -o ${TMPDIR}b${b_lbl}-d${donor}_r${recipient}-$1_1X_S201
      done
    done
  done
}

export -f run_introgression_cases

seq 1 50 | xargs -I{} -P ${NUM_THREADS} -t bash -c "run_introgression_cases {}"
