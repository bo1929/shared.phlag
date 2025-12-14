#!/bin/bash
TMPDIR=${4:-$(mktemp -d)}
if [[ "${4}" == "${TMPDIR}" ]]; then
	TMPDIR=$(mktemp -d)
fi
trap 'rm -rf -- "$TMPDIR"' EXIT
NUM_THREADS=${5:-1}

OUTDIR=${3:-${PWD}}
mkdir -p ${OUTDIR}

GENE_TREES=$(realpath ${2})
SPECIES_TREE=$(realpath ${1})

NUM_GENES=3000

date

echo "Estimating CU branch lengths and simulating gene trees..."
# # astral4_coalescent_unit -C -c ${SPECIES_TREE} -i ${GENE_TREES} -o ${TMPDIR}/cu_tree.tree
# # CUSPECIES_TREE="${TMPDIR}/cu_tree.tree"
python simulate_gene_trees.py -i ${SPECIES_TREE} -g ${GENE_TREES} -n ${NUM_GENES} -o ${OUTDIR} -t ${NUM_THREADS}

ROOT_ND=$(nw_labels -r ${OUTDIR}/labelled_cu_tree.tree)
QQS_CMD="astral4 -C -c labelled_cu_tree.tree \
  --root ${ROOT_ND} \
  -u 3 -i <(sed '{}q;d' simulated.gtrees) > /dev/null 2>&1 \
  && gsed 's/^/{}	/' freqQuad.csv && rm freqQuad.csv"

echo "Computing the QQS distribution..."
(cd ${OUTDIR} && seq 1 ${NUM_GENES} | xargs -S1024 -I{} bash -c "${QQS_CMD}" \
  | cut -f1,2,3,5,6 | grep -e $'\tt1\t' -e $'\tt2\t' \
  > nullDist.tsv)

(cd ${OUTDIR} && seq 1 | xargs -S1024 -I{} bash -c "${QQS_CMD}" \
  | cut -f2,3,4 | grep -e $'\tt1\t' -e $'\tt2\t' > nameMap.tsv)

date
