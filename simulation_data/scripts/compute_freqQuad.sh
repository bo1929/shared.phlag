#!/bin/bash
TMPDIR=${4:-$(mktemp -d)}
if [[ "${4}" == "${TMPDIR}" ]]; then
	TMPDIR=$(mktemp -d)
fi
trap 'rm -rf -- "$TMPDIR"' EXIT

OUTDIR=${3:-${PWD}}
mkdir -p ${OUTDIR}

GENE_TREES=$(realpath ${2})
SPECIES_TREE=$(realpath ${1})

NUM_GENES=$(cat ${GENE_TREES} | wc -l)
ROOT_ND=$(nw_labels -r ${SPECIES_TREE})

(cd ${OUTDIR} && seq 1 ${NUM_GENES} | \
  xargs -S1024 -I{} bash -c "astral4 -C -c ${SPECIES_TREE} --root ${ROOT_ND} \
  -u 3 -i <(sed '{}q;d' ${GENE_TREES}) > /dev/null 2>&1 \
  && gsed 's/^/{}	/' freqQuad.csv && rm freqQuad.csv" \
  | cut -f1,2,3,5,6 | grep -e $'\tt1\t' -e $'\tt2\t' > emissionsQQS.tsv)
