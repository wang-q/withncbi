#!/bin/bash

USAGE="Usage: $0 GENOME_NAME FEATURE_FILE"

if [[ "$#" -lt 2 ]]; then
    echo "$USAGE"
    exit 1
fi

echo "==> parameters <=="
echo "    " $@

GENOME_NAME=$1
FEATURE_FILE=$2
FEATURE_BASE=`basename "${FEATURE_FILE%.*}"`

echo "==> all-gene"
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

runlist stat2 ../feature/all-gene.yml ${FEATURE_FILE} -s ../data/chr.sizes \
    --op intersect --mk --all \
    -o ${GENOME_NAME}.all-gene.${FEATURE_BASE}.csv

