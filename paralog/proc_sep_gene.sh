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

cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat

echo "==> intersect"
# parallel in 1 threads to save memory
for ftr in gene upstream downstream exon CDS intron five_prime_UTR three_prime_UTR; do
    echo ${ftr}
done \
    | parallel -j 1 --keep-order "
        echo \"==> {} ${FEATURE_BASE}\";
        sleep 1;
        jrunlist statop \
            ../data/chr.sizes ../feature/sep-{}.yml ${FEATURE_FILE}  \
            --op intersect --all \
            -o stat.sep-{}.${FEATURE_BASE}.csv;
        cat stat.sep-{}.${FEATURE_BASE}.csv \
            | cut -d ',' -f 1,3,5 \
            > stat.sep-{}.${FEATURE_BASE}.csv.tmp;
    "

echo "==> concat gene"
printf "gene_id," > ${GENOME_NAME}.gene.${FEATURE_BASE}.csv
for ftr in gene upstream downstream; do
    printf "${ftr}_length,${ftr}_${FEATURE_BASE},"
done >> ${GENOME_NAME}.gene.${FEATURE_BASE}.csv
echo >> ${GENOME_NAME}.gene.${FEATURE_BASE}.csv

for ftr in gene upstream downstream; do
    cat stat.sep-${ftr}.${FEATURE_BASE}.csv.tmp
done \
    | grep -v "^key" \
    | perl ~/Scripts/withncbi/util/merge_csv.pl --concat -f 0 -o stdout \
    >> ${GENOME_NAME}.gene.${FEATURE_BASE}.csv

echo "==> concat trans"
printf "trans_id," > ${GENOME_NAME}.trans.${FEATURE_BASE}.csv
for ftr in exon CDS intron five_prime_UTR three_prime_UTR; do
    printf "${ftr}_length,${ftr}_${FEATURE_BASE},"
done >> ${GENOME_NAME}.trans.${FEATURE_BASE}.csv
echo >> ${GENOME_NAME}.trans.${FEATURE_BASE}.csv

for ftr in exon CDS intron five_prime_UTR three_prime_UTR; do
    cat stat.sep-${ftr}.${FEATURE_BASE}.csv.tmp
done \
    | grep -v "^key" \
    | perl ~/Scripts/withncbi/util/merge_csv.pl --concat -f 0 -o stdout \
    >> ${GENOME_NAME}.trans.${FEATURE_BASE}.csv

echo "==> clean"
rm stat.sep-*.${FEATURE_BASE}.csv.tmp
rm stat.sep-*.${FEATURE_BASE}.csv
