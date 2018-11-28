#!/bin/bash

USAGE="Usage: $0 GENOME_NAME"

if [[ "$#" -lt 1 ]]; then
    echo "$USAGE"
    exit 1
fi

echo "==> parameters <=="
echo "    " $@

GENOME_NAME=$1

cd ~/data/alignment/gene-paralog/${GENOME_NAME}/yml
cp ../data/paralog.yml ../yml

echo "==> paralog_adjacent"
runlist span    --op pad    -n 2000 paralog.yml               -o paralog_adjacent.1.yml
runlist compare --op diff  paralog_adjacent.1.yml paralog.yml -o paralog_adjacent.2.yml
runlist span    --op excise -n 100  paralog_adjacent.2.yml    -o paralog_adjacent.3.yml
runlist span    --op fill   -n 100  paralog_adjacent.3.yml    -o paralog_adjacent.4.yml
runlist genome ../data/chr.sizes -o genome.yml
runlist compare --op intersect paralog_adjacent.4.yml genome.yml -o paralog_adjacent.5.yml
mv paralog_adjacent.5.yml paralog_adjacent.yml
rm paralog_adjacent.*.yml

echo "==> paralog_gene"
mkdir -p ~/data/alignment/gene-paralog/${GENOME_NAME}/yml/gene
runlist split ../feature/all-gene.yml -o ~/data/alignment/gene-paralog/${GENOME_NAME}/yml/gene
runlist compare --op union gene/gene.yml gene/upstream.yml -o paralog_gene.1.yml
runlist compare --op union gene/downstream.yml paralog_gene.1.yml -o paralog_gene.2.yml
runlist compare --op intersect paralog_gene.2.yml paralog.yml -o paralog_gene.3.yml
mv paralog_gene.3.yml paralog_gene.yml
rm paralog_gene.*.yml

for ftr in paralog paralog_adjacent paralog_gene; do
    echo "==> ${ftr} coverages"
    sleep 1;
    runlist stat -s ../data/chr.sizes ../yml/${ftr}.yml -o ../stat/${ftr}.yml.csv

    echo "==> ${ftr} stats"
    cd ~/data/alignment/gene-paralog/${GENOME_NAME}/stat
    cat ../yml/repeat.family.txt \
        | parallel -j 1 --keep-order "
            echo \"==> {}\";
            runlist stat2 ../yml/${ftr}.yml ../yml/{}.yml -s ../data/chr.sizes \
                --op intersect --all \
                -o ../stat/${GENOME_NAME}.${ftr}.{}.csv
        "
    echo "key,chr_length,${ftr}_size,key_length,key_size,c1,c2,ratio" \
        > ../stat/${GENOME_NAME}.${ftr}.all-repeat.csv
    cat ../yml/repeat.family.txt \
        | parallel -j 1 --keep-order "
            cat ../stat/${GENOME_NAME}.${ftr}.{}.csv \
                | perl -nl -e '/^chr_length/ and next; print qq({},\$_)'
        " >> ../stat/${GENOME_NAME}.${ftr}.all-repeat.csv
    cat ../yml/repeat.family.txt \
        | parallel -j 1 --keep-order "rm ../stat/${GENOME_NAME}.${ftr}.{}.csv"
done
