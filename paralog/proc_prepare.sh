#!/bin/bash

USAGE="Usage: $0 GENOME_NAME BASE_DIR"

if [[ "$#" -lt 1 ]]; then
    echo >&2 "$USAGE"
    exit 1
fi

echo "==> parameters <=="
echo "    " $@

GENOME_NAME=$1
BASE_DIR=${2:-~/data/alignment/Ensembl}

cd ~/data/alignment/gene-paralog/${GENOME_NAME}/data

echo "==> get genome"
if [[ -f ~/data/alignment/gene-paralog/${GENOME_NAME}/data/genome.fa ]]; then
    echo "genome.fa exists"
else
    find ${BASE_DIR}/${GENOME_NAME} -type f -name "*.fa" |
        sort |
        xargs cat |
        perl -nl -e '/^>/ or $_ = uc; print' \
        > genome.fa
fi

echo "==> run RepeatMasker"
if [[ -f ~/data/alignment/gene-paralog/${GENOME_NAME}/data/genome.fa.out ]]; then
    echo "genome.fa.out exists"
else
    RepeatMasker genome.fa -species Viridiplantae -xsmall --parallel 8
    rm genome.fa.cat.gz genome.fa.masked
    rm -fr RM_*
fi

echo "==> run dustmasker"
dustmasker -in genome.fa -infmt fasta -out - -outfmt interval |
    perl -nl -e '
        BEGIN { $name = qq{}}
        if ( /^>(\S+)/ ) {
            $name = $1;
        } elsif ( /(\d+)\s+\-\s+(\d+)/ ) {
            print qq{$name:$1-$2};
        }
    ' \
    > dustmasker.output.txt
runlist cover dustmasker.output.txt -o dustmasker.yml

echo "==> Convert gff3 to runlists"
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/feature

# For Atha, 0m41.121s. With --clean, 50m49.994s
time perl ~/Scripts/withncbi/util/gff2runlist.pl \
    --file ../data/chr.gff \
    --size ../data/chr.sizes \
    --range 2000
