#!/bin/bash

USAGE="Usage: $0 GENOME_NAME"

if [[ "$#" -lt 1 ]]; then
    echo "$USAGE"
    exit 1
fi

echo "==> parameters <=="
echo "    " $@

GENOME_NAME=$1
cd ~/data/alignment/gene-paralog/${GENOME_NAME}/repeat
find . -type f -name "*.yml" -or -name "*.csv" | parallel rm

echo "==> gff types"
gzip -dcf ../data/chr.gff |
    perl -nla -e '/^#/ and next; print $F[2]' |
    sort |
    uniq -c \
    > gff.type.txt

gzip -dcf ../data/chr.gff |
    perl -nla -e '
        /^#/ and next;
        $F[2] eq q{repeat_region} or next;
        $F[8] =~ /description\=(\w+)/i or next;
        print qq{$F[2]\t$F[1]\t$1};
        ' |
    sort |
    uniq -c \
    > gff.repeat.txt

echo "==> rmout families"
cat ../data/genome.fa.out |
    perl -nla -e '/^\s*\d+/ or next; print $F[10]' |
    sort |
    uniq -c \
    > rmout.family.txt
cp ../data/genome.fa.out ../repeat
cp ../data/genome.fa.tbl ../repeat

echo "==> rmout results"
perl ~/Scripts/withncbi/util/rmout2runlist.pl \
    --file ../data/genome.fa.out \
    --size ../data/chr.sizes

runlist split ../feature/all-repeat.yml -o .
runlist merge *.yml -o all-repeat.1.yml
find . -type f -name "*.yml" -not -name "all-*" | parallel rm

runlist span --op excise -n 10 --mk all-repeat.1.yml -o all-repeat.2.yml   # remove small spans
runlist span --op fill   -n 10 --mk all-repeat.2.yml -o all-repeat.3.yml   # fill small holes
runlist split all-repeat.3.yml -o .
mv all-repeat.3.yml all-repeat.yml
rm all-repeat.*.yml

echo "==> basic stat"
runlist stat all-repeat.yml -s ../data/chr.sizes --mk --all

echo "==> find repeat families large enough"
cat all-repeat.yml.csv |
    perl -nla -F"," -e '
        next if $F[3] =~ /^c/;
        print $F[0] if $F[3] > 0.0005
    ' \
    > repeat.family.txt
cat repeat.family.txt |
    parallel -j 8 "
        cp {}.yml ../yml
    "

echo "==> clean"
mv all-repeat.yml.csv ../stat
cp repeat.family.txt ../yml
